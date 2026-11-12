#!/usr/bin/env python3
"""Resumable, multi-host sharded driver for rpart v2 (whole-plane box search).

  run   --n N --m M --shards S --select SEL --threads T --output DIR [--binary ./rpart] [--nice 10]
  merge DIR [DIR ...] --into OUT

The c range [1, M] is cut into S consecutive shards. SEL is 'r1,r2/m' (shards j with
j % m in {r1, r2}) or 'a-b' (shards a..b) or 'all'. Each finished shard leaves
shard_XXXX.json (the engine's final line) and shard_XXXX.done.json; a restart with the
same --output skips finished shards, so a reclaim loses only the shards that were running.
progress_<host>.json follows the format read by the dashboard probes.
"""
import argparse, hashlib, json, os, pathlib, socket, subprocess, sys, time

SOURCES = ['src/main.rs', 'Cargo.toml', 'run_rpart.py']
COUNTERS = ['c_values', 'partners', 'pairs', 'q_candidates', 'd_tests', 'three', 'solutions', 'divisor_nodes']


def jst(t=None):
    return time.strftime('%Y-%m-%d %H:%M:%S JST', time.gmtime((t or time.time()) + 9 * 3600))


def sha256(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def write_json_atomic(path, obj):
    tmp = pathlib.Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=1) + '\n')
    os.replace(tmp, path)


def intervals(m, shards):
    cuts = [1 + (m * i) // shards for i in range(shards)] + [m + 1]
    return [(cuts[i], cuts[i + 1] - 1) for i in range(shards)]


def select(sel, shards):
    if sel == 'all':
        return list(range(shards))
    if '/' in sel:
        rs, mod = sel.split('/')
        rs = {int(x) for x in rs.split(',')}
        return [j for j in range(shards) if j % int(mod) in rs]
    a, b = map(int, sel.split('-'))
    return list(range(a, b + 1))


def cmd_run(a):
    root = pathlib.Path(__file__).resolve().parent
    out = pathlib.Path(a.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    binary = pathlib.Path(a.binary).resolve()
    plan = {'engine': 'rpart_v2', 'N': a.n, 'M': a.m, 'shards': a.shards, 'lower': 1, 'bound': a.n,
            'intervals': intervals(a.m, a.shards),
            'source_sha256': {s: sha256(root / s) for s in SOURCES}}
    plan = json.loads(json.dumps(plan))
    pf = out / 'plan.json'
    if pf.exists():
        old = json.loads(pf.read_text())
        if old != plan:
            raise SystemExit('plan.json differs from this invocation (different N, M, shards or sources)')
    else:
        write_json_atomic(pf, plan)
    host = socket.gethostname().split('.')[0]
    chosen = select(a.select, a.shards)
    pending = [j for j in chosen if not (out / f'shard_{j:04d}.done.json').exists()]
    started = time.time()
    prog = {'host': host, 'pid': os.getpid(), 'select': a.select, 'selected': len(chosen), 'shards_total': a.shards,
            'N': a.n, 'M': a.m, 'workers': 1, 'threads_per_worker': a.threads, 'started': jst(started),
            'started_epoch': started, 'skipped_already_done': len(chosen) - len(pending),
            'done': len(chosen) - len(pending), 'failed': [], 'running': {}, 'state': 'running',
            'binary_sha256': sha256(binary)}

    def save():
        prog['updated'] = jst(); prog['updated_epoch'] = time.time()
        write_json_atomic(out / f'progress_{host}.json', prog)

    save()
    for j in pending:
        lo, hi = plan['intervals'][j]
        stem = out / f'shard_{j:04d}'
        cmd = ['nice', '-n', str(a.nice), str(binary), '--n', str(a.n), '--m', str(a.m),
               '--c-lo', str(lo), '--c-hi', str(hi), '--threads', str(a.threads)]
        prog['running'] = {str(j): {'c_lo': lo, 'c_hi': hi, 'started': jst()}}
        save()
        t0 = time.time()
        with open(str(stem) + '.part.json', 'w') as fo, open(str(stem) + '.err', 'w') as fe:
            p = subprocess.Popen(cmd, stdout=fo, stderr=fe, stdin=subprocess.DEVNULL)
            _, status, usage = os.wait4(p.pid, 0)
        try:
            if os.waitstatus_to_exitcode(status) != 0:
                raise RuntimeError('exit code %d' % os.waitstatus_to_exitcode(status))
            lines = pathlib.Path(str(stem) + '.part.json').read_text().splitlines()
            res = json.loads(lines[-1])
            sols = [json.loads(x) for x in lines[:-1] if x.startswith('{"solution"')]
            if (res['engine'], res['N'], res['M'], res['c_lo'], res['c_hi'], res['filtered']) != \
                    ('rpart_v2', a.n, a.m, lo, hi, True):
                raise RuntimeError('engine output does not match the shard')
            if res['solutions'] != len(sols):
                raise RuntimeError('solution lines do not match the counter')
            os.replace(str(stem) + '.part.json', str(stem) + '.json')
            write_json_atomic(str(stem) + '.done.json', {
                'shard': j, 'c_lo': lo, 'c_hi': hi, 'result': res, 'solutions': sols, 'host': host,
                'cpu_seconds': usage.ru_utime + usage.ru_stime, 'wall_seconds': time.time() - t0, 'finished': jst()})
            prog['done'] += 1
            if sols:
                prog.setdefault('solutions', []).extend(sols)
        except Exception as e:
            prog['failed'].append({'shard': j, 'error': str(e)})
            save()
            break
        save()
    prog['running'] = {}
    prog['state'] = 'error' if prog['failed'] else 'done'
    prog['finished'] = jst()
    save()
    print(json.dumps({k: prog[k] for k in ('host', 'state', 'done', 'selected', 'failed')}))
    return 1 if prog['failed'] else 0


def cmd_merge(a):
    plans, done = [], {}
    for d in a.inputs:
        d = pathlib.Path(d)
        plans.append(json.loads((d / 'plan.json').read_text()))
        for f in sorted(d.glob('shard_*.done.json')):
            r = json.loads(f.read_text())
            if r['shard'] in done:
                raise SystemExit(f'shard {r["shard"]} appears twice')
            done[r['shard']] = r
    plan = plans[0]
    if any(p != plan for p in plans):
        raise SystemExit('plans differ between inputs')
    S = plan['shards']
    if sorted(done) != list(range(S)):
        raise SystemExit('missing shards: %s' % sorted(set(range(S)) - set(done))[:20])
    prev = 0
    for j in range(S):
        lo, hi = done[j]['c_lo'], done[j]['c_hi']
        if [lo, hi] != plan['intervals'][j] or lo != prev + 1:
            raise SystemExit(f'shard {j} has wrong range')
        prev = hi
    if prev != plan['M']:
        raise SystemExit('c range does not end at M')
    tot = {k: sum(done[j]['result'][k] for j in range(S)) for k in COUNTERS}
    sols = [s for j in range(S) for s in done[j]['solutions']]
    if tot['solutions'] != len(sols):
        raise SystemExit('solution count mismatch')
    for s in sols:  # exact check of every reported solution
        x = s['solution']; q, a_, b, c, d = x['Q'], x['a'], x['b'], x['c'], x['d']
        import math
        ok = all(math.isqrt(u * u + v * v) ** 2 == u * u + v * v for u, v in ((a_, c), (b, c), (a_, d), (b, d)))
        x['verified_exact'] = ok
    outd = pathlib.Path(a.into); outd.mkdir(parents=True, exist_ok=False)
    summary = {'engine': plan['engine'], 'N': plan['N'], 'M': plan['M'],
               'box': 'all points (A/Q, C/Q) with Q <= N and max(|A|,|Q-A|,|C|,|Q-C|) <= M',
               'shards': S, 'all_shards_exactly_once': True, 'c_range_contiguous': True,
               'source_sha256': plan['source_sha256'], 'totals': tot, 'solutions': sols,
               'cpu_seconds': sum(done[j]['cpu_seconds'] for j in range(S)),
               'hosts': sorted({done[j]['host'] for j in range(S)}), 'merged': jst()}
    write_json_atomic(outd / 'summary.json', summary)
    print(json.dumps({k: summary[k] for k in ('N', 'M', 'shards', 'totals', 'cpu_seconds')}))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--n', type=int, required=True); r.add_argument('--m', type=int, required=True)
    r.add_argument('--shards', type=int, required=True); r.add_argument('--select', default='all')
    r.add_argument('--threads', type=int, default=1); r.add_argument('--output', required=True)
    r.add_argument('--binary', default='./rpart2'); r.add_argument('--nice', type=int, default=10)
    m = sp.add_parser('merge'); m.add_argument('inputs', nargs='+'); m.add_argument('--into', required=True)
    a = ap.parse_args()
    return cmd_run(a) if a.cmd == 'run' else cmd_merge(a)


if __name__ == '__main__':
    sys.exit(main())
