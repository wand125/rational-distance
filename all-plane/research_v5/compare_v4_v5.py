#!/usr/bin/env python3
"""Run frozen v4 and v5 on the same denominator intervals and require identical
non-empty parameter sets, diagnostic three-distance points and counters.

Usage: compare_v4_v5.py --bound B --parts P --jobs J --out DIR (DIR must be new)
Binaries are taken from the parent directory: all_plane_gap_v4_local, all_plane_gap_v5_local.
"""
import argparse, hashlib, json, math, os, pathlib, subprocess, sys, time

TIMING = {'generator', 'reconstruction', 'seconds_setup', 'seconds_search', 'seconds_total_excluding_validation'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bound', type=int, required=True)
    ap.add_argument('--parts', type=int, default=4)
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.nice(10)
    root = pathlib.Path(__file__).resolve().parent.parent
    out = pathlib.Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    kmax = args.bound // 24
    ends = [24*math.floor(kmax*math.sqrt((i+1)/args.parts)) for i in range(args.parts)]
    ends[-1] = args.bound
    intervals, lo = [], 1
    for hi in ends:
        intervals.append((lo, hi)); lo = hi+1
    binaries = {'v4': root/'all_plane_gap_v4_local', 'v5': root/'all_plane_gap_v5_local'}
    jobs = [(i, name) for i in range(len(intervals)) for name in ('v4', 'v5')]
    running, results = {}, {}
    started = time.time()
    while jobs or running:
        while jobs and len(running) < args.jobs:
            i, name = jobs.pop(0)
            lo, hi = intervals[i]
            stem = out/f'{name}_{i:02d}'
            cmd = [str(binaries[name]), str(hi), '--min-denominator', str(lo),
                   '--audit-three', str(stem)+'_three.csv', '--audit-parameters', str(stem)+'_par.csv']
            p = subprocess.Popen(cmd, stdout=open(str(stem)+'.json', 'w'), stderr=open(str(stem)+'.err', 'w'))
            running[p.pid] = (i, name, p)
        pid, status, usage = os.wait4(-1, 0)
        i, name, p = running.pop(pid)
        if os.waitstatus_to_exitcode(status) != 0:
            raise RuntimeError(f'{name} interval {i} failed')
        stem = out/f'{name}_{i:02d}'
        stats = json.loads(pathlib.Path(str(stem)+'.json').read_text().splitlines()[-1])
        results[(i, name)] = {'stats': stats, 'cpu_seconds': usage.ru_utime+usage.ru_stime}
        print(json.dumps({'interval': i, 'binary': name, 'cpu': results[(i, name)]['cpu_seconds'],
                          'elapsed': time.time()-started}), flush=True)
    report = {'bound': args.bound, 'intervals': intervals, 'per_interval': [], 'identical': True}
    digests = {}
    for name in ('v4', 'v5'):
        for kind in ('par', 'three'):
            h = hashlib.sha256(); rows = 0
            for i in range(len(intervals)):
                data = pathlib.Path(f'{out}/{name}_{i:02d}_{kind}.csv').read_bytes()
                lines = data.splitlines(keepends=True)
                assert lines[0] in (b'N,p,q\n', b'N,p,q,X,Y,A\n')
                body = b''.join(lines[1:]); rows += len(lines)-1
                h.update(body)
            digests[(name, kind)] = (h.hexdigest(), rows)
    for kind in ('par', 'three'):
        same = digests[('v4', kind)] == digests[('v5', kind)]
        report[f'{kind}_rows'] = digests[('v5', kind)][1]
        report[f'{kind}_sha256'] = digests[('v5', kind)][0]
        report[f'{kind}_identical'] = same
        report['identical'] &= same
    totals = {'v4': {}, 'v5': {}}
    for i in range(len(intervals)):
        a, b = results[(i, 'v4')], results[(i, 'v5')]
        diff = {k: (a['stats'].get(k), b['stats'].get(k)) for k in set(a['stats']) | set(b['stats'])
                if k not in TIMING and a['stats'].get(k) != b['stats'].get(k)}
        report['identical'] &= not diff
        report['per_interval'].append({'interval': intervals[i], 'counter_differences': diff,
                                       'v4_cpu': a['cpu_seconds'], 'v5_cpu': b['cpu_seconds'],
                                       'v4_search': a['stats']['seconds_search'], 'v5_search': b['stats']['seconds_search']})
        for name, r in (('v4', a), ('v5', b)):
            for k, v in r['stats'].items():
                if type(v) is int and k not in ('min_denominator', 'max_minimal_coordinate_denominator'):
                    totals[name][k] = totals[name].get(k, 0)+v
            totals[name]['cpu_seconds'] = totals[name].get('cpu_seconds', 0)+r['cpu_seconds']
    report['totals'] = totals
    report['wall_seconds'] = time.time()-started
    (out/'report.json').write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps({k: report[k] for k in ('bound', 'identical', 'par_rows', 'three_rows', 'wall_seconds')}
                     | {'v4_cpu': totals['v4']['cpu_seconds'], 'v5_cpu': totals['v5']['cpu_seconds'],
                        'solutions': totals['v5']['solutions']}))
    return 0 if report['identical'] else 1


if __name__ == '__main__':
    sys.exit(main())
