#!/usr/bin/env python3
"""Resumable, multi-host sharded driver for the whole-plane v5 search.

run   : execute the selected shards of a fixed partition of [1, bound] into
        many small denominator intervals. Finished shards are recorded with a
        *.done.json marker after their output was checked; a restart skips them,
        so an interruption loses at most the shards that were running.
merge : combine one or more run directories (e.g. copied from several hosts),
        require every shard exactly once, check configuration, coverage and
        counters, and verify every emitted witness with Python integers.

Workers are single-threaded processes. Python must run without -O.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import pathlib
import socket
import subprocess
import sys
import time

SOURCE_NAMES = (
    'all_plane_gap_v5.cpp', 'wide_arith_v5.hpp', 'all_plane_gap_v4.cpp', 'all_plane_gap_v3.cpp',
    'all_plane_sieved.cpp', 'research_v4_local/lifted_projection.hpp',
    'research_v4_math/critical_extension.hpp', 'run_sharded_v5.py',
)
REQUIRED_CONFIGURATION = {
    'generator': 'gap_CRT_v5',
    'reconstruction': 'exact_u256_v5',
    'scope': 'whole_plane',
    'reference_v3': False,
    'gamma_stage': 'late',
    'skip_redundant19': True,
    'additional_filters': True,
    'endpoint_pruning': True,
    'fifth_power_condition': 'early',
    'prime_power_projections': True,
    'critical_extension': True,
    'fast_fifth_mark': True,
}
MAX_BOUND = 10_000_000


def jst(t=None):
    return time.strftime('%Y-%m-%d %H:%M:%S JST', time.gmtime((time.time() if t is None else t) + 9*3600))


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def partition(bound, shards, lower=1):
    """Disjoint intervals covering [lower,bound]; work per n grows ~linearly, so
    cut points follow sqrt of the cumulative work. Every cut is a multiple of 24
    (only 24|N is searched). The first searched side is N=24(k0+1) >= lower."""
    if not 1 <= lower <= bound:
        raise ValueError('need 1 <= lower <= bound')
    k0, kmax = (lower+23) // 24 - 1, bound // 24
    ends = []
    for i in range(shards):
        e = bound if i == shards-1 else 24*math.floor(math.sqrt(k0*k0+(kmax*kmax-k0*k0)*(i+1)/shards))
        ends.append(e)
    intervals, lower0 = [], lower
    for upper in ends:
        if upper < lower or (upper < bound and upper % 24):
            raise ValueError('empty or misaligned partition interval; use fewer shards')
        intervals.append((lower, upper))
        lower = upper+1
    assert lower == bound+1 and intervals[0][0] == lower0
    return intervals


def parse_select(text, shards):
    # "all" | "r1,r2/m" (indices with j % m in {r..}) | "a-b" (inclusive range)
    if text == 'all':
        return list(range(shards))
    if '/' in text:
        residues, modulus = text.split('/')
        m = int(modulus)
        keep = {int(r) for r in residues.split(',')}
        assert all(0 <= r < m for r in keep)
        return [j for j in range(shards) if j % m in keep]
    a, b = text.split('-')
    return list(range(int(a), int(b)+1))


def source_hashes(root):
    hashes = {name: sha256(root/name) for name in SOURCE_NAMES if name != 'run_sharded_v5.py'}
    hashes['run_sharded_v5.py'] = sha256(pathlib.Path(__file__).resolve())
    return hashes


def check_stats(stats, lo, hi):
    for flag, expected in REQUIRED_CONFIGURATION.items():
        actual = stats.get(flag)
        if type(actual) is not type(expected) or actual != expected:
            raise RuntimeError('worker did not use the required v5 configuration: ' + flag)
    if stats['min_denominator'] != lo or stats['max_minimal_coordinate_denominator'] != hi:
        raise RuntimeError('worker interval mismatch')


def read_shard(stem, lo, hi):
    lines = [json.loads(x) for x in stem.with_suffix('.jsonl').read_text().splitlines() if x.strip()]
    stats = lines[-1]
    check_stats(stats, lo, hi)
    solutions = [m['solution'] for m in lines if 'solution' in m]
    return stats, solutions


def write_json_atomic(path, data):
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=1) + '\n')
    os.replace(tmp, path)


def cmd_run(args):
    if args.nice:
        os.nice(args.nice)
    binary = pathlib.Path(args.binary).resolve()
    root = binary.parent
    out = pathlib.Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if not (24 <= args.bound <= MAX_BOUND):
        raise ValueError('bound must be 24..10000000')
    intervals = partition(args.bound, args.shards, args.lower)
    hashes = source_hashes(root)
    plan = json.loads(json.dumps({'lower': args.lower, 'bound': args.bound, 'shards': args.shards,
                                  'intervals': intervals, 'source_sha256': hashes}))
    plan_path = out/'plan.json'
    if plan_path.exists():
        old = json.loads(plan_path.read_text())
        if old != plan:
            raise RuntimeError('existing plan.json differs (bound, shards or sources changed); use a new output directory')
    else:
        write_json_atomic(plan_path, plan)
    selected = parse_select(args.select, args.shards)
    binary_hash = sha256(binary)
    host = socket.gethostname()
    pending = [j for j in selected if not (out/f'shard_{j:04d}.done.json').exists()]
    started = time.time()
    progress = {'host': host, 'lower': args.lower, 'select': args.select, 'selected': len(selected), 'shards_total': args.shards,
                'bound': args.bound, 'workers': args.workers, 'started': jst(started), 'started_epoch': started,
                'skipped_already_done': len(selected)-len(pending), 'done': len(selected)-len(pending),
                'failed': [], 'running': {}, 'cpu_seconds_this_session': 0.0, 'state': 'running',
                'binary_sha256': binary_hash, 'pid': os.getpid()}

    def save_progress():
        progress['updated'] = jst()
        progress['updated_epoch'] = time.time()
        write_json_atomic(out/f'progress_{host}.json', progress)

    save_progress()
    running = {}
    queue = list(pending)
    while queue or running:
        while queue and len(running) < args.workers and not progress['failed']:
            j = queue.pop(0)
            lo, hi = intervals[j]
            stem = out/f'shard_{j:04d}'
            part = out/f'shard_{j:04d}.part'
            stdout = part.with_suffix('.part.jsonl').open('w')
            stderr = part.with_suffix('.part.stderr.txt').open('w')
            command = [str(binary), str(hi), '--min-denominator', str(lo),
                       '--audit-three', str(part.with_suffix('.part.csv'))]
            proc = subprocess.Popen(command, stdout=stdout, stderr=stderr)
            running[proc.pid] = (j, proc, stdout, stderr, time.time())
            progress['running'][str(j)] = {'min': lo, 'max': hi, 'pid': proc.pid, 'started': jst()}
            save_progress()
        if not running:
            break
        pid, status, usage = os.wait4(-1, 0)
        if pid not in running:
            continue
        j, proc, stdout, stderr, t0 = running.pop(pid)
        stdout.close(); stderr.close()
        code = os.waitstatus_to_exitcode(status)
        lo, hi = intervals[j]
        part = out/f'shard_{j:04d}.part'
        stem = out/f'shard_{j:04d}'
        cpu = usage.ru_utime + usage.ru_stime
        progress['cpu_seconds_this_session'] += cpu
        del progress['running'][str(j)]
        try:
            if code != 0:
                raise RuntimeError(f'exit code {code}')
            for suffix in ('.jsonl', '.csv', '.stderr.txt'):
                os.replace(part.with_suffix('.part' + suffix), stem.with_suffix(suffix))
            stats, solutions = read_shard(stem, lo, hi)
            done = {'shard': j, 'min': lo, 'max': hi, 'returncode': code, 'host': host,
                    'binary_sha256': binary_hash, 'source_sha256': hashes,
                    'cpu_seconds': cpu, 'wall_seconds': time.time()-t0, 'finished': jst(),
                    'max_rss_raw': usage.ru_maxrss, 'solutions': solutions,
                    'csv_sha256': sha256(stem.with_suffix('.csv')),
                    'jsonl_sha256': sha256(stem.with_suffix('.jsonl'))}
            write_json_atomic(out/f'shard_{j:04d}.done.json', done)
            progress['done'] += 1
            if solutions:
                progress.setdefault('solutions', []).extend(solutions)
        except Exception as error:  # keep the evidence, stop launching new shards
            progress['failed'].append({'shard': j, 'error': str(error)})
        save_progress()
    if source_hashes(root) != hashes or sha256(binary) != binary_hash:
        progress['failed'].append({'shard': None, 'error': 'source or binary changed during run'})
    progress['state'] = 'error' if progress['failed'] else 'done'
    progress['finished'] = jst()
    save_progress()
    print(json.dumps({k: progress[k] for k in ('host', 'state', 'done', 'selected', 'failed', 'cpu_seconds_this_session')}))
    return 1 if progress['failed'] else 0


def cmd_merge(args):
    dirs = [pathlib.Path(d).resolve() for d in args.inputs]
    plans = [json.loads((d/'plan.json').read_text()) for d in dirs]
    if any(p != plans[0] for p in plans):
        raise RuntimeError('input directories use different plans')
    plan = plans[0]
    lower, bound, intervals = plan['lower'], plan['bound'], [tuple(x) for x in plan['intervals']]
    found = {}
    for d in dirs:
        for marker in sorted(d.glob('shard_*.done.json')):
            done = json.loads(marker.read_text())
            j = done['shard']
            if done['source_sha256'] != plan['source_sha256']:
                raise RuntimeError(f'shard {j} was produced from different sources')
            if j in found:
                raise RuntimeError(f'shard {j} present twice ({found[j][0]} and {d}); keep one copy')
            found[j] = (d, done)
    missing = [j for j in range(len(intervals)) if j not in found]
    if missing:
        raise RuntimeError(f'{len(missing)} shards missing, first {missing[:10]}')
    results, points, solutions, hosts, cpu = [], [], [], {}, 0.0
    for j, (lo, hi) in enumerate(intervals):
        d, done = found[j]
        assert (done['min'], done['max']) == (lo, hi)
        stem = d/f'shard_{j:04d}'
        assert sha256(stem.with_suffix('.csv')) == done['csv_sha256']
        assert sha256(stem.with_suffix('.jsonl')) == done['jsonl_sha256']
        stats, shard_solutions = read_shard(stem, lo, hi)
        results.append(stats)
        solutions.extend(shard_solutions)
        hosts.setdefault(done['host'], {'shards': 0, 'cpu_seconds': 0.0, 'binary_sha256': set()})
        hosts[done['host']]['shards'] += 1
        hosts[done['host']]['cpu_seconds'] += done['cpu_seconds']
        hosts[done['host']]['binary_sha256'].add(done['binary_sha256'])
        cpu += done['cpu_seconds']
        with stem.with_suffix('.csv').open(newline='') as source:
            for row in csv.DictReader(source):
                n, p, q, x, y, a = (int(row[k]) for k in ('N', 'p', 'q', 'X', 'Y', 'A'))
                assert lo <= n <= hi and n % 24 == 0
                assert math.gcd(math.gcd(x, y), n) == 1
                assert x*x+y*y == a*a
                assert (n-x)**2+y*y == (a+p)**2
                assert x*x+(n-y)**2 == (a+q)**2
                points.append((n, p, q, x, y, a))
    assert len(points) == len(set(points)), 'duplicate three-distance point'
    points.sort()
    final = pathlib.Path(args.into).resolve()
    if final.exists() and any(final.iterdir()):
        raise RuntimeError('merge destination is not empty')
    final.mkdir(parents=True, exist_ok=True)
    with (final/'three_candidates.csv').open('w', newline='') as destination:
        writer = csv.writer(destination)
        writer.writerow(('N', 'p', 'q', 'X', 'Y', 'A'))
        writer.writerows(points)
    counters = {}
    excluded = {'min_denominator', 'max_minimal_coordinate_denominator'}
    for key, value in results[0].items():
        if type(value) is int and not isinstance(value, bool) and key not in excluded:
            counters[key] = sum(r[key] for r in results)
    k0, kmax = (lower+23) // 24 - 1, bound // 24
    expected_dense = sum((12*k+1)*(12*k+2)//2 for k in range(k0+1, kmax+1))
    assert counters['dense_even_pairs'] == expected_dense, 'denominator coverage mismatch'
    assert counters['sides'] == kmax-k0
    assert counters['primitive_three_candidates'] == len(points)
    actual_four = []
    for n, p, q, x, y, a in points:
        squared = (n-x)**2+(n-y)**2
        d = math.isqrt(squared)
        if d*d == squared:
            actual_four.append((n, x, y, a, d))
    reported_four = [tuple(int(s[k]) for k in ('N', 'X', 'Y', 'A', 'D')) for s in solutions]
    assert counters['solutions'] == len(actual_four) == len(reported_four)
    assert sorted(actual_four) == sorted(reported_four), 'reported four-distance witness mismatch'
    summary = {
        'generator_version': 'v5', 'scope': 'whole_plane_primitive_minimal_coordinate_denominator',
        'min_denominator': lower, 'max_minimal_coordinate_denominator': bound,
        'shards': len(intervals), 'all_shards_completed_exactly_once': True,
        'exact_interval_coverage_checked': True, 'generator_configuration': REQUIRED_CONFIGURATION,
        'three_distance_witnesses_independently_verified': True,
        'four_distance_witnesses_independently_verified': True,
        'counts': counters, 'worker_cpu_seconds': cpu,
        'hosts': {h: {**v, 'binary_sha256': sorted(v['binary_sha256'])} for h, v in hosts.items()},
        'three_candidates_on_supporting_lines': sum(x in (0, n) or y in (0, n) for n, p, q, x, y, a in points),
        'three_candidates_strictly_inside': sum(0 < x < n and 0 < y < n for n, p, q, x, y, a in points),
        'solutions': solutions, 'source_sha256': plan['source_sha256'], 'inputs': [str(d) for d in dirs],
        'merged': jst(),
    }
    (final/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: summary[k] for k in ('max_minimal_coordinate_denominator', 'shards', 'worker_cpu_seconds')}
                     | {'parameter_pairs': counters['exact_parameter_pairs'],
                        'three': counters['primitive_three_candidates'], 'solutions': counters['solutions']}))
    return 0


def main():
    if sys.flags.optimize:
        raise RuntimeError('Verification requires assertions; run Python without -O or PYTHONOPTIMIZE.')
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    run = sub.add_parser('run')
    run.add_argument('--bound', type=int, required=True)
    run.add_argument('--lower', type=int, default=1)
    run.add_argument('--shards', type=int, required=True)
    run.add_argument('--workers', type=int, default=2)
    run.add_argument('--select', default='all')
    run.add_argument('--binary', default='./all_plane_gap_v5_local')
    run.add_argument('--output', required=True)
    run.add_argument('--nice', type=int, default=10)
    merge = sub.add_parser('merge')
    merge.add_argument('inputs', nargs='+')
    merge.add_argument('--into', required=True)
    args = parser.parse_args()
    if args.command == 'run':
        if not (1 <= args.workers <= 32):
            raise ValueError('workers 1..32')
        return cmd_run(args)
    return cmd_merge(args)


if __name__ == '__main__':
    sys.exit(main())
