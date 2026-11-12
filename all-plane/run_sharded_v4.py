#!/usr/bin/env python3
"""Run disjoint v4 denominator intervals and verify exact coverage and witnesses.

Each worker is single-threaded. The driver requires the complete default v4
configuration, preserves every direct/transitive source hash, and independently
checks emitted distances with Python integers before writing a success summary.
"""
import argparse
import csv
import hashlib
import json
import math
import pathlib
import resource
import subprocess
import sys
import time


SOURCE_NAMES = (
    'all_plane_gap_v4.cpp', 'all_plane_gap_v3.cpp', 'all_plane_sieved.cpp',
    'research_v4_local/lifted_projection.hpp',
    'research_v4_math/critical_extension.hpp', 'run_sharded_v4.py',
)
REQUIRED_CONFIGURATION = {
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


def main():
    if sys.flags.optimize:
        raise RuntimeError('Verification requires assertions; run Python without -O or PYTHONOPTIMIZE.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound', type=int, default=100000)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--binary', default='./all_plane_gap_v4')
    parser.add_argument('--output', default='sweep_v4_100000')
    args = parser.parse_args()
    if not (24 <= args.bound <= 2000000 and 1 <= args.workers <= 32):
        raise ValueError('bound must be 24..2000000, workers 1..32')
    binary = pathlib.Path(args.binary).resolve()
    out = pathlib.Path(args.output).resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError('output directory is not empty; use a new output directory')
    out.mkdir(parents=True, exist_ok=True)
    kmax = args.bound // 24
    workers = min(args.workers, kmax)
    ends = [24*math.floor(kmax*math.sqrt((i+1)/workers)) for i in range(workers)]
    ends[-1] = args.bound
    intervals = []
    lower = 1
    for upper in ends:
        if upper < lower:
            raise ValueError('empty partition interval')
        intervals.append((lower, upper))
        lower = upper+1
    assert lower == args.bound+1
    source_paths = {name: binary.parent/name for name in SOURCE_NAMES}
    # Hash this running driver, even if it was invoked outside binary.parent.
    source_paths['run_sharded_v4.py'] = pathlib.Path(__file__).resolve()
    source_hashes_before = {name: hashlib.sha256(path.read_bytes()).hexdigest()
                            for name, path in source_paths.items()}
    binary_hash_before = hashlib.sha256(binary.read_bytes()).hexdigest()
    started = time.perf_counter()
    children = []
    for i, (lo, hi) in enumerate(intervals):
        stem = out/f'shard_{i:02d}'
        stdout = stem.with_suffix('.jsonl').open('w')
        stderr = stem.with_suffix('.stderr.txt').open('w')
        command = [str(binary), str(hi), '--min-denominator', str(lo),
                   '--audit-three', str(stem.with_suffix('.csv'))]
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        children.append((process, stdout, stderr, stem, lo, hi))
    remaining = set(range(len(children)))
    while remaining:
        for i in list(remaining):
            process, stdout, stderr, stem, lo, hi = children[i]
            code = process.poll()
            if code is None:
                continue
            stdout.close()
            stderr.close()
            remaining.remove(i)
            print(json.dumps({'shard': i, 'min': lo, 'max': hi, 'returncode': code,
                              'wall_elapsed': time.perf_counter()-started}), flush=True)
        if remaining:
            time.sleep(0.2)
    wall = time.perf_counter()-started
    results, points, solutions = [], [], []
    for process, stdout, stderr, stem, lo, hi in children:
        if process.returncode != 0:
            raise RuntimeError(f'worker failed; see {stem.with_suffix(".stderr.txt")}')
        messages = [json.loads(line) for line in stem.with_suffix('.jsonl').read_text().splitlines()
                    if line.strip()]
        stats = messages[-1]
        if stats.get('generator') != 'gap_CRT_v4' or stats.get('scope') != 'whole_plane':
            raise RuntimeError('Worker output is not from the expected whole-plane v4 generator.')
        for flag, expected in REQUIRED_CONFIGURATION.items():
            actual = stats.get(flag)
            if type(actual) is not type(expected) or actual != expected:
                raise RuntimeError('Worker did not use the required v4 configuration: ' + flag)
        assert stats['min_denominator'] == lo
        assert stats['max_minimal_coordinate_denominator'] == hi
        solutions.extend(m['solution'] for m in messages if 'solution' in m)
        results.append(stats)
        with stem.with_suffix('.csv').open(newline='') as source:
            for row in csv.DictReader(source):
                point = tuple(int(row[k]) for k in ('N', 'p', 'q', 'X', 'Y', 'A'))
                n, p, q, x, y, a = point
                assert lo <= n <= hi
                assert math.gcd(math.gcd(x, y), n) == 1
                assert x*x+y*y == a*a
                assert (n-x)**2+y*y == (a+p)**2
                assert x*x+(n-y)**2 == (a+q)**2
                points.append(point)
    assert len(points) == len(set(points)), 'duplicate three-distance point'
    points.sort()
    with (out/'three_candidates.csv').open('w', newline='') as destination:
        writer = csv.writer(destination)
        writer.writerow(('N', 'p', 'q', 'X', 'Y', 'A'))
        writer.writerows(points)
    counters = {}
    excluded = {'min_denominator', 'max_minimal_coordinate_denominator', 'wheel_modulus'}
    for key, value in results[0].items():
        if type(value) is int and key not in excluded:
            counters[key] = sum(r[key] for r in results)
    expected_dense = sum((12*k+1)*(12*k+2)//2 for k in range(1, kmax+1))
    assert counters['dense_even_pairs'] == expected_dense, 'denominator coverage mismatch'
    assert counters['primitive_three_candidates'] == len(points)
    actual_four = []
    for n, p, q, x, y, a in points:
        squared = (n-x)**2+(n-y)**2
        d = math.isqrt(squared)
        if d*d == squared:
            actual_four.append((n, x, y, a, d))
    reported_four = [tuple(s[k] for k in ('N', 'X', 'Y', 'A', 'D')) for s in solutions]
    assert counters['solutions'] == len(actual_four) == len(reported_four)
    assert sorted(actual_four) == sorted(reported_four), 'reported four-distance witness mismatch'
    sources = {name: hashlib.sha256(path.read_bytes()).hexdigest()
               for name, path in source_paths.items()}
    assert sources == source_hashes_before, 'source changed during computation'
    assert hashlib.sha256(binary.read_bytes()).hexdigest() == binary_hash_before, 'binary changed during computation'
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    summary = {
        'generator_version': 'v4', 'scope': 'whole_plane_primitive_minimal_coordinate_denominator',
        'min_denominator': 1, 'max_minimal_coordinate_denominator': args.bound,
        'workers': workers, 'threads_per_worker': 1, 'intervals': intervals,
        'all_intervals_completed': True, 'exact_interval_coverage_checked': True,
        'generator_configuration_checked': True, 'generator_configuration': REQUIRED_CONFIGURATION,
        'three_distance_witnesses_independently_verified': True,
        'four_distance_witnesses_independently_verified': True,
        'counts': counters, 'wall_seconds_including_startup': wall,
        'worker_cpu_seconds': usage.ru_utime+usage.ru_stime,
        'three_candidates_on_supporting_lines': sum(x in (0, n) or y in (0, n) for n, p, q, x, y, a in points),
        'three_candidates_strictly_inside': sum(0<x<n and 0<y<n for n, p, q, x, y, a in points),
        'solutions': solutions, 'source_sha256': sources, 'binary_sha256': binary_hash_before,
        'source_and_binary_unchanged_during_run': True,
    }
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
