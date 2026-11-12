#!/usr/bin/env python3
"""Run disjoint denominator intervals and audit complete range coverage.

Each worker runs the exact C++ gap-first search, using one thread. The driver
does not weaken any mathematical condition. All emitted three-distance points
are checked independently with Python integer arithmetic.
"""
import argparse
import csv
import hashlib
import json
import math
import pathlib
import resource
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bound', type=int, default=1000000)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--binary', default='./all_plane_gap_v2')
    parser.add_argument('--output', default='sweep_v2_2000000')
    args = parser.parse_args()
    if not (24 <= args.bound <= 2000000 and 1 <= args.workers <= 32):
        raise ValueError('bound must be 24..2000000, workers 1..32')
    binary = pathlib.Path(args.binary).resolve()
    out = pathlib.Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    # Approximately equal work for a quadratic cumulative candidate count.
    kmax = args.bound // 24
    workers = min(args.workers, kmax)
    ends = [24*math.floor(kmax*math.sqrt((i+1)/workers))
            for i in range(workers)]
    ends[-1] = args.bound
    intervals = []
    lower = 1
    for upper in ends:
        if upper < lower:
            raise ValueError('empty partition interval')
        intervals.append((lower,upper))
        lower = upper+1
    assert lower == args.bound+1
    source_names = ['all_plane_gap_v2.cpp', 'all_plane_sieved.cpp', 'run_sharded_v2.py']
    source_hashes_before = {name: hashlib.sha256((binary.parent/name).read_bytes()).hexdigest()
                            for name in source_names}
    binary_hash_before = hashlib.sha256(binary.read_bytes()).hexdigest()
    started = time.perf_counter()
    children = []
    for i,(lo,hi) in enumerate(intervals):
        stem = out/f'shard_{i:02d}'
        stdout = stem.with_suffix('.jsonl').open('w')
        stderr = stem.with_suffix('.stderr.txt').open('w')
        command = [str(binary), str(hi), '--min-denominator', str(lo),
                   '--audit-three', str(stem.with_suffix('.csv'))]
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        children.append((process,stdout,stderr,stem,lo,hi,command))
    remaining = set(range(len(children)))
    while remaining:
        for i in list(remaining):
            process,stdout,stderr,stem,lo,hi,command = children[i]
            code = process.poll()
            if code is None:
                continue
            stdout.close()
            stderr.close()
            remaining.remove(i)
            print(json.dumps({'shard':i,'min':lo,'max':hi,'returncode':code,
                              'wall_elapsed':time.perf_counter()-started}),flush=True)
        if remaining:
            time.sleep(0.2)
    wall = time.perf_counter()-started
    results, points, solutions = [], [], []
    for process,stdout,stderr,stem,lo,hi,command in children:
        if process.returncode != 0:
            raise RuntimeError(f'worker failed; see {stem.with_suffix(".stderr.txt")}')
        messages = [json.loads(line) for line in stem.with_suffix('.jsonl').read_text().splitlines() if line.strip()]
        stats = messages[-1]
        assert stats['min_denominator'] == lo
        assert stats['max_minimal_coordinate_denominator'] == hi
        solutions.extend(m['solution'] for m in messages if 'solution' in m)
        results.append(stats)
        with stem.with_suffix('.csv').open() as f:
            for row in csv.DictReader(f):
                point = tuple(int(row[k]) for k in ['N','p','q','X','Y','A'])
                n,p,q,x,y,a = point
                assert lo <= n <= hi
                assert math.gcd(math.gcd(x,y),n) == 1
                assert x*x+y*y == a*a
                assert (n-x)**2+y*y == (a+p)**2
                assert x*x+(n-y)**2 == (a+q)**2
                points.append(point)
    assert len(points) == len(set(points)), 'duplicate three-distance point'
    points.sort()
    with (out/'three_candidates.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['N','p','q','X','Y','A'])
        writer.writerows(points)
    counters = {}
    # The programs report only scalar metadata and integer additive counters.
    excluded = {'min_denominator','max_minimal_coordinate_denominator','wheel_modulus'}
    for key,value in results[0].items():
        if type(value) is int and key not in excluded:
            counters[key] = sum(r[key] for r in results)
    expected_dense = sum((12*k+1)*(12*k+2)//2 for k in range(1,kmax+1))
    assert counters['dense_even_pairs'] == expected_dense, 'denominator coverage mismatch'
    assert counters['primitive_three_candidates'] == len(points)
    actual_four = [p for p in points if math.isqrt((p[0]-p[3])**2+(p[0]-p[4])**2)**2
                   == (p[0]-p[3])**2+(p[0]-p[4])**2]
    assert counters['solutions'] == len(actual_four) == len(solutions)
    resource_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    summary = {'scope':'whole_plane_primitive_minimal_coordinate_denominator',
               'min_denominator':1,'max_minimal_coordinate_denominator':args.bound,
               'workers':workers,'threads_per_worker':1,'intervals':intervals,
               'all_intervals_completed':True,'exact_interval_coverage_checked':True,
               'three_distance_witnesses_independently_verified':True,
               'counts':counters,'wall_seconds_including_startup':wall,
               'worker_cpu_seconds':resource_usage.ru_utime+resource_usage.ru_stime,
               'three_candidates_on_supporting_lines':sum(x in (0,n) or y in (0,n) for n,p,q,x,y,a in points),
               'three_candidates_strictly_inside':sum(0<x<n and 0<y<n for n,p,q,x,y,a in points),
               'solutions':solutions}
    sources = {}
    for filename in source_names:
        path = binary.parent/filename
        if path.exists():
            sources[filename] = hashlib.sha256(path.read_bytes()).hexdigest()
    assert sources == source_hashes_before, 'source changed during computation'
    assert hashlib.sha256(binary.read_bytes()).hexdigest() == binary_hash_before, 'binary changed during computation'
    summary['source_sha256'] = sources
    summary['binary_sha256'] = binary_hash_before
    summary['source_and_binary_unchanged_during_run'] = True
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
