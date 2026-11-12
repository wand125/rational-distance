#!/usr/bin/env python3
"""Independent small-bound audit, using coordinate loops and math.isqrt.

No Euclidean triple generator, CSR graph, residue buckets, or floating point
are used to generate the reference set. The C++ audit set is compared exactly.
"""
import argparse
import csv
import json
import math
import pathlib
import subprocess
import tempfile
import time


def is_square(n):
    return math.isqrt(n) ** 2 == n


def good(n):
    """Trial division, independent of the production prime sieve."""
    if n <= 0 or n % 2 == 0:
        return False
    p = 3
    while p * p <= n:
        if n % p == 0:
            if p % 8 in (3, 5):
                return False
            while n % p == 0:
                n //= p
        p += 2
    return n % 8 not in (3, 5)


def reference_wedges(bound):
    result = set()
    for n in range(4, bound + 1, 4):
        for a in range(1, n // 2 + 1, 2):
            for c in range(4, n // 2 + 1, 4):
                if math.gcd(math.gcd(a, c), n) != 1:
                    continue
                if is_square(a*a + c*c) and is_square(a*a + (n-c)**2):
                    result.add((n, a, c))
    return result


def good_wedge(t):
    n, a, c = t
    b, d = n-a, n-c
    return all(good(x) for x in
               (math.gcd(a, c), math.gcd(a, d), b+c, b+d))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', default='./square_search')
    parser.add_argument('--limit', type=int, default=1200)
    args = parser.parse_args()
    binary = str(pathlib.Path(args.binary).resolve())
    started = time.perf_counter()
    reference = reference_wedges(args.limit)
    report = {'limit': args.limit, 'method': 'coordinate loops plus math.isqrt',
              'reference_primitive_two_distance_wedges': len(reference),
              'configurations': []}
    with tempfile.TemporaryDirectory() as directory:
        audit = pathlib.Path(directory) / 'wedges.csv'
        for step, use_good in [(4, False), (12, False), (24, False),
                               (4, True), (24, True)]:
            command = [binary, '--limit', str(args.limit), '--side-multiple',
                       str(step), '--no-sieve', '--audit-wedges', str(audit)]
            if not use_good:
                command.append('--no-good')
            actual_stats = json.loads(subprocess.check_output(command, text=True))
            with audit.open() as f:
                rows = [tuple(map(int, (r['N'], r['X'], r['Y'])))
                        for r in csv.DictReader(f)]
            actual = set(rows)
            expected = {t for t in reference if t[0] % step == 0
                        and (not use_good or good_wedge(t))}
            assert len(rows) == len(actual), 'duplicate audit rows'
            assert actual == expected, {
                'step': step, 'good': use_good,
                'missing': sorted(expected-actual)[:20],
                'unexpected': sorted(actual-expected)[:20]}
            exact_s = sum(is_square((n-a)**2+c*c) for n,a,c in expected)
            exact_t = sum(is_square((n-a)**2+(n-c)**2) for n,a,c in expected)
            three = sum(is_square((n-a)**2+c*c) or
                        is_square((n-a)**2+(n-c)**2) for n,a,c in expected)
            four = sum(is_square((n-a)**2+c*c) and
                       is_square((n-a)**2+(n-c)**2) for n,a,c in expected)
            assert actual_stats['exact_square_S'] == exact_s
            assert actual_stats['exact_square_T'] == exact_t
            assert actual_stats['three_distances_after_sieves'] == three
            assert actual_stats['solutions'] == four
            report['configurations'].append({
                'side_multiple': step, 'good_filter': use_good,
                'exact_candidate_set_match': True, 'wedges': len(expected),
                'three_distances': three, 'four_distances': four})
    report['seconds'] = time.perf_counter() - started
    report['passed'] = True
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
