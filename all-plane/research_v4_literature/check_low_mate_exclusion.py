#!/usr/bin/env python3
"""Independent finite factor-pair checks for the low-mate family theorem."""
import argparse
import itertools
import json
import math
from pathlib import Path


def factors(value):
    out = []
    prime = 2
    while prime * prime <= value:
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        if exponent:
            out.append((prime, exponent))
        prime = 3 if prime == 2 else prime + 2
    if value > 1:
        out.append((value, 1))
    return out


def predicted_count(k):
    bad = 0
    good_choices = 1
    for prime, exponent in factors(k):
        if prime % 8 in (3, 5):
            bad += 1
        else:
            good_choices *= 2 * exponent + 1
    if bad:
        return 2 ** (bad - 1) * good_choices
    return (good_choices - 1) // 2


def direct_mates(k):
    """Enumerate every positive factor a < k of k^2, not the count formula."""
    bad_primes = [p for p, _ in factors(k) if p % 8 in (3, 5)]
    mates = set()
    for a in range(1, k, 2):
        if k * k % a:
            continue
        c = k * k // a
        z, w = (c - a) // 2, (c + a) // 2
        if any(z % p == 0 for p in bad_primes):
            continue
        if not z > 0 or k*k + z*z != w*w:
            raise RuntimeError('Factor-pair identity failed')
        mates.add(z)
    return mates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-leg', type=int, default=5001)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('low_mate_certificate.json'))
    args = parser.parse_args()
    all_sets = {}
    tested_mates = 0
    small_sets = []
    for k in range(1, args.max_leg + 1, 2):
        actual = direct_mates(k)
        predicted = predicted_count(k)
        if len(actual) != predicted:
            raise RuntimeError('Count mismatch at ' + str(k))
        if k > 1 and max(actual) != (k*k - 1)//2:
            raise RuntimeError('Maximum mismatch at ' + str(k))
        all_sets[k] = actual
        tested_mates += len(actual)
        if len(actual) <= 2 and k <= 501:
            small_sets.append(k)
    pairs = 0
    nonempty = 0
    for k, ell in itertools.combinations(small_sets, 2):
        intersection = all_sets[k] & all_sets[ell]
        if len(intersection) >= 2:
            raise RuntimeError('Two common mates for distinct low-mate legs')
        pairs += 1
        nonempty += bool(intersection)
    # An illustrative primitive-denominator row, not a complete sweep.
    side = 99984
    odd_columns = list(range(1, side//2, 2))
    counts = {
        'odd_columns_below_midline': len(odd_columns),
        'single_leg_count_below_two': 0,
        'additional_both_counts_exactly_two': 0,
    }
    for x in odd_columns:
        left, right = predicted_count(x), predicted_count(side-x)
        if min(left, right) < 2:
            counts['single_leg_count_below_two'] += 1
        elif left == right == 2:
            counts['additional_both_counts_exactly_two'] += 1
    record = {
        'status': 'passed',
        'max_odd_leg': args.max_leg,
        'odd_legs_checked': len(all_sets),
        'positive_admissible_mates_checked': tested_mates,
        'count_formula_and_maximum_match': True,
        'small_low_mate_distinct_leg_pairs_checked': pairs,
        'pairs_with_one_common_mate': nonempty,
        'pairs_with_two_common_mates': 0,
        'examples': {str(k): sorted(all_sets[k]) for k in [7, 15, 45, 49, 75, 143, 225]},
        'illustrative_fixed_side': side,
        'illustrative_column_counts': counts,
        'caveat': 'Finite checks support the written general proof. Column counts precede other sieves and are not a measured search speedup.',
    }
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
