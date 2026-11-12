#!/usr/bin/env python3
"""Bounded certificate for the next digit of critical-content projections.

Everything here uses exact modular integer arithmetic.  The direct checker
enumerates (a,r) modulo p^2 and tests the uneliminated three-radius equations;
it does not use the quartic, a modular square root of -1, or its denominator.
"""
import argparse
import json
import math
import time
from pathlib import Path


def quartic(modulus, iota, r):
    t = r*(r-1) % modulus
    return ((2*r-1)*(2*t-1)-2*iota*t*(t+2)) % modulus


def derivative(p, iota, r):
    t = r*(r-1) % p
    return (12*t-4*iota*(t+1)*(2*r-1)) % p


def projection(p, exponent):
    modulus = p**exponent
    squares = {a*a % modulus for a in range(modulus)}
    iotas = [i for i in range(modulus) if i*i % modulus == modulus-1]
    assert len(iotas) == 2
    allowed = set()
    for r in range(modulus):
        if r % p in (0, 1):
            continue
        if any((r-1-iota*r) % p and quartic(modulus, iota, r) in squares
               for iota in iotas):
            allowed.add(r)
    return allowed, iotas, squares


def direct_normalized_projection(p):
    """O(p^4) independent reference, deliberately limited to p<=41."""
    modulus = p*p
    squares = {a*a % modulus for a in range(modulus)}
    allowed = set()
    for r in range(modulus):
        if r % p in (0, 1):
            continue
        for a in range(modulus):
            # b=a+2r, c=a+2(r-1).  Use their square differences directly.
            first = 4*r*(a+r) % modulus
            second = 4*(r-1)*(a+r-1) % modulus
            if first % p == 0 or second % p == 0:
                continue
            if (first*first+second*second) % modulus:
                continue
            d_square = (a*a+first+second) % modulus
            if d_square in squares:
                allowed.add(r)
                break
    return allowed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(Path(__file__).with_name('critical_lift_certificate.json')))
    args = parser.parse_args()
    start = time.process_time()
    rows = []
    for p in (5, 13, 17, 29, 37, 41, 53, 61):
        old, small_iotas, small_squares = projection(p, 1)
        new, iotas, _ = projection(p, 2)
        assert {r % p for r in new} <= old
        old_lifts = {r for r in range(p*p) if r % p in old}
        reference_checked = p <= 41
        if reference_checked:
            reference = direct_normalized_projection(p)
            assert reference == new, (p, sorted(reference ^ new))
        branches = []
        for r0 in sorted(old):
            possible = [i for i in small_iotas
                        if (r0-1-i*r0) % p and quartic(p, i, r0) in small_squares]
            unit_branch = any(quartic(p, i, r0) != 0 for i in possible)
            lifts = sorted(r for r in new if r % p == r0)
            if unit_branch:
                assert len(lifts) == p
            else:
                zero_branches = []
                for i in iotas:
                    if (r0-1-i*r0) % p and quartic(p, i, r0) == 0:
                        fp = derivative(p, i, r0)
                        assert fp != 0
                        k = (-quartic(p*p, i, r0)//p*pow(fp, -1, p)) % p
                        root = r0+p*k
                        assert quartic(p*p, i, root) == 0
                        zero_branches.append({'iota_mod_p2': i, 'derivative_mod_p': fp,
                                              'lift_digit': k, 'root_mod_p2': root})
                assert set(lifts) == {x['root_mod_p2'] for x in zero_branches}
                branches.append({'r_mod_p': r0, 'allowed_lifts': lifts,
                                 'zero_branches': zero_branches})
        rows.append({'prime': p, 'old_residue_count': len(old),
                     'old_lift_count': len(old_lifts), 'new_lift_count': len(new),
                     'additional_rejected': len(old_lifts-new),
                     'uneliminated_enumeration_matches': reference_checked,
                     'branch_only_residues': branches})
    row41 = next(row for row in rows if row['prime'] == 41)
    assert row41['additional_rejected'] == 240
    assert {x['r_mod_p']: x['allowed_lifts'] for x in row41['branch_only_residues']} == {
        7: [48], 12: [750], 15: [1327], 27: [355], 30: [932], 35: [1634]}
    old_local_primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 43)
    orders = {p: next(e for e in range(1, p) if pow(41, e, p) == 1)
              for p in old_local_primes}
    period = math.lcm(*orders.values())
    assert period == 55440
    r_example = 384*pow(119, -1, 1681) % 1681
    assert r_example == 314 and r_example % 41 == 27
    data = {'scope': 'whole_plane_critical_content_next_digit',
            'necessary_modulus': 'p^(2s), where v_p(n)=2s and v_p(G)=s>0',
            'zero_square_roots_included': True, 'certificates': rows,
            'v3_strict_example': {'n': 20172, 'u': 15744, 'v': 10865,
                                  'r_mod_1681': r_example, 'required_r_mod_1681': 355},
            'infinite_old_filter_family': {
                's': '1+55440k, k>=0', 'n': '12*41^(2s)',
                'u': '384*41^s', 'v': '265*41^s',
                'orders_of_41_at_other_old_local_primes': orders,
                'period': period,
                'scope': 'passes v3 arithmetic parameter predicates over integers; not a reconstructed point; not a machine-range claim'},
            'cpu_seconds': time.process_time()-start,
            'mathematical_novelty_claimed': False}
    with open(args.output, 'w', encoding='utf-8') as out:
        json.dump(data, out, indent=2)
        out.write('\n')
    print(json.dumps({'output': args.output, 'cpu_seconds': data['cpu_seconds'],
                      'direct_reference_primes': [5, 13, 17, 29, 37, 41],
                      'critical_41_rejected_lifts': 240,
                      'all_assertions_passed': True}, indent=2))


if __name__ == '__main__':
    main()
