#!/usr/bin/env python3
"""Construct critical masks with one quartic test per field residue.

No modular inversions or square-root sign enumeration are needed by the fast
method. The test driver compares it with an eliminated-radius reference and
the existing square-class pattern method, including zero square roots.
"""
import argparse
import json
import time

from critical_projection import critical_projection, primes


def quadratic_table(p):
    squares = [False] * p
    iota = None
    for a in range(p):
        square = a*a % p
        squares[square] = True
        if square == p-1 and iota is None:
            iota = a
    if p % 4 != 1 or iota is None:
        raise ValueError('an odd prime p == 1 mod4 is required')
    return squares, iota


def quartic(p, iota, r):
    t = r*(r-1) % p
    # Equivalent to 4r^3-6r^2+1 - 2*iota*r(r-1)*(r^2-r+2).
    return ((2*r-1)*(2*t-1) - 2*iota*t*(t+2)) % p


def fast_projection(p):
    squares, iota = quadratic_table(p)
    allowed = set()
    for r in range(2, p):
        denominator = (r-1-iota*r) % p
        if denominator == 0:
            continue
        if squares[quartic(p, iota, r)]:
            allowed.add(r)
            # F_i(1-r)=-F_{-i}(r), with -1 a square; this covers -i.
            allowed.add((1-r) % p)
    return allowed


def radius_reference(p):
    squares, chosen = quadratic_table(p)
    allowed = set()
    zero_radius_examples = []
    invalid_denominator_cases = []
    for iota in (chosen, (-chosen) % p):
        for r in range(2, p):
            u, v = r, r-1
            denominator = (v-iota*u) % p
            numerator = (iota*u*u-v*v) % p
            f = quartic(p, iota, r)
            if denominator == 0:
                # The affine equation has nonzero RHS and no solution here.
                assert numerator != 0
                assert f == numerator*numerator % p
                invalid_denominator_cases.append({'iota': iota, 'r': r})
                continue
            a = numerator*pow(denominator, -1, p) % p
            b, c = (a+2*u) % p, (a+2*v) % p
            beta = (b*b-a*a) % p
            assert beta == 4*u*v*pow(denominator, -1, p) % p
            assert beta != 0
            assert (c*c-a*a) % p == iota*beta % p
            d_square = (b*b+c*c-a*a) % p
            assert d_square == (a*a+4*a*(2*r-1)+4*(2*r*r-2*r+1)) % p
            assert f == d_square*denominator*denominator % p
            if squares[d_square]:
                allowed.add(r)
                if not a or not b or not c or not d_square:
                    zero_radius_examples.append({'iota': iota, 'r': r,
                                                 'first_three_roots': [a, b, c],
                                                 'fourth_square': d_square})
    return allowed, zero_radius_examples, invalid_denominator_cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=257)
    parser.add_argument('--output', default='fast_critical_projection_certificate.json')
    args = parser.parse_args()
    if not 5 <= args.limit <= 10000:
        parser.error('limit must be5..10000; default257 is the documented verification')
    started = time.process_time()
    results = []
    for p in primes(args.limit):
        if p % 4 != 1:
            continue
        fast = fast_projection(p)
        reference, zero, invalid = radius_reference(p)
        old = set(critical_projection(p)['allowed_u_over_gap'])
        assert fast == reference == old
        results.append({
            'prime': p, 'all_three_methods_match': True,
            'allowed_count': len(fast), 'allowed_u_over_gap': sorted(fast),
            'mask_lsb_is_residue_zero': hex(sum(1 << r for r in fast)),
            'zero_radius_case_count': len(zero),
            'zero_radius_example': zero[0] if zero else None,
            'denominator_zero_cases_rejected': invalid,
        })
    assert any(row['allowed_count'] for row in results) if args.limit >= 17 else True
    if args.limit >= 41:
        assert any(row['zero_radius_case_count'] for row in results)
    data = {
        'scope': 'whole_plane_critical_content_projection',
        'formula': 'F_i(r)=(2r-1)(2r(r-1)-1)-2i*r(r-1)*(r(r-1)+2)',
        'conditions': 'r!=0,1; i*i=-1; r-1-i*r!=0; F_i(r) square including zero',
        'one_sign_symmetry': 'test one i; include r and1-r',
        'complexity': 'O(p) field operations and O(p) table space for each prime; not an O(p) bit-complexity claim',
        'limit': args.limit, 'prime_count': len(results), 'certificates': results,
        'cpu_seconds': time.process_time()-started,
        'old_pattern_certificate_unchanged': True,
        'mathematical_novelty_claimed': False,
    }
    with open(args.output, 'w', encoding='utf-8') as output:
        json.dump(data, output, ensure_ascii=False, indent=2)
        output.write('\n')
    print(json.dumps({'output': args.output, 'all_methods_match': True,
                      'prime_count': len(results), 'cpu_seconds': data['cpu_seconds'],
                      'mask_53': next((c['mask_lsb_is_residue_zero'] for c in results if c['prime']==53), None),
                      'mask_61': next((c['mask_lsb_is_residue_zero'] for c in results if c['prime']==61), None)}, indent=2))


if __name__ == '__main__':
    main()
