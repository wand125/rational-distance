#!/usr/bin/env python3
"""Small exact certificate for the infinite family surviving v4 minus 41^2.

This checks finite congruence witnesses and periods. It never constructs the
enormous integers with exponent 1+3603600k, and does not run a global search.
"""
import argparse
import json
from math import gcd, lcm
from pathlib import Path


def radius_witnesses(modulus, n, u, v):
    du = u*pow(n, -1, modulus) % modulus
    dv = v*pow(n, -1, modulus) % modulus
    inv2 = pow(2, -1, modulus)
    square_roots = {}
    for a in range(modulus):
        square_roots.setdefault(a*a % modulus, []).append(a)
    witnesses = []
    # Independently reconstruct coordinates from the first signed radius.
    for a in range(modulus):
        x = (1-du*du-2*a*du)*inv2 % modulus
        y = (1-dv*dv-2*a*dv)*inv2 % modulus
        d2 = ((1-x)**2+(1-y)**2) % modulus
        if (x*x+y*y-a*a) % modulus or d2 not in square_roots:
            continue
        b, c, d = (a+du) % modulus, (a+dv) % modulus, square_roots[d2][0]
        radicands = (x*x+y*y, (1-x)**2+y*y,
                     x*x+(1-y)**2, (1-x)**2+(1-y)**2)
        assert all((radius*radius-radicand) % modulus == 0
                   for radius, radicand in zip((a,b,c,d), radicands))
        assert (b-a-du) % modulus == (c-a-dv) % modulus == 0
        witnesses.append({'x': x, 'y': y, 'signed_radii': [a,b,c,d]})
    return du, dv, witnesses


def gamma(n, t):
    g = gcd(n, t)
    return gcd(g, n//g)


def quartic(modulus, iota, r):
    t = r*(r-1) % modulus
    return ((2*r-1)*(2*t-1)-2*iota*t*(t+2)) % modulus


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(Path(__file__).with_name('critical41_v4_family_certificate.json')))
    args = parser.parse_args()
    n, u, v = 20172, 15744, 10865
    h = u-v
    assert 0 < v < u < n
    assert gamma(n, u) == gamma(n, v) == gamma(n, h) == 41
    assert gcd(n, gcd(u, v)) == 41
    U, V, H = u//41, v//41, h//41
    assert (U, V, H) == (384, 265, 119)
    period = 55440
    checks = []
    for modulus in (49, 121, 169, 125):
        assert gcd(n, modulus) == 1
        du, dv, witnesses = radius_witnesses(modulus, n, u, v)
        assert witnesses
        order = next(k for k in range(1, modulus+1) if pow(41, k, modulus) == 1)
        period = lcm(period, order)
        checks.append({'modulus': modulus, 'normalized_u': du, 'normalized_v': dv,
                       'radius_witness_count': len(witnesses),
                       'witness': witnesses[0], 'multiplicative_order_of_41': order})
    assert period == 3603600 and period % 55440 == 0
    old_local_primes = (5,7,11,13,17,19,23,29,31,37,43)
    for modulus in (*old_local_primes, 49, 121, 169, 125):
        assert pow(41, period, modulus) == 1
    fifth = (u*pow(n, -1, 5) % 5, v*pow(n, -1, 5) % 5)
    assert fifth == (2, 0)
    r = U*pow(H, -1, 1681) % 1681
    assert r == 314 and r % 41 == 27
    first, second = quartic(1681, 1303, r), quartic(1681, 378, r)
    assert first == 533 == 41*13 and second % 41 == 11
    assert pow(11, 20, 41) == 40
    data = {
        'scope': 'integer parameter family satisfying the implemented v4 arithmetic predicates except the critical 41^2 predicate',
        'family': {'s': '1+3603600k, k>=0', 'n': '12*41^(2s)',
                   'u': '384*41^s', 'v': '265*41^s'},
        'base': {'n': n, 'u': u, 'v': v},
        'old_v3_period': 55440, 'new_period': period,
        'unit_projection_checks': checks,
        'fifth_condition_normalized_residues': list(fifth),
        'critical_41_square_obstruction': {'r_mod_1681': r, 'required_r_mod_1681': 355,
                                          'F_1303_mod_1681': first,
                                          'F_378_mod_41': second % 41},
        'not_a_machine_bound_extension': True,
        'not_a_family_of_reconstructed_points': True,
        'mathematical_novelty_claimed': False,
    }
    with open(args.output, 'w', encoding='utf-8') as out:
        json.dump(data, out, indent=2)
        out.write('\n')
    print(json.dumps({'output': args.output, 'period': period,
                      'all_four_unit_projections_have_witnesses': True,
                      'all_assertions_passed': True}, indent=2))


if __name__ == '__main__':
    main()
