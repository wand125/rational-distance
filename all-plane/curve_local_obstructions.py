#!/usr/bin/env python3
"""Certificates for whole rational ellipse/hyperbola families.

For prime ell=3 mod4, two branches cover every primitive configuration:
  ell does not divide N: exact affine four-distance projection;
  ell divides N and R is ell-integral: 2-R^2 must be a square modulo ell.
R is either the sum or the difference of distances to opposite vertices.
The second branch is a necessary condition, not a sufficiency claim.
The complete proof is in research_families/family_exclusions.md.
"""
import argparse
import json
import math
from local_obstructions import certificate


def curve_certificate(ell):
    if ell % 4 != 3:
        raise ValueError('the proof requires a prime congruent to3 mod4')
    data = certificate(ell)
    finite = {(u-v) % ell for u,v in data['projected_difference_pairs']}
    squares = {a*a % ell for a in range(ell)}
    infinity_necessary = {r for r in range(ell) if (2-r*r) % ell in squares}
    forbidden = set(range(ell))-finite-infinity_necessary
    coefficients = sorted({r*r % ell for r in forbidden if r})
    return {'prime':ell,
            'affine_coordinate_classes':data['coordinate_classes'],
            'affine_allowed_sum_or_difference':sorted(finite),
            'N_divisible_branch_necessary_residues':sorted(infinity_necessary),
            'forbidden_residues_when_reduced_denominator_is_a_unit':sorted(forbidden),
            'nonzero_polynomial_coefficients_c':coefficients,
            'polynomial_condition':'for every listed c, ell does not divide h*h-c*k*k, for reduced R=h/k',
            'zero_residue_also_excluded':0 in forbidden}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primes',nargs='+',type=int,default=[11,19])
    args = parser.parse_args()
    certs = [curve_certificate(p) for p in args.primes]
    examples = []
    available = {c['prime']:c for c in certs}
    for kind,h,k,p in [('difference',1,2,11),('difference',7,8,11),
                       ('sum',7,3,11),('sum',7,4,19),('difference',1,6,19)]:
        if p not in available:
            continue
        assert math.gcd(h,k) == 1 and k % p
        residue = h*pow(k,-1,p) % p
        assert residue in available[p]['forbidden_residues_when_reduced_denominator_is_a_unit']
        nondegenerate = h*h > 2*k*k if kind=='sum' else 0 < h*h < 2*k*k
        assert nondegenerate
        examples.append({'kind':kind,'numerator':h,'denominator':k,
                         'obstructing_prime':p,'residue':residue,
                         'nondegenerate_real_curve':nondegenerate})
    print(json.dumps({'certificates':certs,'excluded_curve_examples':examples},indent=2))


if __name__ == '__main__':
    main()
