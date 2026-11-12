#!/usr/bin/env python3
"""Finite certificates for the unit-side four-distance equations modulo primes.

If prime ell does not divide N, divide coordinates and distances by N in F_ell.
The returned coordinate and (p/N,q/N) lists are exhaustive necessary classes.
Every square root, including zero, is enumerated. There are no heuristics.
"""
import argparse
import json
import math


def certificate(ell):
    if ell < 3 or any(ell % d == 0 for d in range(2, math.isqrt(ell)+1)):
        raise ValueError('an odd prime is required')
    roots = [[] for _ in range(ell)]
    for a in range(ell):
        roots[a*a % ell].append(a)
    points = []
    pairs = set()
    for x in range(ell):
        for y in range(ell):
            values = [x*x+y*y, (1-x)**2+y*y,
                      x*x+(1-y)**2, (1-x)**2+(1-y)**2]
            rs = [roots[z % ell] for z in values]
            if not all(rs):
                continue
            points.append([x,y])
            for a in rs[0]:
                for b in rs[1]:
                    for c in rs[2]:
                        pairs.add(((b-a) % ell, (c-a) % ell))
    allowed_rows = sorted({p for p,q in pairs})
    midpoint = (ell+1)//2
    edge_midpoint_classes = sorted([[0,midpoint], [1,midpoint],
                                    [midpoint,0], [midpoint,1]])
    return {'prime':ell, 'condition':'prime does not divide N',
            'coordinate_classes':points,
            'exactly_four_edge_midpoint_classes':points == edge_midpoint_classes,
            'projected_difference_pairs':[list(t) for t in sorted(pairs)],
            'allowed_difference_rows':allowed_rows,
            'forbidden_difference_rows':[i for i in range(ell) if i not in allowed_rows],
            'projection_count':len(pairs), 'all_pair_count':ell*ell}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primes', nargs='+', type=int, default=[3,5,7,11,13,17,19])
    args = parser.parse_args()
    data = [certificate(p) for p in args.primes]
    print(json.dumps({'method':'exhaustive finite-field coordinates and signed distance roots',
                      'certificates':data}, indent=2))


if __name__ == '__main__':
    main()
