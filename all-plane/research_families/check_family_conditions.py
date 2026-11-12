#!/usr/bin/env python3
"""Independent finite checks of the algebraic candidate-generation conditions.

The CRT set is compared against direct enumeration of A modulo n. The positive
examples are previously generated exact THREE-distance points whose coordinate
sum/difference conditions satisfy GOOD; they are not four-distance solutions.
"""
import csv
import itertools
import json
import math
from pathlib import Path


def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def valuation(n, p):
    n = abs(n)
    if n == 0:
        return math.inf
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def good(n):
    return n > 0 and n % 2 == 1 and all(p % 8 in (1, 7) for p in factor(n))


def local_family_certificate(p):
    roots = {s: [r for r in range(p) if r*r % p == s] for s in range(p)}
    points, plus, minus = [], set(), set()
    for X in range(p):
        for Y in range(p):
            squares = [(X*X+Y*Y) % p, ((1-X)**2+Y*Y) % p,
                       (X*X+(1-Y)**2) % p, ((1-X)**2+(1-Y)**2) % p]
            if not all(roots[s] for s in squares):
                continue
            points.append([X, Y])
            for B in roots[squares[1]]:
                for C in roots[squares[2]]:
                    plus.add((B+C) % p)
                    minus.add((B-C) % p)
    H = plus | minus
    J = {r for r in range(p) if roots[(2-r*r) % p]}
    return {'prime': p, 'N_unit_coordinate_points': points,
            'N_unit_sum_residues': sorted(plus),
            'N_unit_difference_residues': sorted(minus),
            'N_divisible_integral_R_necessary_residues': sorted(J),
            'forbidden_integral_R_residues': sorted(set(range(p))-(H | J))}


def main():
    count = allowed_count = 0
    for n in range(1, 61):
        for u in range(n + 1):
            m = n // math.gcd(n, u)
            powers = [p**e for p, e in factor(m).items()]
            choices = [sorted({0, u % pe}) for pe in powers]
            if m == 1:
                generated = {0}
            else:
                weights = [(m//pe) * pow(m//pe, -1, pe) for pe in powers]
                generated = {sum(a*w for a, w in zip(residues, weights)) % m
                             for residues in itertools.product(*choices)}
            for v in range(n + 1):
                by_direct_A = any(u*(A+u) % n == 0 and v*(A+v) % n == 0
                                  for A in range(n))
                by_generated_set = (v % m) in generated
                assert by_direct_A == by_generated_set, (n, u, v)
                count += 1
                allowed_count += by_direct_A

    root = Path(__file__).resolve().parents[1]
    rows = list(csv.DictReader((root / 'near_without_modular_100000000.csv').open()))
    checked = 0
    for row in rows:
        N, X, Y = (int(row[k]) for k in ('N', 'X', 'Y'))
        # Put the missing corner at (N,N), retaining three exact distances.
        if row['missing'] == 'S':
            Y = N-Y
        values = [X*X+Y*Y, (N-X)**2+Y*Y, X*X+(N-Y)**2]
        A, B, C = [math.isqrt(s) for s in values]
        assert [A*A, B*B, C*C] == values
        assert math.gcd(math.gcd(X, Y), N) == 1
        assert N % 24 == 0 and A % 2 and B % 2 and C % 2
        assert good(abs(X-Y))
        n, u, v = N//2, (B-A)//2, (C-A)//2
        delta, S = abs(u-v), A+u+v
        assert delta > 0 and S > 0
        assert delta*S == n*abs(Y-X)
        assert good(delta // math.gcd(delta, n))
        assert good(S // math.gcd(S, n))
        for p, e in factor(delta).items():
            if p % 8 in (3, 5):
                assert e <= valuation(n, p)
            if p % 8 == 3:
                assert e == valuation(n, p)
        e2 = valuation(n, 2)
        assert valuation(delta, 2) in (0, e2)
        for z in (u, v):
            vz = valuation(z, 2)
            assert vz in (0, e2) or vz >= e2+2
        for p, e in factor(n).items():
            if p == 2 or p % 4 == 3:
                for z in (u, v):
                    vz = valuation(z, p)
                    assert vz == 0 or vz >= e
        g = math.gcd(n, math.gcd(u, v))
        n0, u0, v0 = n//g, u//g, v//g
        Q0 = u0*u0+v0*v0-n0*n0
        assert (n0*Q0) % g == 0
        if g % 2 == 0:
            assert Q0 % (2**(valuation(g, 2)+1)) == 0
        checked += 1

    # This is a nontrivial positive example for the valid general scaling rule
    # and a negative control for the unjustified replacement g | Q0.
    N, X, Y, A, B, C = 20280, 4795, 8772, 9997, 17797, 12467
    assert X*X+Y*Y == A*A
    assert (N-X)**2+Y*Y == B*B
    assert X*X+(N-Y)**2 == C*C
    n, u, v = N//2, (B-A)//2, (C-A)//2
    g = math.gcd(n, math.gcd(u, v))
    n0, u0, v0 = n//g, u//g, v//g
    Q0 = u0*u0+v0*v0-n0*n0
    assert (n0*Q0) % g == 0 and Q0 % g != 0
    local_certificates = [local_family_certificate(p) for p in (11, 19)]
    assert local_certificates[0]['forbidden_integral_R_residues'] == [0, 5, 6]
    assert local_certificates[1]['forbidden_integral_R_residues'] == [0, 3, 16]
    print(json.dumps({
        'passed': True,
        'crt_moduli': 'n=1..60',
        'crt_coordinate_triples_checked': count,
        'crt_allowed_triples': allowed_count,
        'crt_reference': 'direct enumeration of A modulo n',
        'existing_good_three_distance_examples_checked': checked,
        'scaling_negative_control': {
            'N': N, 'X': X, 'Y': Y, 'A': A, 'B': B, 'C': C,
            'g': g, 'n0': n0, 'Q0': Q0,
            'g_divides_n0_Q0': True, 'g_divides_Q0': False},
        'prime_11_and_19_curve_family_certificates': local_certificates,
        'scope': 'finite implementation checks supplement the algebraic proofs; '
                 'the positive examples are not four-distance solutions'
    }, indent=2))


if __name__ == '__main__':
    main()
