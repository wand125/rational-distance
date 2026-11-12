#!/usr/bin/env python3
"""Linear-residue construction of all critical next-digit exceptions <=281."""
import argparse
import json
import time
from pathlib import Path
from critical_lift import quartic, derivative, projection


def gaussian_mul(x, y):
    return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def gaussian_pow(x, exponent):
    answer = (1, 0)
    for _ in range(exponent):
        answer = gaussian_mul(answer, x)
    return answer


def quartic_discriminant():
    # Universal quartic discriminant formula, evaluated in Z[i], with no floats.
    a, b, c, d, e = (0, -2), (4, 4), (-6, -6), (0, 4), (1, 0)
    terms = [
        (256, (3,0,0,0,3)), (-192, (2,1,0,1,2)),
        (-128, (2,0,2,0,2)), (144, (2,0,1,2,1)),
        (-27, (2,0,0,4,0)), (144, (1,2,1,0,2)),
        (-6, (1,2,0,2,1)), (-80, (1,1,2,1,1)),
        (18, (1,1,1,3,0)), (16, (1,0,4,0,1)),
        (-4, (1,0,3,2,0)), (-27, (0,4,0,0,2)),
        (18, (0,3,1,1,1)), (-4, (0,3,0,3,0)),
        (-4, (0,2,3,0,1)), (1, (0,2,2,2,0)),
    ]
    answer = (0, 0)
    for factor, exponents in terms:
        term = (factor, 0)
        for coefficient, exponent in zip((a,b,c,d,e), exponents):
            term = gaussian_mul(term, gaussian_pow(coefficient, exponent))
        answer = (answer[0]+term[0], answer[1]+term[1])
    assert answer == (16384, 0)
    return answer


def sparse_data(p):
    squares = {a*a % p for a in range(p)}
    iotas = [a for a in range(p) if a*a % p == p-1]
    old, unit, exceptions = set(), set(), []
    for r in range(2, p):
        viable = [i for i in iotas if (r-1-i*r) % p
                  and quartic(p, i, r) in squares]
        if not viable:
            continue
        old.add(r)
        if any(quartic(p, i, r) != 0 for i in viable):
            unit.add(r)
            continue
        roots = []
        for i in viable:
            fp = derivative(p, i, r)
            assert fp != 0
            i2 = i+p*((-((i*i+1)//p)*pow(2*i, -1, p)) % p)
            assert (i2*i2+1) % (p*p) == 0
            digit = (-quartic(p*p, i2, r)//p*pow(fp, -1, p)) % p
            root = r+p*digit
            assert quartic(p*p, i2, root) == 0
            roots.append(root)
        exceptions.append({'residue_mod_p': r, 'roots_mod_p2': sorted(set(roots))})
    return old, unit, exceptions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(Path(__file__).with_name('critical_sparse_certificate.json')))
    args = parser.parse_args()
    start = time.process_time()
    discriminant = quartic_discriminant()
    certificates = []
    for p in range(5, 282, 4):
        if any(p % d == 0 for d in range(2, int(p**0.5)+1)):
            continue
        old, unit, exceptions = sparse_data(p)
        branch_count = sum(len(row['roots_mod_p2']) for row in exceptions)
        assert branch_count <= 8
        new_count = len(unit)*p+branch_count
        independently_expanded = bool(exceptions)
        if exceptions:
            full, _, _ = projection(p, 2)
            expected = {r for r in range(p*p) if r % p in unit}
            expected.update(r for row in exceptions for r in row['roots_mod_p2'])
            assert full == expected
            assert len(full) == new_count
        certificates.append({'prime': p, 'old_residue_count': len(old),
                             'old_lift_count': len(old)*p, 'new_lift_count': new_count,
                             'extra_rejected': len(old)*p-new_count,
                             'exceptions': exceptions,
                             'full_quartic_p2_projection_matches': independently_expanded,
                             'minimum_primitive_N': 24*p*p})
    active = [row['prime'] for row in certificates if row['exceptions']]
    assert active == [41, 113, 137, 257]
    data = {'limit': 281, 'prime_count': len(certificates),
            'quartic_discriminant_gaussian_integer': list(discriminant),
            'construction_complexity': 'O(p) residue/field operations per prime; inverse computations have their usual logarithmic cost',
            'active_next_digit_primes': active, 'certificates': certificates,
            'cpu_seconds': time.process_time()-start,
            'mathematical_novelty_claimed': False}
    with open(args.output, 'w', encoding='utf-8') as out:
        json.dump(data, out, indent=2)
        out.write('\n')
    print(json.dumps({'output': args.output, 'prime_count': len(certificates),
                      'active_next_digit_primes': active,
                      'cpu_seconds': data['cpu_seconds'], 'all_assertions_passed': True}, indent=2))


if __name__ == '__main__':
    main()
