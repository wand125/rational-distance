#!/usr/bin/env python3
"""Exact finite-field certificates for primitive four-distance content.

The proof and scope are in critical_family_refinement.md. This is a necessary
condition certificate, not a search for global rational four-distance points.
All arithmetic is integral; the direct radius enumeration is independent of
the square-class reduction used to produce the implementation masks.
"""
import argparse
import json
import math
import time


def primes(limit):
    return [p for p in range(3, limit + 1, 2)
            if all(p % d for d in range(2, math.isqrt(p) + 1))]


def square_roots(p):
    result = [[] for _ in range(p)]
    for a in range(p):
        result[a*a % p].append(a)
    return result


def critical_projection(p):
    if p % 4 != 1 or p not in primes(p):
        raise ValueError('critical content requires a prime p == 1 mod4')
    roots = square_roots(p)
    iota = next(i for i in range(1, p) if i*i % p == p-1)
    nonsquare = next(i for i in range(1, p) if not roots[i])
    allowed = set()
    patterns = []
    prefix_cases = []
    # Rescaling all four radii by a nonzero field element rescales the
    # square increment by a square; 1 and nonsquare represent both classes.
    for shift in (1, nonsquare):
        for c in range(p):
            values = [c, (c+shift) % p, (c+iota*shift) % p,
                      (c+(1+iota)*shift) % p]
            if roots[values[0]] and roots[values[1]]:
                prefix_cases.append({'shift': shift, 'c': c, 'squares': values,
                                     'all_four_squares': all(roots[v] for v in values)})
            if not all(roots[v] for v in values):
                continue
            patterns.append({'shift': shift, 'c': c, 'squares': values})
            for a in roots[values[0]]:
                for b in roots[values[1]]:
                    for croot in roots[values[2]]:
                        u, v = (b-a) % p, (croot-a) % p
                        # The critical squared differences are all nonzero.
                        assert u and v and (u-v) % p
                        z = u*pow((u-v) % p, -1, p) % p
                        allowed.add(z)
                        # Swapping adjacent corners changes iota to -iota.
                        allowed.add((1-z) % p)
    return {
        'prime': p, 'iota': iota, 'nonsquare_representative': nonsquare,
        'square_class_patterns': patterns,
        'prefix_cases': prefix_cases if p in (5, 13) else [],
        'allowed_u_over_gap': sorted(allowed),
        'forbidden_u_over_gap': sorted(set(range(p))-allowed),
        'new_forbidden_beyond_gamma_units': sorted(set(range(2, p))-allowed),
        'allowed_count': len(allowed),
        'mask_lsb_is_residue_zero': hex(sum(1 << z for z in allowed)) if p < 64 else None,
    }


def direct_radius_projection(p):
    """Enumerate every signed (a,b,c) without square-class normalization.

At a critical prime, b^2-a^2 and c^2-a^2 are nonzero and their squares
sum to zero; the fourth radius square is b^2+c^2-a^2.
"""
    roots = square_roots(p)
    allowed = set()
    triples = 0
    for a in range(p):
        aa = a*a % p
        for b in range(p):
            bb = b*b % p
            first = (bb-aa) % p
            if first == 0:
                continue
            for c in range(p):
                cc = c*c % p
                second = (cc-aa) % p
                if (first*first+second*second) % p:
                    continue
                if not roots[(bb+cc-aa) % p]:
                    continue
                u, v = (b-a) % p, (c-a) % p
                assert u and v and (u-v) % p
                allowed.add(u*pow((u-v) % p, -1, p) % p)
                triples += 1
    return allowed, triples


def diagonal_gap_projection(p):
    """Recheck the previous affine H_p plus N-divisible J_p criterion."""
    roots = square_roots(p)
    affine = set()
    coordinate_count = 0
    for x in range(p):
        for y in range(p):
            values = [(x*x+y*y) % p, ((1-x)**2+y*y) % p,
                      (x*x+(1-y)**2) % p, ((1-x)**2+(1-y)**2) % p]
            if not all(roots[v] for v in values):
                continue
            coordinate_count += 1
            for b in roots[values[1]]:
                for c in roots[values[2]]:
                    affine.add((b-c) % p)
    infinity = {r for r in range(p) if roots[(2-r*r) % p]}
    forbidden = set(range(p))-affine-infinity
    return {'prime': p, 'affine_coordinate_count': coordinate_count,
            'forbidden_residues': sorted(forbidden),
            'nonzero_polynomial_coefficients': sorted({r*r % p for r in forbidden if r})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='critical_projection_certificate.json')
    args = parser.parse_args()
    started = time.process_time()
    selected = [5, 13, 17, 29, 37, 41, 53, 61]
    certs = []
    for p in selected:
        cert = critical_projection(p)
        direct, triples = direct_radius_projection(p)
        assert direct == set(cert['allowed_u_over_gap'])
        cert['independent_radius_projection_matches'] = True
        cert['direct_signed_radius_triple_count'] = triples
        certs.append(cert)
    assert certs[0]['allowed_count'] == certs[1]['allowed_count'] == 0
    assert certs[2]['allowed_u_over_gap'] == list(range(2, 17, 2))
    critical_survey = [critical_projection(p) for p in primes(257) if p % 4 == 1]
    assert [c['prime'] for c in critical_survey if not c['allowed_count']] == [5, 13]
    previous_criterion = [diagonal_gap_projection(p) for p in primes(251) if p % 4 == 3]
    assert {c['prime']: c['nonzero_polynomial_coefficients'] for c in previous_criterion
            if c['nonzero_polynomial_coefficients']} == {11: [3], 19: [9]}
    data = {
        'scope': 'whole_plane_primitive_four_distance_necessary_conditions',
        'source': 'https://zenodo.org/records/23638745',
        'source_section': '08_critical_content.tex; critical pattern, thm:crit-basic, thm:crit-17',
        'variables': 'n=N/2,u=(B-A)/2,v=(C-A)/2,h=u-v',
        'applicability': 'p^(2s)||n and p^s||h with s>=1; then z=(u/p^s)/(h/p^s) modp must be allowed',
        'certificates': certs,
        'critical_survey_through_257': [
            {'prime': c['prime'], 'allowed_count': c['allowed_count']}
            for c in critical_survey],
        'old_diagonal_criterion_survey_through_251': previous_criterion,
        'cpu_seconds': time.process_time()-started,
        'global_rational_solution_found': False,
        'mathematical_novelty_claimed': False,
    }
    with open(args.output, 'w', encoding='utf-8') as destination:
        json.dump(data, destination, ensure_ascii=False, indent=2)
        destination.write('\n')
    print(json.dumps({'output': args.output,
                      'critical_masks': {str(c['prime']): c['mask_lsb_is_residue_zero'] for c in certs},
                      'all_independent_projections_match': True,
                      'cpu_seconds': data['cpu_seconds']}, indent=2))


if __name__ == '__main__':
    main()
