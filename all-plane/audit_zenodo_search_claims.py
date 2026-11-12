#!/usr/bin/env python3
"""Small independent audit of two claims in Zenodo record 19613981.

This does not execute or copy the third-party search, and is not a search for
four-distance points. It checks a loop counter and a Pythagorean-index bound.
"""

import argparse
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path


def mobius_values(limit):
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes = []
    composite = bytearray(limit + 1)
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for p in primes:
            if n * p > limit:
                break
            composite[n * p] = 1
            if n % p == 0:
                mu[n * p] = 0
                break
            mu[n * p] = -mu[n]
    return mu


def required_counter(q_max):
    """Count 0 < r < p < q <= q_max after mod-3 and gcd filters.

    For q=3m, there are 2m(m-1) unordered pairs with exactly one coordinate
    divisible by 3. Mobius inversion removes common divisors of p,r,q;
    such a divisor cannot be divisible by 3. Summing over m gives this formula.
    """
    limit = q_max // 3
    mu = mobius_values(limit)
    total = 0
    for d in range(1, limit + 1):
        if d % 3 and mu[d]:
            k = limit // d
            total += mu[d] * (2 * (k**3 - k) // 3)
    return total


def direct_counter(q_max):
    return sum(
        1
        for q in range(3, q_max + 1, 3)
        for p in range(1, q)
        for r in range(1, p)
        if (p % 3 == 0) != (r % 3 == 0) and gcd(gcd(p, r), q) == 1
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    small_checks = []
    for bound in (30, 300):
        formula, direct = required_counter(bound), direct_counter(bound)
        assert formula == direct
        small_checks.append({"q_max": bound, "formula": formula, "direct": direct})

    bound, m, n = 5_000_000, 2238, 925
    a, b, c = m*m - n*n, 2*m*n, m*m + n*n
    cap = isqrt(bound) + 1
    assert gcd(m, n) == 1 and (m-n) % 2 == 1
    assert gcd(gcd(a, b), c) == 1
    assert 0 < a <= bound and 0 < b <= bound
    assert a*a + b*b == c*c and m > cap

    result = {
        "source_record": "https://zenodo.org/records/19613981",
        "source_zip": "https://zenodo.org/records/19613981/files/reports_v1_2.zip?download=1",
        "source_zip_published_md5": "ddc41e2c31d8a5a80b45c2179662b43f",
        "denominator_counter": {
            "q_max": 100_000,
            "stage": "total_tested increment: after mod-3 and gcd, before mod-5/7 and square tests",
            "required_total_tested": required_counter(100_000),
            "report_combinations_tested": 1_476_656_036,
            "numbers_agree": required_counter(100_000) == 1_476_656_036,
            "caveat": "Comparison identifies an inconsistency if the report's combinations field is the supplied script's total_tested field. No search-run log was found in the ZIP.",
            "formula_brute_checks": small_checks,
        },
        "leg_index_counterexample": {
            "leg_bound": bound,
            "published_m_cap": cap,
            "primitive_euclid_parameters": [m, n],
            "primitive_triple": [a, b, c],
            "both_legs_within_bound": True,
            "m_beyond_published_cap": True,
            "conclusion": "The published generator does not index every Pythagorean triple with both legs within the stated bound.",
            "caveat": "This is not a four-distance point and does not show that the search missed an actual solution.",
        },
    }
    downloaded = Path(__file__).resolve().parent.parent / "rational_distance_sources" / "reports_v1_2.zip"
    if downloaded.exists():
        result["source_zip_observed_md5"] = hashlib.md5(downloaded.read_bytes()).hexdigest()
        assert result["source_zip_observed_md5"] == result["source_zip_published_md5"]
    output = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
