# Rational distances to the four vertices of the unit square

Is there a point P in the plane whose distances to all four vertices of the unit square
(0,0), (1,0), (0,1), (1,1) are rational? This is an open problem. This repository contains two
computations that find no such point. Neither one solves the problem.

Site: <https://wand125.github.io/rational-distance/> (overview, the paper draft in English and Japanese,
and a work-in-progress page on extending the circles to height 10000).

Japanese version: [README.ja.md](README.ja.md).

## Results

**1. Whole plane, bounded denominator (exhaustive search, `all-plane/`).**
There is no point P in the whole plane (inside or outside the square) whose four vertex distances
are rational and whose coordinates have least common denominator at most 10,000,000.

- The search enumerates candidates through necessary arithmetic conditions, then reconstructs every
  surviving candidate with exact integer arithmetic.
- Witnesses with three rational distances are written out and re-checked with Python integers.
  Each 2,000,000-wide denominator block produced about 100 of them.
- This is a finite computation, not a proof of non-existence.

**1b. A bounded box, larger denominators (partner search, `partner-search/`).**
There is no such point in the box [−9, 10]² with common denominator at most 10¹¹.

- For each |C| the program lists the Pythagorean partners from the divisors of C², then tests pairs.
  The idea of enumerating partners comes from the GitHub repository Mulatram-UwU/rational-distance-problem
  (interior only); we extended it to signed coordinates and wrote our own Rust code.
- It agrees exactly with independent brute force on small boxes, and its interior counts agree with the
  counts published in that repository for N = 10⁴ to 10⁷.
- N = 10⁹, 10¹⁰, 10¹¹ took 11,957, 160,343 and 1,836,020 CPU seconds. Again a finite computation.

**2. Circles centred at a vertex (proofs per circle, `circles/`).**
Any such point lies on a circle |P − V| = r about some vertex V with r rational. By the symmetry of the square,
it suffices to take V = (0,0). We took every rational r ≠ 1 with max(|numerator|, denominator) ≤ 1000:
608,382 circles in all. For every one of them we prove that the circle carries no point with four rational
vertex distances. The circle r = 1 is handled by hand.

| Method | Circles excluded |
|---|---:|
| Common square class descent with p-adic local tests (p ≤ 2000) | 606,491 |
| A genus-1 quotient of rank 0, followed by an exact check of all its rational points | 588 + 7 |
| Refined 2-descent over Q(i) (primes up to 10⁶) | 1,140 + 31 |
| Mordell–Weil sieve on a genus-1 quotient | 116 + 4 |
| A genus-1 quotient of rank 0 over Q(i), followed by an exact check of its points | 5 |
| Open | 0 |

Each circle is counted once, in the first row that excludes it. The Mordell–Weil sieve excluded
151 circles in all; 31 of them are also excluded by the refined descent. The "+ 4" and the Q(i) row
are the last 9 circles to be closed (r = 607/120, 863/360, 967/360, 167/480, 727/720, 343/960, 887/960, 911/960, 383/984).

### What the circle proofs rely on

- **Descents and local tests.** A "no point" verdict is issued only where the local square class is
  determined. Undetermined discs count as "may contain a point", so only the "no point" verdicts are
  used as proofs. These verdicts were tested on circles that carry real points with rational distances to some of
  the vertices: 372 circles for the first descent, and 2,579 cases for the refined descent.
  No real point was lost. An error found on 11/08 was fixed and every affected verdict recomputed
  (see `circles/README.ja.md`).
- **Mordell–Weil sieve.**
  - Generators come from PARI `ellrank` (lower bound = upper bound) and are saturated.
  - For the last 4 circles `ellrank` proved rank 1 but returned no generator (expected heights 80–130).
    The generators came from a Heegner point (343/960), from rational points on 2-coverings of isogenous
    curves (167/480, 383/984), and from rational points on 4-coverings built by splitting a 2-covering
    quartic into two quadratics (887/960, `foursearch.gp`). The 4-covering search reproduces the known
    generators of 167/480 and 383/984.
  - Every kill was checked in three ways: positive controls with planted points, a direct cross-check
    of small multiples, and an independent second system (Sage with eclib: j-invariant, torsion,
    saturation index 1 with the height bound, and the rank).
  - The rank in that second system comes from mwrank. Where mwrank could not decide, it comes instead
    from Gross–Zagier–Kolyvagin, using the root number and a rigorously bounded L(E,1) or L′(E,1) ≠ 0.
- **Rank 0 quotients.** The 7 circles marked "+ 7" use L(E,1) ≠ 0 with an error bound (Sage
  `lseries().at1`), so the rank is 0 by Kolyvagin.
- **Rank 0 over Q(i).** The quartics of these quotients have roots in Q(i), so the curve is an elliptic
  curve over Q(i). When the Jacobian J has rank 0 over Q (and, for 3 circles, also its twist by −1),
  the rational points are found exactly from the torsion over Q(i) and the Galois action. The ranks come from
  L(J,1) ≠ 0 with an error bound.
- **One circle, r = 769/360,** relies on a single system for the rank. The sieve curve has rank 2, the
  upper bound comes only from PARI `ellrank`, and mwrank could not prove it on any curve of the isogeny class.

## Contents

- `all-plane/`: the whole-plane search.
  - It includes the earlier package (frozen versions v2–v4, their proofs of the arithmetic conditions,
    and the 2,000,000 run) and the versions v5–v8 used to reach 10⁷.
  - The progress report is `all-plane/progress_v5_2026-11-05.md`, and the integer bounds are in
    `all-plane/arithmetic_bounds_v5.md`.
  - Results are in the `sweep_*` directories (`summary.json`, `three_candidates.csv`).
  - v6 is a Rust port of v5 written with OpenAI Codex and accepted by exact comparison with v5.
- `partner-search/`: the box search (versions v1–v5, brute-force checks in `check/`, merged `summary.json` for each N; README in Japanese).
- `docs/`: the project site and a paper draft in English and Japanese (`docs/paper.html`, `docs/paper.ja.html`).
- `circles/`: the circle computation. It contains the PARI/GP and Python code, every verdict file,
  and a detailed report in Japanese (`circles/README.ja.md`) with SHA-256 prefixes of the files it cites.

## Software

- PARI/GP 2.17.2 and PARI 2.15.4 (cypari2).
- Sage 10.9 and passagemath 10.6 with eclib.
- A C++17 compiler.
- Rust 1.86.

## Notes

**What appears new, as far as our literature survey goes** (to November 2026; Section 2.6 of `docs/paper.html`):

- We found no earlier work that fixes the distance to one vertex and proves the circles empty one by one.
- We found no whole-plane search, with method and output we could check, reaching denominator 10⁷.
  Larger whole-plane bounds in forum posts could not be verified.
- The box [−9, 10]² contains the square, so result 1b covers the interior to 10¹¹. The largest interior
  search we found (Mulatram-UwU, 5×10¹⁰) is the source of the partner idea.

The individual methods are standard, and both searches are finite computations.

License: MIT (see `LICENSE`). Author: wand125.
