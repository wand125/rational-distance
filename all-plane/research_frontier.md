# Public computational reports: scope and audit

This project note records primary sources checked during the November 2026 research session. It is not a claim to identify a world record. A numerical maximum is meaningful only with its search domain, normalization, coverage proof, and execution evidence.

## Findings at a glance

| Source | What the author reports | Scope and interpretation | Verification here |
|---|---|---|---|
| Ed Pegg, 2 June 2016, [Wolfram Community](https://community.wolfram.com/t/rational-distance-problem/8686) | Analysis of primitive Heronian triangles with longest side at most 1,410,000; 2,877 canonical three-distance points | A Heronian-triangle construction, including exterior examples. The bound is not the common denominator of a square point. The same post says data through 6,000,000 are available and requests volunteers; it does not claim that larger search was completed. | Original author's post read in full. No complete search program or run log was supplied in that post. |
| William McCloskey, [arXiv:1904.12097v3](https://arxiv.org/html/1904.12097v3), 5 June 2020, section 3 | Solutions to the distance equations exist modulo powers of 3 through 3^100 | A local-solubility computation. It does not exclude rational points with denominators below 3^100, nor establish a rational solution. | Full primary text inspected. |
| Yuan Si, [Zenodo 19613981](https://zenodo.org/records/19613981), April 2026, v2 | No interior solutions with common denominator at most 100,000; no solutions in a Pythagorean-leg search through 5,000,000 | Interior only. The leg cutoff is not an exhaustive denominator cutoff of twice that size. | Original ZIP downloaded, checksum matched, report and code inspected. Two issues below prevent treating its full advertised coverage as independently verified. |
| Yimin Rong, [Math StackExchange](https://math.stackexchange.com/questions/1379191/simplified-rational-distance-problem), July 2015 | A search with integer side L < 2^31 and near-misses a+b=L±1 | The question imposes a+b=L on two Pythagorean triples. This is a point on an edge, a much narrower subproblem. | Original post read. It is not an all-plane four-distance denominator record. |
| Bremner and Guy, [Nu-configurations in tiling the square](https://doi.org/10.1090/S0025-5718-1992-1134716-2), 1992 | Rational points on 772 curves of an elliptic pencil | A different four-triangle tiling configuration. It is not a report of a uniform search for a point at rational distances from all four square vertices. | [Primary author-institution record](https://asu.elsevierpure.com/en/publications/nu-configurations-in-tiling-the-square/) inspected; full publisher PDF unavailable through the retrieval tool. |

No uniform all-plane common-denominator frontier was verified from these public sources. This statement describes the evidence inspected; it does not establish that no larger computation exists.

## Why the cutoffs cannot be compared directly

For a rational square point write x=X/N and y=Y/N with gcd(X,Y,N)=1. An exhaustive all-plane N cutoff permits X and Y outside [0,N]. A denominator sweep restricted to 0<Y<X<N covers the interior only.

An interior leg-index search often uses positive p,r,s,t with p+s=r+t=q and p²+r², s²+r², p²+t², s²+t² square. The condition max(p,r,s,t)≤B allows q as large as 2B, but does not cover every interior point with q≤2B. It does cover all interior points with q≤B if the leg index itself is complete and the search has actually completed. Neither implication converts such a search into an all-plane sweep.

A bound on a primitive Heronian triangle's longest side is a bound on a different normalized object. A formula relating that normalization to N would be needed before converting Pegg's number into a square-denominator statement. The post does not supply such a conversion.

Finally, finding points modulo 3^k is a necessary-condition exercise. Rational points need compatible conditions at every relevant place and must satisfy the exact equations. McCloskey's 3^100 computation is especially useful as a warning that deeper 3-adic checking does not automatically finish the problem.

## Audit of the April 2026 Zenodo archive

The exact public archive is [reports_v1_2.zip](https://zenodo.org/records/19613981/files/reports_v1_2.zip?download=1), published MD5 `ddc41e2c31d8a5a80b45c2179662b43f`. The downloaded bytes match. The archive states CC-BY 4.0 for text and MIT for code. Its provided output files concern algebra/CAS work; no completed search-run output log was found.

The following are direct checks of the supplied artifacts, not claims about the author's intentions or about the existence of an actual four-distance point.

### 1. Pythagorean generator misses triples within its advertised leg bound

In `scripts/search/search_5M.py`, line 17 sets B=5,000,000; lines 23–25 bound the Euclidean parameter by m≤floor(sqrt(B))+1=2237. The generator then uses a=m²−n², b=2mn and integer multiples. `two_triple_search.py` uses the same bound pattern.

Take m=2238, n=925. These parameters are coprime and have opposite parity. They give the primitive triple

\[
(a,b,c)=(4,153,019,\ 4,140,300,\ 5,864,269).
\]

Both legs are below five million, but m exceeds the loop cap. Uniqueness of the primitive Euclidean parametrization means an alternative multiplier cannot recover this missing triple. The index therefore does not contain every leg pair promised by the description.

This counterexample refutes the completeness of the published leg generator. It is not itself a point at four rational distances and does not demonstrate that an actual four-distance solution was missed.

### 2. Denominator-sweep counter and report do not agree under the natural field mapping

Report 2, section 2.4, Run 1 associates `direct_search.py` with q≤100,000, and gives 1,476,656,036 as its combinations count. The supplied script sets `Q_MAX=100000` at line 63. Its `total_tested` counter increments at line 139 after the mod-3 and coprimality filters, before the mod-5/7 filter and any square tests. It prints that counter as the total number of tested triples at line 191.

That counter can be calculated without running the expensive search. Let M=floor(Q/3). For q=3m there are 2m(m−1) unordered positive pairs p,r<q with exactly one coordinate divisible by 3. Applying Möbius inversion to gcd(p,r,q)=1 and summing over m gives

\[
C(Q)=\sum_{\substack{d\le M\\3\nmid d}}\mu(d)\frac{2}{3}
\left(\left\lfloor\frac{M}{d}\right\rfloor^3-
\left\lfloor\frac{M}{d}\right\rfloor\right).
\]

At Q=100,000 this is

\[
C(100000)=21,330,687,731,144,
\]

which differs from 1,476,656,036. The formula was checked against literal loop enumeration at Q=30 and Q=300. This comparison is specifically to `total_tested`, not to a count after a later filter.

If the report's combinations field is intended to mean some different stage, the archived report does not provide the run log or mapping needed to verify that interpretation. Therefore the supplied report and code do not establish a reproducible completed q≤100,000 run on their own. This conclusion is narrower than saying that the mathematical absence claim is false.

### Reproduction and preserved sources

`audit_zenodo_search_claims.py` independently reproduces the counter formula, two small brute-force comparisons, and the omitted primitive triple. Its output is `public_search_audit.json`. It does not run the third-party searches or alter any project search implementation or benchmark log.

The original downloaded archive and extracted copies are kept outside the deliverable project under `../rational_distance_sources/`. Relevant paths are:

- `reports_v1_2.zip`
- `zenodo_v2/rational_distance_research_archive.pdf`
- `zenodo_v2/report.txt` (PDF text extraction, Run 1 at lines 383–393)
- `zenodo_v2/scripts/search/direct_search.py`
- `zenodo_v2/scripts/search/search_5M.py`
- `zenodo_v2/scripts/search/two_triple_search.py`

## Status and nearby results

Corvaja, Turchet, and Zannier's peer-reviewed 2025 article, [Rational distances from given rational points in the plane](https://link.springer.com/article/10.1007/s10711-025-01019-0), explicitly treats the unit-square case as having no known rational-distance point. Its algebraic-geometric results do not supply a finite denominator bound or a nonexistence proof for every point. Recent repository items must be evaluated separately rather than treated as a change to the accepted status merely because they have a DOI.

[Elkies and Ulas's primary MathOverflow answers](https://mathoverflow.net/questions/460755/a-rational-distance-problem-with-possibly-multiple-solutions), December 2023, describe infinite rational-distance constructions for other quadrilaterals and a dense set of rational rectangle aspect ratios. They do not give a solution or a computational exclusion frontier for the unit square. In particular, a dense set of successful rectangle aspect ratios need not include aspect ratio 1.

For reporting this project's own computation, retain its exact normalization, complete parameter range, primitive-point convention, symmetry reduction, proof of coverage, source version, and execution logs. Compare only computations that use a provably comparable domain and cutoff.
