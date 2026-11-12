# Task: faithful Rust port of the whole-plane search v5 (`rsq6`)

Working directory: the parent of this file (`rational_square_search_v4/`). Write ONLY inside `rust_v6/`.
Do not modify any other file. In particular the C++ sources `all_plane_gap_v5.cpp`, `wide_arith_v5.hpp`,
`all_plane_gap_v4.cpp`, `all_plane_gap_v3.cpp`, `all_plane_sieved.cpp`, `research_v4_local/lifted_projection.hpp`,
`research_v4_math/critical_extension.hpp` are frozen references. Read them completely; they are the specification.

## What to build

A cargo crate `rust_v6/` (edition 2021, **no external dependencies**, must compile with rustc 1.85 and newer,
`[profile.release] opt-level=3, lto="fat", codegen-units=1, panic="abort"`) producing binary `rsq6`.

It must reproduce exactly the DEFAULT configuration of the C++ binary `all_plane_gap_v5_local`
(gamma stage `late`, skip_redundant19=true, additional filters on: fast early five-adic mark, prime-power
projections 49/121/169/125, critical extension, endpoint pruning), and the `--no-new-filters` configuration
(NoAdditionalSieve: no additional sieve, endpoints NOT trimmed). Other gamma stages and other flags are not needed.

Pipeline per side n = N/2 (n = 12k, N = 24k, N from min_denominator..bound): exactly `run_v5` in
`all_plane_gap_v5.cpp`: permitted_valuations, CriticalContentRows, additional.mark, GammaPowerPlan,
CriticalProjectionPlan, gapwheel.generate over h, endpoint test, permitted_difference, critical_rows,
the CRT/wheel generators (`generate_field_gap_crt_range` when endpoints are trimmed, `generate_field_gap_crt`
otherwise), then per pair: permitted, local_without19, additional.accepts, permitted_scale, late gamma test
(gamma_value), critical_projection, then the exact reconstruction `evaluate_v5` from `wide_arith_v5.hpp`
(port the 256-bit arithmetic exactly: U256 limbs, mul128, add/sub/shr, isqrt256 restoring method,
residue prefilters 64/63/65/11 via one reduction mod 45045, i128 coordinates), then the accept step
(fourth distance D^2 via U256 square test). Port every table builder (ArithmeticTables with spf/bad_part/
CRTPlan idempotents, LocalSieve PrimeProjection<P> for P in 5..43, GapWheel, ProgressionWheel<19>,
LiftedProjection<M,P>, FifthFastMarker, FifthPowerCondition, CriticalExtension with its lift tables,
critical_mask) exactly. Find each definition in the C++ files (some live in all_plane_gap_v3.cpp and
all_plane_sieved.cpp; read them fully, including gamma_value, permitted_scale, permitted_difference).
Use u32/u64/i128/u128 exactly where the C++ does; never use floating point in any decision.

## Parallelism and memory

One process, `--threads T` (default 1) worker threads via `std::thread::scope`. Share the read-only tables
(ArithmeticTables up to bound/2, wheels, projections) across threads (build once). Distribute the sides n
dynamically (atomic counter over k; give out small chunks of consecutive k, larger k cost more). Each thread keeps
its own per-n buffers and its own counters. Output must be deterministic and independent of T:
collect per-n results (parameter rows, three-distance rows, solutions, counters), and emit them in increasing n
order (e.g. a reorder buffer keyed by k, flushed in order). Counters are summed.

## CLI and output (must match the C++ v5 format)

`rsq6 [B<=10000000] [--min-denominator L] [--threads T] [--audit-three CSV] [--audit-parameters CSV] [--no-new-filters] [--validate]`

- Bounds checks as in v5 (1..10000000, minimum <= bound).
- `--audit-parameters`: header `N,p,q` then one row per pair reaching evaluate (written where C++ calls parameter()).
- `--audit-three`: header `N,p,q,X,Y,A` then one row per accepted three-distance point (decimal i128).
- Rows must be byte-identical to the C++ files for the same arguments (same order: increasing n; within n,
  the C++ generation order: h in gapwheel order, residues in CRT order, u increasing within a progression).
- Solutions: print `{"solution":{"N":..,"X":..,"Y":..,"A":..,"D":..}}` lines on stdout as C++ does (order by n).
- Final stdout line: one JSON object with EXACTLY the keys and integer/boolean/string values of the C++ v5
  output, in the same key order, except: `"generator":"gap_CRT_v6_rust"`, `"reconstruction":"exact_u256_v6_rust"`,
  an added key `"threads":T`, and the timing values. All counters must be equal to the C++ ones.
- `--validate`: port the C++ checks that do not need the frozen C++ evaluate: validate_gap, validate_field_wheel,
  validate_v4_generators, AdditionalSieve::verify (LiftedProjection::verify, FifthFastMarker::verify), the
  PrimeProjection::verify of LocalSieve, validate_wide_arithmetic, and part (a) brute force N=1..12 and the scaled-seed
  and extreme-parameter parts of validate_v5_reconstruction (seeds from the Rust evaluate itself at N<=12 plus
  (8400,6720,5040)). Panics/aborts on mismatch; prints a pass line per check on stderr.

## Acceptance tests you must run (on this machine: each command must finish in under 60 seconds,
## at most 2 processes at the same time, never a bound above 200000)

1. `cargo build --release` (offline is fine; there are no dependencies).
2. `./rust_v6/target/release/rsq6 1000 --validate` passes.
3. For each of these argument sets, run the C++ reference `./all_plane_gap_v5_local` and `rsq6` with
   `--audit-three` and `--audit-parameters` to files under `rust_v6/accept/` and require: both CSV files
   byte-identical (`cmp`), and all non-timing JSON values equal except generator/reconstruction/threads:
   - `100000`
   - `200000 --min-denominator 150001`
   - `100000 --no-new-filters`
   - rsq6 with `--threads 1` and `--threads 4` for `200000` (compare to C++ `200000`).
   Write a small script `rust_v6/accept.sh` that performs these comparisons and prints PASS/FAIL per case;
   run it and keep its output in `rust_v6/accept/ACCEPT_LOG.txt`.
4. Report the measured `seconds_search` of C++ and rsq6 (threads 1) for bound 200000.

Write a short `rust_v6/README.md`: build, CLI, how the port maps to the C++ functions, and the acceptance results.
If something cannot be made identical, stop and explain precisely in README.md instead of approximating.
