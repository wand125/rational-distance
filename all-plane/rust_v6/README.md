# rsq6

Dependency-free Rust port of the frozen whole-plane C++ v5 search. Edition 2021,
minimum Rust 1.85. Release settings are `opt-level=3`, fat LTO, one codegen unit,
and abort on panic.

## Build and run

From the repository's parent-of-`rust_v6` working directory:

```sh
(cd rust_v6 && cargo build --release --offline -j 1)
./rust_v6/target/release/rsq6 1000 --validate
./rust_v6/target/release/rsq6 200000 --threads 4 \
  --audit-three rust_v6/accept/example.three.csv \
  --audit-parameters rust_v6/accept/example.parameters.csv
```

CLI:

```text
rsq6 [B<=10000000] [--min-denominator L] [--threads T]
     [--audit-three CSV] [--audit-parameters CSV]
     [--no-new-filters] [--validate]
```

Defaults: B=10000, L=1, T=1. Bounds must be 1..10000000, L<=B, and T>0.
The search uses late gamma filtering and skips the redundant mod-19 test.
The default includes all additional filters; `--no-new-filters` disables their
marking/projections/critical extension and restores untrimmed endpoints.

## Mapping to the frozen sources

| Rust source | C++ definitions ported |
| --- | --- |
| `src/main.rs` | `run_v5<late,true>`, counters, parameter/point callbacks, fourth-distance accept, CLI and JSON |
| `src/arithmetic.rs` | `ArithmeticTables`, `CRTPlan`, valuation/difference/scale predicates, gamma plan/value, critical rows/masks/projection |
| `src/generator.rs` | `GapWheel`, `ProgressionWheel<19>`, u-first/gap/field/range CRT generators; gamma labels/rules/lifting for validation |
| `src/sieve.rs` | `PrimeProjection`, `LocalSieve`, `LiftedProjection`, fifth-adic predicates/fast marker, `CriticalExtension` and every lift constant |
| `src/wide.rs` | `wide_arith_v5.hpp`: U256 limbs, exact operations, restoring square root, mod-45045 residue filters, `evaluate_v5` |
| `src/validation.rs` | gap/wheel/v4 checks, independent projection checks, wide arithmetic checks, reconstruction brute force and scaled/extreme cases |

Arithmetic decisions use integers only, with the reference's u32/u64/i128/u128
widths; Rust `usize` is used for indexing. Tables are built once and borrowed by
scoped worker threads. An atomic counter assigns two consecutive k values per
chunk. Workers reuse their own marking/active-state buffers and return per-side
counters and rows. A bounded scheduling window, channel, and ordered buffer emit
results in increasing n while retaining the original h/CRT/progression order.
`None` additional state implements the reference's `NoAdditionalSieve` behavior.

`--validate` includes all required checks that do not depend on the frozen C++
`evaluate`. Scaled seeds come from Rust reconstruction at N<=12 plus
(8400,6720,5040). The mandated arithmetic-only scaling/extreme checks reach
N=10000000; no search was run with a bound above 200000.

## Acceptance results

Release builds passed with the installed rustc 1.99.0 and 1.86.0 compilers;
both binaries passed `1000 --validate`. Rust 1.85 is the declared minimum and
the source uses APIs available by that version; 1.85 itself is not installed.
Validation includes 724800 wide-arithmetic checks, 28 brute-force points, and
116 scaled seed families. CLI empty-range/rounding/error cases also passed.

Run `./rust_v6/accept.sh` to repeat every reference comparison. It runs one child
at a time and enforces a 54-second deadline for the whole invocation. Individual
cases can be selected by name: `100000`, `200000_min150001`,
`100000_no_new_filters`, or `200000_threads`.

All required cases passed: 100000; 200000 with minimum 150001; 100000 without new
filters; and 200000 with one and four Rust threads. Both audit CSVs are
byte-identical to C++, including row order. Solution output and every non-timing
JSON value match, with only the specified generator/reconstruction/threads
differences. JSON key order is also checked. One- and four-thread CSVs match.

Measured `seconds_search` at bound 200000, with both audits enabled:

| Implementation | Threads | Seconds |
| --- | ---: | ---: |
| C++ v5 | 1 | 6.33156 |
| rsq6 | 1 | 6.276821416 |
| rsq6 | 4 | 1.807789833 |

Full output is retained in `accept/ACCEPT_LOG.txt`; reference and Rust
CSVs/stdout/stderr and validation logs are retained in `accept/`.

All task-created files are inside `rust_v6/`. No frozen C++ source was modified.
