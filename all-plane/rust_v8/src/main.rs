mod arithmetic;
mod fast;
mod generator;
mod sieve;
mod validation;
mod wide;

use arithmetic::*;
use generator::*;
use sieve::*;
use std::collections::BTreeMap;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::{
    atomic::{AtomicU32, Ordering},
    mpsc, Condvar, Mutex,
};
use std::time::Instant;
use wide::*;

macro_rules! counters {
    ($name:ident { $($field:ident),* $(,)? })=>{
        #[derive(Default,Clone,Debug,PartialEq,Eq)]
        pub struct $name {$(pub $field:u64),*}
        impl $name {fn add(&mut self,other:&Self) {$(self.$field+=other.$field;)*}}
    }
}
counters!(Stats {
    sides,
    parameter_pairs,
    linear_pairs,
    square_discriminants,
    integral_a,
    integer_coordinates,
    primitive_three_candidates,
    solutions
});
counters!(SearchCounters {
    dense_even_pairs,
    gap_rows,
    local_gap_rows,
    valuation_gap_rows,
    difference_gap_rows,
    critical_gap_rows,
    endpoint_empty_rows,
    crt_residues,
    generated_pairs,
    valuation_pairs,
    gamma_early_rejected,
    gamma_early_passed,
    gamma_lifted_progressions,
    gamma_lift_residues,
    examined_generation_rejections,
    local_pairs,
    additional_pairs,
    additional_rejected,
    scale_pairs,
    gamma_pairs,
    gamma_rejected,
    critical_projection_pairs,
    critical_projection_rejected,
});
#[derive(Default)]
struct V4Stats {
    c: SearchCounters,
    exact: Stats,
}
impl V4Stats {
    fn add(&mut self, other: &Self) {
        self.c.add(&other.c);
        self.exact.add(&other.exact);
    }
}
struct Tables {
    arithmetic: ArithmeticTables,
    local: LocalSieve,
    gap: GapWheel,
    wheel: ProgressionWheel<19>,
    additional: Option<AdditionalSieve>,
    masks: fast::LocalMasks,
}
#[derive(Default)]
struct WorkerBuffers {
    permitted: Vec<u8>,
    active: Vec<(usize, u32)>,
}
#[derive(Default)]
struct SideResult {
    stats: V4Stats,
    parameters: String,
    three: String,
    solutions: String,
}

// Exactly run_v5<late,true> for one n. Only read-only tables are shared.
fn run_v5(
    n: u32,
    t: &Tables,
    b: &mut WorkerBuffers,
    audit_parameters: bool,
    audit_three: bool,
) -> SideResult {
    let mut out = SideResult::default();
    let n2 = 2 * n;
    let st = &mut out.stats;
    st.exact.sides += 1;
    st.c.dense_even_pairs += (n as u64 + 1) * (n as u64 + 2) / 2;
    st.c.gap_rows += n as u64 + 1;
    let plan = &t.arithmetic.plans[n as usize];
    permitted_valuations(n, plan, &mut b.permitted);
    let critical_rows = CriticalContentRows::new(plan);
    critical_rows.mark(&mut b.permitted, n);
    if let Some(additional) = &t.additional {
        additional.mark(n, &mut b.permitted, &mut b.active);
    }
    let gamma = GammaPowerPlan::new(plan);
    let critical_projection = CriticalProjectionPlan::new(plan);
    let trim = t.additional.is_some();
    t.gap.generate(n, |h| {
        st.c.local_gap_rows += 1;
        st.c.valuation_gap_rows += 1;
        if trim && h >= n - 1 {
            st.c.endpoint_empty_rows += 1;
            return;
        }
        if !permitted_difference(n, h, 0, &t.arithmetic) {
            return;
        }
        st.c.difference_gap_rows += 1;
        if !critical_rows.accepts(h) {
            return;
        }
        st.c.critical_gap_rows += 1;
        let gcd_nh = gcd(n, h);
        let visit = |u: u32| {
            let v = u - h;
            st.c.generated_pairs += 1;
            if b.permitted[u as usize] == 0 || b.permitted[v as usize] == 0 {
                return;
            }
            st.c.valuation_pairs += 1;
            st.c.gamma_early_passed += 1;
            if !t.local.without19(n, u, v) {
                return;
            }
            st.c.local_pairs += 1;
            if let Some(additional) = &t.additional {
                if !additional.accepts(n, u, v, &b.active) {
                    st.c.additional_rejected += 1;
                    return;
                }
            }
            st.c.additional_pairs += 1;
            if !permitted_scale(n, u, v, gcd_nh) {
                return;
            }
            st.c.scale_pairs += 1;
            if gamma.relevant > 1 {
                let gh = gamma_value(gamma.relevant, h);
                if gh != gamma_value(gamma.relevant, u) || gh != gamma_value(gamma.relevant, v) {
                    st.c.gamma_rejected += 1;
                    return;
                }
            }
            st.c.gamma_pairs += 1;
            if !critical_projection.accepts(u, h) {
                st.c.critical_projection_rejected += 1;
                return;
            }
            st.c.critical_projection_pairs += 1;
            if audit_parameters {
                writeln!(out.parameters, "{n2},{},{}", 2 * u, 2 * v).unwrap();
            }
            let mut solutions = 0;
            evaluate_v5(n2, 2 * u, 2 * v, true, &mut st.exact, |x, y, a| {
                if audit_three {
                    writeln!(out.three, "{n2},{},{},{x},{y},{a}", 2 * u, 2 * v).unwrap();
                }
                let d2 = add256(square_i128(n2 as i128 - x), square_i128(n2 as i128 - y));
                if let Some(d) = square256(d2) {
                    solutions += 1;
                    writeln!(
                        out.solutions,
                        "{{\"solution\":{{\"N\":{n2},\"X\":{x},\"Y\":{y},\"A\":{a},\"D\":{d}}}}}"
                    )
                    .unwrap();
                }
            });
            st.exact.solutions += solutions;
        };
        st.c.crt_residues += if trim {
            generate_field_gap_crt_range(n, h, &t.arithmetic, &t.wheel, visit, gcd_nh, h + 1, n - 1)
        } else {
            generate_field_gap_crt(n, h, &t.arithmetic, &t.wheel, visit, gcd_nh)
        };
    });
    out
}

struct Options {
    bound: u32,
    minimum: u32,
    threads: usize,
    validate: bool,
    new_filters: bool,
    reference: bool,
    three: Option<String>,
    parameters: Option<String>,
}
fn parse() -> Result<Option<Options>, String> {
    let mut opt = Options {
        bound: 10000,
        minimum: 1,
        threads: 1,
        validate: false,
        new_filters: true,
        reference: false,
        three: None,
        parameters: None,
    };
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || {
            args.next()
                .ok_or_else(|| "missing option value".to_string())
        };
        let range = |s: String, name: &str| -> Result<u32, String> {
            let v = s.parse::<u64>().map_err(|_| format!("{name}1..10000000"))?;
            if !(1..=MAX_DENOMINATOR as u64).contains(&v) {
                return Err(format!("{name}1..10000000"));
            }
            Ok(v as u32)
        };
        match arg.as_str() {
            "--validate" => opt.validate = true,
            "--no-new-filters" => opt.new_filters = false,
            "--reference-v6" => opt.reference = true,
            "--min-denominator" => opt.minimum = range(value()?, "minimum")?,
            "--threads" => {
                opt.threads = value()?
                    .parse()
                    .map_err(|_| "threads must be positive".to_string())?;
                if opt.threads == 0 {
                    return Err("threads must be positive".to_string());
                }
            }
            "--audit-three" => opt.three = Some(value()?),
            "--audit-parameters" => opt.parameters = Some(value()?),
            "--help" => {
                println!("rsq8 [B<=10000000] [--min-denominator L] [--threads T] [--audit-three CSV] [--audit-parameters CSV] [--no-new-filters] [--reference-v6] [--validate]");
                return Ok(None);
            }
            _ => opt.bound = range(arg, "bound")?,
        }
    }
    if opt.minimum > opt.bound {
        return Err("minimum exceeds bound".to_string());
    }
    Ok(Some(opt))
}
fn audit_file(path: &Option<String>, header: &[u8]) -> io::Result<Option<BufWriter<File>>> {
    path.as_ref()
        .map(|p| {
            let mut f = BufWriter::new(File::create(p)?);
            f.write_all(header)?;
            Ok(f)
        })
        .transpose()
}

fn search(
    opt: &Options,
    t: &Tables,
    parameters: &mut Option<BufWriter<File>>,
    three: &mut Option<BufWriter<File>>,
    stdout: &mut impl Write,
) -> io::Result<V4Stats> {
    let first = (opt.minimum + 23) / 24;
    let end = opt.bound / 24 + 1;
    let next = AtomicU32::new(first);
    // Limit outstanding sides, including completed results, to O(T). A slow
    // early side cannot cause the reorder buffer to grow with the full bound.
    let progress = (Mutex::new(first), Condvar::new());
    let window = opt.threads.saturating_mul(4).min(u32::MAX as usize) as u32;
    let (tx, rx) = mpsc::sync_channel(opt.threads.saturating_mul(2));
    let mut total = V4Stats::default();
    let mut error = None;
    std::thread::scope(|scope| {
        for _ in 0..opt.threads {
            let tx = tx.clone();
            let next = &next;
            let progress = &progress;
            scope.spawn(move || {
                let mut buffers = WorkerBuffers::default();
                loop {
                    let k = next.fetch_add(2, Ordering::Relaxed);
                    if k >= end {
                        break;
                    }
                    let mut emitted = progress.0.lock().unwrap();
                    while k.saturating_sub(*emitted) >= window {
                        emitted = progress.1.wait(emitted).unwrap();
                    }
                    drop(emitted);
                    for side in k..(k + 2).min(end) {
                        let run = if opt.reference { run_v5 } else { fast::run_fast };
                        let result = run(
                            12 * side,
                            t,
                            &mut buffers,
                            opt.parameters.is_some(),
                            opt.three.is_some(),
                        );
                        if tx.send((side, result)).is_err() {
                            return;
                        }
                    }
                }
            });
        }
        drop(tx);
        let mut pending = BTreeMap::new();
        let mut expected = first;
        for (k, result) in rx {
            assert!(pending.insert(k, result).is_none());
            while let Some(result) = pending.remove(&expected) {
                total.add(&result.stats);
                // Continue draining after an I/O failure so scoped workers
                // cannot deadlock on the channel or scheduling window.
                if error.is_none() {
                    let written = (|| {
                        if let Some(f) = parameters.as_mut() {
                            f.write_all(result.parameters.as_bytes())?;
                        }
                        if let Some(f) = three.as_mut() {
                            f.write_all(result.three.as_bytes())?;
                        }
                        if !result.solutions.is_empty() {
                            stdout.write_all(result.solutions.as_bytes())?;
                            stdout.flush()?;
                        }
                        Ok::<(), io::Error>(())
                    })();
                    if let Err(e) = written {
                        error = Some(e);
                    }
                }
                expected += 1;
                *progress.0.lock().unwrap() = expected;
                progress.1.notify_all();
            }
        }
        assert!(pending.is_empty());
        assert_eq!(expected, end.max(first));
    });
    if let Some(e) = error {
        Err(e)
    } else {
        Ok(total)
    }
}

fn final_json(opt: &Options, stats: &V4Stats, setup: f64, search: f64) -> String {
    let mut fields = Vec::new();
    macro_rules! raw {
        ($k:literal,$v:expr) => {
            fields.push(format!("\"{}\":{}", $k, $v));
        };
    }
    macro_rules! string {
        ($k:literal,$v:expr) => {
            fields.push(format!("\"{}\":\"{}\"", $k, $v));
        };
    }
    string!("scope", "whole_plane");
    string!("generator", "gap_CRT_v8_rust");
    string!("reconstruction", "exact_u256_v8_rust");
    string!("engine", if opt.reference { "reference_v6" } else { "fast_v8" });
    raw!("max_minimal_coordinate_denominator", opt.bound);
    raw!("min_denominator", opt.minimum);
    raw!("reference_v3", false);
    string!("gamma_stage", "late");
    raw!("skip_redundant19", true);
    raw!("additional_filters", opt.new_filters);
    raw!("endpoint_pruning", opt.new_filters);
    string!(
        "fifth_power_condition",
        if opt.new_filters { "early" } else { "off" }
    );
    raw!("prime_power_projections", opt.new_filters);
    raw!("critical_extension", opt.new_filters);
    raw!("fast_fifth_mark", opt.new_filters);
    let c = &stats.c;
    let e = &stats.exact;
    raw!("sides", e.sides);
    raw!("dense_even_pairs", c.dense_even_pairs);
    raw!("gap_rows", c.gap_rows);
    raw!("local_gap_rows", c.local_gap_rows);
    raw!("difference_gap_rows", c.difference_gap_rows);
    raw!("critical_gap_rows", c.critical_gap_rows);
    raw!("endpoint_empty_rows", c.endpoint_empty_rows);
    raw!("crt_residue_choices", c.crt_residues);
    raw!("generated_pairs", c.generated_pairs);
    raw!("valuation_pairs", c.valuation_pairs);
    raw!("gamma_early_rejected", c.gamma_early_rejected);
    raw!("gamma_lifted_progressions", c.gamma_lifted_progressions);
    raw!("gamma_lift_residues", c.gamma_lift_residues);
    raw!(
        "examined_generation_rejections",
        c.examined_generation_rejections
    );
    raw!("local_pairs", c.local_pairs);
    raw!("additional_pairs", c.additional_pairs);
    raw!("additional_rejected", c.additional_rejected);
    raw!("scale_pairs", c.scale_pairs);
    raw!("gamma_pairs", c.gamma_pairs);
    raw!("gamma_rejected", c.gamma_rejected);
    raw!("critical_projection_pairs", c.critical_projection_pairs);
    raw!(
        "critical_projection_rejected",
        c.critical_projection_rejected
    );
    raw!("exact_parameter_pairs", e.parameter_pairs);
    raw!("square_discriminants", e.square_discriminants);
    raw!("linear_pairs", e.linear_pairs);
    raw!("primitive_three_candidates", e.primitive_three_candidates);
    raw!("solutions", e.solutions);
    raw!("threads", opt.threads);
    raw!("seconds_setup", setup);
    raw!("seconds_search", search);
    raw!("seconds_total_excluding_validation", setup + search);
    format!("{{{}}}", fields.join(","))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let Some(opt) = parse().map_err(io::Error::other)? else {
        return Ok(());
    };
    let started = Instant::now();
    let arithmetic = ArithmeticTables::new(1000.max(opt.bound / 2));
    let local = LocalSieve::new();
    let gap = GapWheel::new(&local);
    let wheel = ProgressionWheel::new(&local);
    let additional = opt.new_filters.then(AdditionalSieve::new);
    let masks = fast::LocalMasks::new(&local);
    let tables = Tables {
        arithmetic,
        local,
        gap,
        wheel,
        additional,
        masks,
    };
    let setup = started.elapsed().as_secs_f64();
    if opt.validate {
        validation::validate(&tables);
        tables.masks.verify(&tables.local);
        fast::verify_small_valuation(opt.bound);
        validate_engines(&tables, opt.bound);
    }
    let mut three = audit_file(&opt.three, b"N,p,q,X,Y,A\n")?;
    let mut parameters = audit_file(&opt.parameters, b"N,p,q\n")?;
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let search_started = Instant::now();
    let stats = search(&opt, &tables, &mut parameters, &mut three, &mut stdout)?;
    let seconds_search = search_started.elapsed().as_secs_f64();
    if let Some(f) = three.as_mut() {
        f.flush()?;
    }
    if let Some(f) = parameters.as_mut() {
        f.flush()?;
    }
    writeln!(
        stdout,
        "{}",
        final_json(&opt, &stats, setup, seconds_search)
    )?;
    Ok(())
}

// Every side n=12k for k<=K0 and a spread of larger sides up to bound/2:
// the fast engine must return exactly the reference v6 result (counters,
// parameter rows, three-distance rows, solutions).
fn validate_engines(t: &Tables, bound: u32) {
    let kmax = bound / 24;
    let mut ks: Vec<u32> = (1..=kmax.min(700)).collect();
    if kmax > 700 {
        let step = ((kmax - 700) / 300).max(1);
        ks.extend((701..=kmax).step_by(step as usize));
        ks.extend(kmax.saturating_sub(5).max(701)..=kmax);
    }
    ks.sort_unstable();
    ks.dedup();
    let (mut a, mut b) = (WorkerBuffers::default(), WorkerBuffers::default());
    let mut pairs = 0u64;
    for &k in &ks {
        let n = 12 * k;
        let r = run_v5(n, t, &mut a, true, true);
        let f = fast::run_fast(n, t, &mut b, true, true);
        assert!(
            r.stats.c == f.stats.c
                && r.stats.exact == f.stats.exact
                && r.parameters == f.parameters
                && r.three == f.three
                && r.solutions == f.solutions,
            "fast engine differs from reference v6 at n={n}"
        );
        pairs += r.stats.c.generated_pairs;
    }
    eprintln!(
        "V8 engine validation passed: {} sides up to n={}, {pairs} generated pairs identical to reference v6.",
        ks.len(),
        12 * ks.last().copied().unwrap_or(0)
    );
}
