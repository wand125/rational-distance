//! rpart v2 (v1 + mod-12 bucketed pair loop, identical counters): whole-plane box search for points at rational distance from the four
//! vertices of the unit square, by Pythagorean-partner enumeration.
//!
//! A point (x,y) = (A/Q, C/Q) with B = Q-A, D = Q-C (signed integers, Q >= 1).
//! The four distances are rational iff A²+C², B²+C², A²+D², B²+D² are squares.
//! Write a=|A|, b=|B|, c=|C|, d=|D|. Then {a,b} ⊂ T(c) ∩ T(d) where
//! T(c) = { v >= 1 : c²+v² is a square }, and
//!   Q = a+b (0<x<1) or Q = |a-b| (x outside [0,1]),
//!   d = |Q-c| (0<=y<=1 or y>1) or d = Q+c (y<0).
//! The search covers every (Q, A, C) with Q <= N and max(a,b,c,d) <= M,
//! i.e. for M = K*N every point of the box [1-K, K]² with a common denominator <= N.
//! Cases a = 0 or b = 0 (point on x=0 or x=1) are covered through the symmetry
//! (x,y) <-> (y,x) together with d = 0 (point on y=0 or y=1). a = b only for
//! x = 1/2, covered through the same symmetry unless also y = 1/2 (the centre,
//! whose distances are √2/2).
//!
//! Filters (mode "filtered"): a solution with minimal denominator Q has
//! gcd(A,B,C,D) = 1, hence (mod 8 / mod 3 lemmas, valid for signed integers):
//!   c ≡ 2 (mod 4) impossible; c odd => a,b ≡ 0 (mod 4); 4|c => a,b odd; 12 | Q.
//! Mode "all" uses none of these (for validation against brute force).
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

const SEG: u64 = 1 << 18;
const MAXF: usize = 15;

#[derive(Default, Clone, Debug)]
struct Stats {
    c_values: u64,
    partners: u64,
    pairs: u64,
    q_candidates: u64,
    d_tests: u64,
    three: u64,
    solutions: u64,
    divisor_nodes: u64,
    // Mulatram-compatible interior counters (only in --interior-stats):
    // partners with v <= N, pairs a<b in T(c) with c < a+b <= N.
    int_partners: u64,
    int_candidates: u64,
}
impl Stats {
    fn add(&mut self, o: &Stats) {
        self.c_values += o.c_values;
        self.partners += o.partners;
        self.pairs += o.pairs;
        self.q_candidates += o.q_candidates;
        self.d_tests += o.d_tests;
        self.three += o.three;
        self.solutions += o.solutions;
        self.divisor_nodes += o.divisor_nodes;
        self.int_partners += o.int_partners;
        self.int_candidates += o.int_candidates;
    }
}

struct Opt {
    n: u64,
    m: u64,
    c_lo: u64,
    c_hi: u64,
    threads: usize,
    filtered: bool,
    dump_three: Option<String>,
    interior_stats: bool,
}

fn parse() -> Opt {
    let a: Vec<String> = std::env::args().collect();
    let mut o = Opt {
        n: 0,
        m: 0,
        c_lo: 1,
        c_hi: 0,
        threads: 1,
        filtered: true,
        dump_three: None,
        interior_stats: false,
    };
    let mut i = 1;
    while i < a.len() {
        let v = |i: usize| a.get(i + 1).expect("missing value").clone();
        match a[i].as_str() {
            "--n" => { o.n = v(i).parse().unwrap(); i += 1; }
            "--m" => { o.m = v(i).parse().unwrap(); i += 1; }
            "--c-lo" => { o.c_lo = v(i).parse().unwrap(); i += 1; }
            "--c-hi" => { o.c_hi = v(i).parse().unwrap(); i += 1; }
            "--threads" => { o.threads = v(i).parse().unwrap(); i += 1; }
            "--all" => o.filtered = false,
            "--dump-three" => { o.dump_three = Some(v(i)); i += 1; }
            "--interior-stats" => o.interior_stats = true,
            x => panic!("unknown argument {x}"),
        }
        i += 1;
    }
    assert!(o.n >= 1 && o.m >= 1, "--n and --m are required");
    assert!(o.m <= 1u64 << 40, "M too large for the u128 squares used");
    if o.c_hi == 0 { o.c_hi = o.m; }
    assert!(o.c_lo >= 1 && o.c_lo <= o.c_hi && o.c_hi <= o.m);
    o
}

fn small_primes(limit: u64) -> Vec<u64> {
    let l = limit as usize;
    let mut s = vec![true; l + 1];
    let mut p = Vec::new();
    for i in 2..=l {
        if s[i] {
            p.push(i as u64);
            let mut j = i * i;
            while j <= l { s[j] = false; j += i; }
        }
    }
    p
}

fn isqrt_u64(x: u64) -> u64 {
    let mut r = (x as f64).sqrt() as u64;
    while r.saturating_mul(r) > x { r -= 1; }
    while (r + 1).saturating_mul(r + 1) <= x { r += 1; }
    r
}

/// Exact square test of x < 2^82 (squares are computed exactly in u128).
#[inline]
fn is_square(x: u128) -> bool {
    // quadratic residues mod 64 (exact: computed from all residues)
    const Q64: u64 = {
        let mut m = 0u64;
        let mut i = 0;
        while i < 64 { m |= 1 << ((i * i) % 64); i += 1; }
        m
    };
    if (Q64 >> (x as u64 & 63)) & 1 == 0 { return false; }
    let r63 = (x % 63) as u32;
    const Q63: u64 = {
        let mut m = 0u64;
        let mut i = 0;
        while i < 63 { m |= 1 << ((i * i) % 63); i += 1; }
        m
    };
    if (Q63 >> r63) & 1 == 0 { return false; }
    let mut r = (x as f64).sqrt() as u128;
    while r * r > x { r -= 1; }
    while (r + 1) * (r + 1) <= x { r += 1; }
    r * r == x
}

const fn qr_mask(m: u64) -> u128 {
    let mut mask = 0u128;
    let mut i = 0;
    while i < m { mask |= 1u128 << ((i * i) % m); i += 1; }
    mask
}
const QR64: u128 = qr_mask(64);
const QR63: u128 = qr_mask(63);
const QR65: u128 = qr_mask(65);
const QR11: u128 = qr_mask(11);

/// Squares of v modulo 64, 63, 65, 11 (packed), for the residue prefilter.
#[derive(Clone, Copy)]
struct Sq { m64: u8, m63: u8, m65: u8, m11: u8 }
#[inline(always)]
fn sq_res(v: u64) -> Sq {
    let (r63, r65, r11) = (v % 63, v % 65, v % 11);
    Sq { m64: (v.wrapping_mul(v) & 63) as u8, m63: (r63 * r63 % 63) as u8, m65: (r65 * r65 % 65) as u8, m11: (r11 * r11 % 11) as u8 }
}
/// Exact test that a² + d² is a perfect square (a, d < 2^40). The residue tests
/// mod 64, 63, 65, 11 use tables built from all residues, so they only reject
/// non-squares; survivors get the exact u128 test.
#[inline(always)]
fn sum_sq_is_square(a: u64, sa: Sq, d: u64, sd: Sq) -> bool {
    if (QR64 >> ((sa.m64 + sd.m64) & 63)) & 1 == 0 { return false; }
    let mut x = (sa.m63 + sd.m63) as u32; if x >= 63 { x -= 63; }
    if (QR63 >> x) & 1 == 0 { return false; }
    let mut y = (sa.m65 + sd.m65) as u32; if y >= 65 { y -= 65; }
    if (QR65 >> y) & 1 == 0 { return false; }
    let mut z = (sa.m11 + sd.m11) as u32; if z >= 11 { z -= 11; }
    if (QR11 >> z) & 1 == 0 { return false; }
    is_square(a as u128 * a as u128 + d as u128 * d as u128)
}

/// Factorisation of every value in a segment [lo, lo+len).
struct SegFactors {
    cnt: Vec<u8>,
    p: Vec<u64>,
    e: Vec<u8>,
}
fn factor_segment(lo: u64, len: u64, primes: &[u64], f: &mut SegFactors) {
    let len = len as usize;
    f.cnt.clear(); f.cnt.resize(len, 0);
    f.p.resize(len * MAXF, 0);
    f.e.resize(len * MAXF, 0);
    let mut rem: Vec<u64> = (0..len as u64).map(|i| lo + i).collect();
    let hi = lo + len as u64;
    for &p in primes {
        if p * p >= hi { break; }
        let mut j = ((lo + p - 1) / p) * p;
        while j < hi {
            let i = (j - lo) as usize;
            let mut e = 0u8;
            while rem[i] % p == 0 { rem[i] /= p; e += 1; }
            let k = f.cnt[i] as usize;
            f.p[i * MAXF + k] = p;
            f.e[i * MAXF + k] = e;
            f.cnt[i] += 1;
            j += p;
        }
    }
    for i in 0..len {
        if rem[i] > 1 {
            let k = f.cnt[i] as usize;
            f.p[i * MAXF + k] = rem[i];
            f.e[i * MAXF + k] = 1;
            f.cnt[i] += 1;
        }
    }
}

/// T(c) ∩ [1, m]: v = (w-u)/2 for u*w = c², u < w, u ≡ w (mod 2).
/// Primes are visited largest first; w_partial (product of the complementary
/// powers of the primes assigned so far) only grows, so a node with
/// w_partial > 2m + c cannot reach a leaf with v <= m (v = (w-u)/2 <= m needs
/// w <= 2m + u < 2m + c).
fn partners(c: u64, ps: &[(u64, u8)], m: u64, out: &mut Vec<u64>, nodes: &mut u64) {
    out.clear();
    let cap = 2 * m as u128 + c as u128;
    let even = c % 2 == 0;
    fn rec(
        k: usize, u: u128, w: u128, ps: &[(u64, u8)], c: u128, cap: u128, even: bool,
        m: u64, out: &mut Vec<u64>, nodes: &mut u64,
    ) {
        *nodes += 1;
        if k == ps.len() {
            if u < c && (w - u) % 2 == 0 {
                let v = (w - u) / 2;
                if v >= 1 && v <= m as u128 { out.push(v as u64); }
            }
            return;
        }
        let (p, e) = ps[k];
        let p = p as u128;
        let e = e as u32;
        let (flo, fhi) = if even && p == 2 { (1, 2 * e - 1) } else { (0, 2 * e) };
        // f = exponent in u; complementary exponent 2e-f in w.
        // Iterate f from high to low so w grows; stop once w exceeds cap.
        let mut uu = u * p.pow(flo);
        for f in flo..=fhi {
            let wexp = 2 * e - f;
            let ww = w * p.pow(wexp);
            if ww <= cap && uu < c {
                rec(k + 1, uu, ww, ps, c, cap, even, m, out, nodes);
            }
            if uu >= c { break; }
            uu *= p;
        }
    }
    rec(0, 1, 1, ps, c as u128, cap, even, m, out, nodes);
}

fn main() {
    let o = parse();
    let started = std::time::Instant::now();
    let primes = small_primes(isqrt_u64(o.c_hi) + 1);
    let next = AtomicU64::new(o.c_lo);
    let total = Mutex::new(Stats::default());
    let three_rows: Mutex<Vec<(u64, u64, u64, u64, u64)>> = Mutex::new(Vec::new());
    let solutions: Mutex<Vec<String>> = Mutex::new(Vec::new());
    let dump = o.dump_three.is_some();
    std::thread::scope(|s| {
        for _ in 0..o.threads {
            s.spawn(|| {
                let mut st = Stats::default();
                let mut f = SegFactors { cnt: vec![], p: vec![], e: vec![] };
                let mut t: Vec<u64> = Vec::new();
                let mut ps: Vec<(u64, u8)> = Vec::new();
                let mut local_three = Vec::new();
                let mut buckets: Vec<Vec<u32>> = (0..12).map(|_| Vec::new()).collect();
                let mut ts: Vec<Sq> = Vec::new();
                let (n, m) = (o.n, o.m);
                loop {
                    let lo = next.fetch_add(SEG, Ordering::Relaxed);
                    if lo > o.c_hi { break; }
                    let len = SEG.min(o.c_hi + 1 - lo);
                    factor_segment(lo, len, &primes, &mut f);
                    for i in 0..len as usize {
                        let c = lo + i as u64;
                        if o.filtered && c % 4 == 2 { continue; }
                        st.c_values += 1;
                        ps.clear();
                        for k in 0..f.cnt[i] as usize {
                            ps.push((f.p[i * MAXF + k], f.e[i * MAXF + k]));
                        }
                        ps.sort_unstable_by(|x, y| y.0.cmp(&x.0));
                        partners(c, &ps, m, &mut t, &mut st.divisor_nodes);
                        if o.interior_stats && c <= n {
                            t.sort_unstable();
                            let tn: Vec<u64> = t.iter().copied().filter(|&v| v <= n).collect();
                            st.int_partners += tn.len() as u64;
                            for x in 0..tn.len() {
                                for y in x + 1..tn.len() {
                                    let s2 = tn[x] + tn[y];
                                    if s2 > n { break; }
                                    if s2 > c { st.int_candidates += 1; }
                                }
                            }
                        }
                        if o.filtered {
                            if c % 2 == 1 { t.retain(|&v| v % 4 == 0); } else { t.retain(|&v| v % 2 == 1); }
                        }
                        st.partners += t.len() as u64;
                        t.sort_unstable();
                        ts.clear();
                        ts.extend(t.iter().map(|&v| sq_res(v)));
                        if o.filtered {
                            // v2: pairs counted by two pointers; only b ≡ ±a (mod 12) can give 12 | Q
                            // (a+b ≡ 0 or b-a ≡ 0), so only those buckets are visited. Same counters as v1.
                            for v in buckets.iter_mut() { v.clear(); }
                            let mut hi = 0usize;
                            for x in 0..t.len() {
                                if hi < x + 1 { hi = x + 1; }
                                while hi < t.len() && t[hi] - t[x] <= n { hi += 1; }
                                st.pairs += (hi - x - 1) as u64;
                            }
                            for (i, &v) in t.iter().enumerate() { buckets[(v % 12) as usize].push(i as u32); }
                            for x in 0..t.len() {
                                let a = t[x];
                                let r1 = (a % 12) as usize;
                                let r2 = ((12 - r1) % 12) as usize;
                                for (k, &r) in [r1, r2].iter().enumerate() {
                                    if k == 1 && r2 == r1 { break; }
                                    let lst = &buckets[r];
                                    let mut j = lst.partition_point(|&i| (i as usize) <= x);
                                    while j < lst.len() {
                                        let yi = lst[j] as usize;
                                        let b = t[yi];
                                        if b - a > n { break; }
                                        j += 1;
                                        let qs = [a + b, b - a];
                                        for (qi, &q) in qs.iter().enumerate() {
                                            if q > n { continue; }
                                            if o.filtered && q % 12 != 0 { continue; }
                                            let _ = qi;
                                            st.q_candidates += 1;
                                            let ds = [q + c, if q >= c { q - c } else { c - q }];
                                            for &d in &ds {
                                                if d > m { continue; }
                                                st.d_tests += 1;
                                                let sd = sq_res(d);
                                                let sa = sum_sq_is_square(a, ts[x], d, sd);
                                                let sb = sum_sq_is_square(b, ts[yi], d, sd);
                                                if sa || sb {
                                                    st.three += 1;
                                                    if dump { local_three.push((q, a, b, c, d)); }
                                                }
                                                if sa && sb {
                                                    st.solutions += 1;
                                                    solutions.lock().unwrap().push(format!(
                                                        "{{\"solution\":{{\"Q\":{q},\"a\":{a},\"b\":{b},\"c\":{c},\"d\":{d}}}}}"
                                                    ));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            for x in 0..t.len() {
                                let a = t[x];
                                for y in x + 1..t.len() {
                                    let b = t[y];
                                    if b - a > n { break; }
                                    st.pairs += 1;
                                    if o.filtered && (a % 3 == 0) != (b % 3 == 0) { continue; }
                                    let qs = [a + b, b - a];
                                    for (qi, &q) in qs.iter().enumerate() {
                                        if q > n { continue; }
                                        if o.filtered && q % 12 != 0 { continue; }
                                        let _ = qi;
                                        st.q_candidates += 1;
                                        let ds = [q + c, if q >= c { q - c } else { c - q }];
                                        for &d in &ds {
                                            if d > m { continue; }
                                            st.d_tests += 1;
                                            let sd = sq_res(d);
                                            let sa = sum_sq_is_square(a, ts[x], d, sd);
                                            let sb = sum_sq_is_square(b, ts[y], d, sd);
                                            if sa || sb {
                                                st.three += 1;
                                                if dump { local_three.push((q, a, b, c, d)); }
                                            }
                                            if sa && sb {
                                                st.solutions += 1;
                                                solutions.lock().unwrap().push(format!(
                                                    "{{\"solution\":{{\"Q\":{q},\"a\":{a},\"b\":{b},\"c\":{c},\"d\":{d}}}}}"
                                                ));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                total.lock().unwrap().add(&st);
                three_rows.lock().unwrap().extend(local_three);
            });
        }
    });
    let st = total.into_inner().unwrap();
    for s in solutions.into_inner().unwrap() { println!("{s}"); }
    if let Some(path) = &o.dump_three {
        let mut rows = three_rows.into_inner().unwrap();
        rows.sort_unstable();
        rows.dedup();
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        writeln!(w, "Q,a,b,c,d").unwrap();
        for r in rows { writeln!(w, "{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4).unwrap(); }
    }
    println!(
        "{{\"engine\":\"rpart_v2\",\"N\":{},\"M\":{},\"c_lo\":{},\"c_hi\":{},\"filtered\":{},\"threads\":{},\
\"c_values\":{},\"partners\":{},\"pairs\":{},\"q_candidates\":{},\"d_tests\":{},\"three\":{},\"solutions\":{},\
\"divisor_nodes\":{},\"int_partners\":{},\"int_candidates\":{},\"seconds\":{:.3}}}",
        o.n, o.m, o.c_lo, o.c_hi, o.filtered, o.threads, st.c_values, st.partners, st.pairs,
        st.q_candidates, st.d_tests, st.three, st.solutions, st.divisor_nodes, st.int_partners,
        st.int_candidates, started.elapsed().as_secs_f64()
    );
}
