//! v7 fast side search. Same pairs, same order, same counters as run_v5
//! (rust_v6, itself byte-identical to C++ v5); only the arithmetic used to
//! produce them changes:
//!  * rows with gcd(n,h)=1 (m=n): every CRT residue yields at most one u, decided
//!    without division; the shifts h*e_i mod n use a per-n Barrett constant and
//!    the test h % p^e is skipped (it is never 0 when gcd(n,h)=1);
//!  * rows with gcd>1: the lower end is split once per row (no division per residue);
//!  * gcd by the binary method;
//!  * the 11 local projections are tested through per-row bit masks.
use crate::arithmetic::*;
use crate::sieve::LocalSieve;
use crate::wide::*;
use crate::{SideResult, Tables, WorkerBuffers};
use std::fmt::Write as _;

/// Same order as LocalSieve::without19 (the order does not affect counters).
pub const LOCAL_PRIMES: [u32; 11] = [5, 11, 7, 29, 37, 13, 17, 23, 31, 41, 43];

pub struct LocalMasks {
    /// masks[i][(n%P)*P + h%P] has bit (u%P) set iff the projection mod P
    /// accepts (n, u, v=u-h).
    masks: Vec<Vec<u64>>,
}
impl LocalMasks {
    pub fn new(local: &LocalSieve) -> Self {
        let accept = |p: u32, n: u32, u: u32, v: u32| -> bool {
            match p {
                5 => local.p5.accepts(n, u, v),
                7 => local.p7.accepts(n, u, v),
                11 => local.p11.accepts(n, u, v),
                13 => local.p13.accepts(n, u, v),
                17 => local.p17.accepts(n, u, v),
                23 => local.p23.accepts(n, u, v),
                29 => local.p29.accepts(n, u, v),
                31 => local.p31.accepts(n, u, v),
                37 => local.p37.accepts(n, u, v),
                41 => local.p41.accepts(n, u, v),
                43 => local.p43.accepts(n, u, v),
                _ => unreachable!(),
            }
        };
        let masks = LOCAL_PRIMES
            .iter()
            .map(|&p| {
                let mut m = vec![0u64; (p * p) as usize];
                for n in 0..p {
                    for h in 0..p {
                        let mut bits = 0u64;
                        for u in 0..p {
                            if accept(p, n, u, (u + p - h) % p) {
                                bits |= 1 << u;
                            }
                        }
                        m[(n * p + h) as usize] = bits;
                    }
                }
                m
            })
            .collect();
        Self { masks }
    }
    /// Exhaustive check against LocalSieve::without19 on all residues.
    pub fn verify(&self, local: &LocalSieve) {
        let modulus: u32 = 5 * 7 * 11 * 13;
        let mut checked = 0u64;
        for n in [0u32, 1, 12, 24, 60, 1001, 5005, 9240, 46189, 3_999_996, 4_999_992] {
            for h in (0..modulus).step_by(7) {
                let row = self.row(n, h);
                for u in (h..h + 3 * modulus).step_by(13) {
                    assert_eq!(
                        accepts_row(&row, u),
                        local.without19(n, u, u - h),
                        "local mask mismatch n={n} h={h} u={u}"
                    );
                    checked += 1;
                }
            }
        }
        eprintln!("V7 local mask validation passed: {checked} comparisons.");
    }
    #[inline(always)]
    pub fn row(&self, n: u32, h: u32) -> [u64; 11] {
        macro_rules! at {
            ($i:expr, $p:expr) => {
                self.masks[$i][((n % $p) * $p + h % $p) as usize]
            };
        }
        [
            at!(0, 5),
            at!(1, 11),
            at!(2, 7),
            at!(3, 29),
            at!(4, 37),
            at!(5, 13),
            at!(6, 17),
            at!(7, 23),
            at!(8, 31),
            at!(9, 41),
            at!(10, 43),
        ]
    }
}
#[inline(always)]
pub fn accepts_row(row: &[u64; 11], u: u32) -> bool {
    (row[0] >> (u % 5)) & 1 != 0
        && (row[1] >> (u % 11)) & 1 != 0
        && (row[2] >> (u % 7)) & 1 != 0
        && (row[3] >> (u % 29)) & 1 != 0
        && (row[4] >> (u % 37)) & 1 != 0
        && (row[5] >> (u % 13)) & 1 != 0
        && (row[6] >> (u % 17)) & 1 != 0
        && (row[7] >> (u % 23)) & 1 != 0
        && (row[8] >> (u % 31)) & 1 != 0
        && (row[9] >> (u % 41)) & 1 != 0
        && (row[10] >> (u % 43)) & 1 != 0
}

#[inline]
pub fn binary_gcd(mut a: u32, mut b: u32) -> u32 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        b >>= b.trailing_zeros();
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        b -= a;
        if b == 0 {
            return a << shift;
        }
    }
}

/// Exact a mod d for a < 2^64 with one 64x64->128 multiply (Barrett).
#[derive(Clone, Copy)]
pub struct Barrett {
    d: u64,
    inv: u64,
}
impl Barrett {
    pub fn new(d: u32) -> Self {
        Self { d: d as u64, inv: u64::MAX / d as u64 }
    }
    #[inline(always)]
    pub fn rem(self, a: u64) -> u64 {
        let q = ((a as u128 * self.inv as u128) >> 64) as u64;
        let mut r = a - q * self.d;
        while r >= self.d {
            r -= self.d;
        }
        r
    }
}

pub fn run_fast(
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
    let barrett_n = Barrett::new(n.max(1));
    let plan_count = plan.count;
    let mut idem_n = [0u64; 7];
    for i in 0..plan_count {
        idem_n[i] = plan.idempotent[i] as u64;
    }
    let permitted = &b.permitted;
    let active = &b.active;
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
        let gcd_nh = binary_gcd(n, h);
        let (lower, upper) = if trim { (h + 1, n - 1) } else { (h, n) };
        if lower > upper {
            return; // generate_field_gap_crt_range returns 0 residues here
        }
        let masks = t.masks.row(n, h);
        let mut visit = |u: u32| {
            let v = u - h;
            st.c.generated_pairs += 1;
            if permitted[u as usize] == 0 || permitted[v as usize] == 0 {
                return;
            }
            st.c.valuation_pairs += 1;
            st.c.gamma_early_passed += 1;
            if !accepts_row(&masks, u) {
                return;
            }
            st.c.local_pairs += 1;
            if let Some(additional) = &t.additional {
                if !additional.accepts(n, u, v, active) {
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
        let mut residues = [0u32; 128];
        let mut count = 1usize;
        if gcd_nh == 1 {
            // m = n. h shares no prime with n, so no prime power divides h.
            let m = n;
            for i in 0..plan_count {
                let shift = barrett_n.rem(h as u64 * idem_n[i]) as u32;
                let old = count;
                assert!(2 * old <= 128, "CRT capacity");
                for j in 0..old {
                    let mut r = residues[j] + shift;
                    if r >= m {
                        r -= m;
                    }
                    residues[count] = r;
                    count += 1;
                }
            }
            let row = t.wheel.row(n, h, m);
            let top_ok = upper == n && row[(n % 19) as usize] == 0;
            for &r in &residues[..count] {
                // The progression r, r+n, ... meets [lower,upper] in r (r>=lower)
                // or in n (r=0, upper=n); a nonzero wheel skip leaves the range.
                if r >= lower {
                    if row[(r % 19) as usize] == 0 {
                        visit(r);
                    }
                } else if r == 0 && top_ok {
                    visit(n);
                }
            }
        } else {
            let m = n / gcd_nh;
            let mplan = &t.arithmetic.plans[m as usize];
            let bm = Barrett::new(m);
            for i in 0..mplan.count {
                if h % mplan.power[i] == 0 {
                    continue;
                }
                let shift = bm.rem(h as u64 * mplan.idempotent[i] as u64) as u32;
                let old = count;
                assert!(2 * old <= 128, "CRT capacity");
                for j in 0..old {
                    let mut r = residues[j] + shift;
                    if r >= m {
                        r -= m;
                    }
                    residues[count] = r;
                    count += 1;
                }
            }
            let row = t.wheel.row(n, h, m);
            let (ql, rl) = (lower / m, lower % m);
            for &r in &residues[..count] {
                let first = if r >= rl { ql * m + r } else { (ql + 1) * m + r };
                let mut u = first;
                while u <= upper {
                    let skip = row[(u % 19) as usize];
                    if skip == 255 {
                        break;
                    }
                    u += skip as u32 * m;
                    if u > upper {
                        break;
                    }
                    visit(u);
                    u += m;
                }
            }
        }
        st.c.crt_residues += count as u64;
    });
    out
}
