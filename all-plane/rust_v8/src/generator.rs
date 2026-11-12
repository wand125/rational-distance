//! The v3 CRT/wheel generator and v4 validation generators, preserving order.
use crate::arithmetic::*;
use crate::sieve::{LocalSieve, PrimeProjection};

pub struct GapWheel {
    residues: Vec<Vec<u32>>,
}
impl GapWheel {
    pub const MODULUS: u32 = 5 * 11 * 19;
    pub fn new(local: &LocalSieve) -> Self {
        fn projection<const P: usize>(p: &PrimeProjection<P>) -> Vec<u8> {
            let mut possible = vec![0; P * P];
            for n in 0..P as u32 {
                for h in 0..P as u32 {
                    for u in 0..P as u32 {
                        possible[n as usize * P + h as usize] |=
                            u8::from(p.accepts(n, u, (u + P as u32 - h) % P as u32));
                    }
                }
            }
            possible
        }
        let (p5, p11, p19) = (
            projection(&local.p5),
            projection(&local.p11),
            projection(&local.p19),
        );
        let mut residues = vec![Vec::new(); Self::MODULUS as usize];
        for n in 0..Self::MODULUS {
            for h in 0..Self::MODULUS {
                if p5[(n % 5 * 5 + h % 5) as usize] != 0
                    && p11[(n % 11 * 11 + h % 11) as usize] != 0
                    && p19[(n % 19 * 19 + h % 19) as usize] != 0
                {
                    residues[n as usize].push(h);
                }
            }
        }
        Self { residues }
    }
    pub fn generate(&self, n: u32, mut emit: impl FnMut(u32)) {
        let row = &self.residues[(n % Self::MODULUS) as usize];
        for block in (0..=n).step_by(Self::MODULUS as usize) {
            for &r in row {
                let h = block + r;
                if h > n {
                    break;
                }
                emit(h);
            }
        }
    }
}
pub struct ProgressionWheel<const W: usize> {
    next: Vec<u8>,
}
impl<const W: usize> ProgressionWheel<W> {
    pub fn new(local: &LocalSieve) -> Self {
        let mut s = Self {
            next: vec![255; W * W * W * W],
        };
        for n in 0..W {
            for h in 0..W {
                let mut allowed = vec![false; W];
                for (u, yes) in allowed.iter_mut().enumerate() {
                    let v = (u + W - h) % W;
                    *yes = (W % 5 != 0 || local.p5.accepts(n as u32, u as u32, v as u32))
                        && (W % 11 != 0 || local.p11.accepts(n as u32, u as u32, v as u32))
                        && (W % 19 != 0 || local.p19.accepts(n as u32, u as u32, v as u32));
                }
                for step in 0..W {
                    let offset = ((n * W + h) * W + step) * W;
                    let cycles = gcd(step as u32, W as u32) as usize;
                    let len = W / cycles;
                    for first in 0..cycles {
                        let mut nearest = -1i32;
                        for j in (0..2 * len).rev() {
                            let u = (first + j * step) % W;
                            if allowed[u] {
                                nearest = j as i32;
                            }
                            if j < len {
                                s.next[offset + u] = if nearest < 0 {
                                    255
                                } else {
                                    (nearest - j as i32) as u8
                                };
                            }
                        }
                    }
                }
            }
        }
        s
    }
    pub fn row(&self, n: u32, h: u32, step: u32) -> &[u8] {
        let offset = (((n as usize % W) * W + h as usize % W) * W + step as usize % W) * W;
        &self.next[offset..offset + W]
    }
}
// Shared exact CRT construction; append the shifted old residues in C++ order.
fn crt_residues(
    n: u32,
    t: u32,
    tables: &ArithmeticTables,
    known_gcd: u32,
) -> (u32, [u32; 128], usize) {
    let m = n / if known_gcd != 0 { known_gcd } else { gcd(n, t) };
    let plan = &tables.plans[m as usize];
    let mut residues = [0u32; 128];
    let mut count = 1;
    for i in 0..plan.count {
        if t % plan.power[i] == 0 {
            continue;
        }
        let shift = (t as u64 * plan.idempotent[i] as u64 % m as u64) as u32;
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
    (m, residues, count)
}
pub fn generate_crt(
    n: u32,
    u: u32,
    tables: &ArithmeticTables,
    mut emit: impl FnMut(u32),
    known_gcd: u32,
) -> u64 {
    let (m, residues, count) = crt_residues(n, u, tables, known_gcd);
    for &r in &residues[..count] {
        for v in (r..=u).step_by(m as usize) {
            emit(v);
        }
    }
    count as u64
}
pub fn generate_gap_crt(
    n: u32,
    h: u32,
    tables: &ArithmeticTables,
    mut emit: impl FnMut(u32),
    known_gcd: u32,
) -> u64 {
    let (m, residues, count) = crt_residues(n, h, tables, known_gcd);
    for &r in &residues[..count] {
        let mut first = r;
        if first < h {
            first += ((h - first + m - 1) / m) * m;
        }
        for u in (first..=n).step_by(m as usize) {
            emit(u);
        }
    }
    count as u64
}
pub fn generate_field_gap_crt_range(
    n: u32,
    h: u32,
    tables: &ArithmeticTables,
    wheel: &ProgressionWheel<19>,
    mut emit: impl FnMut(u32),
    known_gcd: u32,
    lower: u32,
    upper: u32,
) -> u64 {
    if lower > upper {
        return 0;
    }
    let (m, residues, count) = crt_residues(n, h, tables, known_gcd);
    let row = wheel.row(n, h, m);
    for &r in &residues[..count] {
        let mut first = r;
        if first < lower {
            first += ((lower - first + m - 1) / m) * m;
        }
        let mut u = first as u64;
        while u <= upper as u64 {
            let skip = row[(u % 19) as usize];
            if skip == 255 {
                break;
            }
            u += skip as u64 * m as u64;
            if u > upper as u64 {
                break;
            }
            emit(u as u32);
            u += m as u64;
        }
    }
    count as u64
}
pub fn generate_field_gap_crt(
    n: u32,
    h: u32,
    tables: &ArithmeticTables,
    wheel: &ProgressionWheel<19>,
    emit: impl FnMut(u32),
    known_gcd: u32,
) -> u64 {
    generate_field_gap_crt_range(n, h, tables, wheel, emit, known_gcd, h, n)
}

// Gamma labels/rules/lifting are used only by the required v4 self-check.
pub struct GammaLabels {
    pub value: Vec<u32>,
}
impl GammaLabels {
    pub fn new(n: u32, plan: &GammaPowerPlan) -> Self {
        if n == 0 || plan.entries.is_empty() {
            return Self { value: Vec::new() };
        }
        let mut value = vec![1; n as usize + 1];
        for &(p, _, e) in &plan.entries {
            let mut pk = 1;
            for k in 1..e {
                pk *= p;
                let next = pk * p;
                let factor = p.pow(k.min(e - k));
                for t in (pk..=n).step_by(pk as usize) {
                    if t % next != 0 {
                        value[t as usize] *= factor;
                    }
                }
            }
        }
        Self { value }
    }
}
struct Rule {
    kind: u32,
    p: u32,
    modulus: u32,
    lower_next: u32,
    high: u32,
}
pub struct GammaPowerRow {
    h: u32,
    rules: Vec<Rule>,
}
impl GammaPowerRow {
    pub fn new(h: u32, plan: &GammaPowerPlan) -> Self {
        let mut rules = Vec::new();
        for &(p, power, e) in &plan.entries {
            let mut k = 0;
            let mut hk = 1;
            while k < e && h % (hk * p) == 0 {
                hk *= p;
                k += 1;
            }
            if k == 0 {
                continue;
            }
            rules.push(if k == e {
                Rule {
                    kind: 1,
                    p,
                    modulus: power,
                    lower_next: 0,
                    high: 0,
                }
            } else if 2 * k <= e {
                Rule {
                    kind: 0,
                    p,
                    modulus: p.pow(e - k + 1),
                    lower_next: 0,
                    high: 0,
                }
            } else {
                Rule {
                    kind: 2,
                    p,
                    modulus: hk * p,
                    lower_next: p.pow(e - k + 1),
                    high: hk,
                }
            });
        }
        Self { h, rules }
    }
    fn accepts_rule(&self, i: usize, u: u32) -> bool {
        let r = &self.rules[i];
        match r.kind {
            0 => u % r.modulus != 0 && (u - self.h) % r.modulus != 0,
            1 => u % r.p != 0 || u % r.modulus == 0,
            _ => {
                u % r.lower_next != 0
                    || (u % r.high == 0 && u % r.modulus != 0 && (u - self.h) % r.modulus != 0)
            }
        }
    }
    pub fn accepts(&self, u: u32) -> bool {
        self.accepts_mask(u, 0)
    }
    fn accepts_mask(&self, u: u32, mask: u32) -> bool {
        (0..self.rules.len()).all(|i| mask & (1 << i) != 0 || self.accepts_rule(i, u))
    }
}
#[derive(Default)]
pub struct GenerationStats {
    pub lifted_progressions: u64,
    pub lift_residues: u64,
    pub examined_rejections: u64,
}
// This recursive helper mirrors v4's refine lambda; the late search never calls it.
struct GammaGenerator<'a, F> {
    n: u32,
    h: u32,
    upper: u32,
    wheel: &'a ProgressionWheel<19>,
    gamma: &'a GammaPowerRow,
    lift: bool,
    emit: F,
    stats: &'a mut GenerationStats,
}
impl<F: FnMut(u32)> GammaGenerator<'_, F> {
    fn refine(&mut self, first: u32, step: u32, index: usize, mask: u32) {
        if index == self.gamma.rules.len() {
            let row = self.wheel.row(self.n, self.h, step);
            let all = (1u32 << self.gamma.rules.len()) - 1;
            let mut u = first as u64;
            while u <= self.upper as u64 {
                let skip = row[(u % 19) as usize];
                if skip == 255 {
                    break;
                }
                u += skip as u64 * step as u64;
                if u > self.upper as u64 {
                    break;
                }
                if mask != all && !self.gamma.accepts_mask(u as u32, mask) {
                    self.stats.examined_rejections += 1;
                } else {
                    (self.emit)(u as u32);
                }
                u += step as u64;
            }
            return;
        }
        let modulus = self.gamma.rules[index].modulus;
        let period = modulus / gcd(modulus, step);
        let length = (self.upper - first) / step + 1;
        if !self.lift || period > 61 || length < 2 * period {
            self.refine(first, step, index + 1, mask);
            return;
        }
        self.stats.lifted_progressions += 1;
        let new_step = step * period;
        for j in 0..period {
            if first as u64 + j as u64 * step as u64 > self.upper as u64 {
                break;
            }
            let start = first + j * step;
            if !self.gamma.accepts_rule(index, start) {
                continue;
            }
            self.stats.lift_residues += 1;
            self.refine(start, new_step, index + 1, mask | (1 << index));
        }
    }
}
pub fn generate_gamma_field_crt(
    n: u32,
    h: u32,
    tables: &ArithmeticTables,
    wheel: &ProgressionWheel<19>,
    gamma: &GammaPowerRow,
    lift: bool,
    emit: impl FnMut(u32),
    stats: &mut GenerationStats,
    trim: bool,
) -> u64 {
    let (lower, upper) = (h + u32::from(trim), n - u32::from(trim));
    if lower > upper {
        return 0;
    }
    let (m, residues, count) = crt_residues(n, h, tables, gcd(n, h));
    let mut generator = GammaGenerator {
        n,
        h,
        upper,
        wheel,
        gamma,
        lift,
        emit,
        stats,
    };
    for &r in &residues[..count] {
        let mut first = r;
        if first < lower {
            first += ((lower - first + m - 1) / m) * m;
        }
        if first <= upper {
            generator.refine(first, m, 0, 0);
        }
    }
    count as u64
}
