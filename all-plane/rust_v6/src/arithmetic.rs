//! ArithmeticTables and whole-plane necessary predicates from v2/v3/v4.
pub fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
pub fn inverse_mod(a: u32, modulus: u32) -> u32 {
    let (mut r, mut s, mut t, mut w) = (modulus as i64, a as i64, 0i64, 1i64);
    while s != 0 {
        let q = r / s;
        (r, s) = (s, r - q * s);
        (t, w) = (w, t - q * w);
    }
    assert_eq!(r, 1, "inverse of nonunit");
    t.rem_euclid(modulus as i64) as u32
}
#[derive(Clone, Default)]
pub struct CRTPlan {
    pub prime: [u32; 7],
    pub power: [u32; 7],
    pub idempotent: [u32; 7],
    pub count: usize,
}
pub struct ArithmeticTables {
    pub spf: Vec<u32>,
    pub bad_part: Vec<u32>,
    pub plans: Vec<CRTPlan>,
}
impl ArithmeticTables {
    pub fn new(bound: u32) -> Self {
        let mut s = Self {
            spf: vec![0; bound as usize + 1],
            bad_part: vec![1; bound as usize + 1],
            plans: vec![CRTPlan::default(); bound as usize + 1],
        };
        for p in 2..=bound {
            if s.spf[p as usize] != 0 {
                continue;
            }
            s.spf[p as usize] = p;
            if p as u64 * p as u64 <= bound as u64 {
                for j in (p * p..=bound).step_by(p as usize) {
                    if s.spf[j as usize] == 0 {
                        s.spf[j as usize] = p;
                    }
                }
            }
        }
        s.bad_part[0] = 0;
        for m in 2..=bound {
            let p = s.spf[m as usize];
            s.bad_part[m as usize] =
                s.bad_part[(m / p) as usize] * if p % 8 == 3 || p % 8 == 5 { p } else { 1 };
            let mut rem = m;
            let plan = &mut s.plans[m as usize];
            while rem > 1 {
                let prime = s.spf[rem as usize];
                let mut power = 1;
                loop {
                    power *= prime;
                    rem /= prime;
                    if rem % prime != 0 {
                        break;
                    }
                }
                let i = plan.count;
                assert!(i < 7, "CRT factor capacity");
                plan.count += 1;
                plan.prime[i] = prime;
                plan.power[i] = power;
                let other = m / power;
                plan.idempotent[i] =
                    (other as u64 * inverse_mod(other % power, power) as u64 % m as u64) as u32;
            }
        }
        s
    }
}
pub fn permitted_valuations(n: u32, plan: &CRTPlan, permitted: &mut Vec<u8>) {
    permitted.clear();
    permitted.resize(n as usize + 1, 1);
    for i in 0..plan.count {
        let (p, power) = (plan.prime[i], plan.power[i]);
        if !(p == 2 || p % 4 == 3) || power == p {
            continue;
        }
        for base in (0..=n).step_by(power as usize) {
            let mut j = p;
            while j < power && base + j <= n {
                permitted[(base + j) as usize] = 0;
                j += p;
            }
        }
    }
    if n % 4 == 0 {
        let two = n & n.wrapping_neg();
        for u in (2 * two..=n).step_by((4 * two) as usize) {
            permitted[u as usize] = 0;
        }
    }
}
pub fn permitted_difference(n: u32, u: u32, v: u32, tables: &ArithmeticTables) -> bool {
    let delta = u - v;
    if delta == 0 {
        return false;
    }
    if delta & 1 == 0 {
        let two = n & n.wrapping_neg();
        if delta & (2 * two - 1) != two {
            return false;
        }
    }
    let bad = tables.bad_part[delta as usize];
    bad == 1 || n % bad == 0
}
pub fn permitted_scale(n: u32, u: u32, v: u32, gcd_nu: u32) -> bool {
    if gcd_nu == 1 {
        return true;
    }
    let g = gcd(gcd_nu, v);
    if g == 1 {
        return true;
    }
    let (n0, u0, v0) = ((n / g) as i64, (u / g) as i64, (v / g) as i64);
    let q0 = u0 * u0 + v0 * v0 - n0 * n0;
    let required = g / gcd(g, n0 as u32);
    if required > 1 && q0 % required as i64 != 0 {
        return false;
    }
    if g & 1 == 0 {
        let two = g & g.wrapping_neg();
        if q0 % (2 * two) as i64 != 0 {
            return false;
        }
    }
    true
}
pub fn gamma_value(n: u32, t: u32) -> u32 {
    let g = gcd(n, t);
    gcd(g, n / g)
}
pub struct CriticalContentRows {
    midpoint: [u32; 2],
    next: [u32; 2],
    count: usize,
}
impl CriticalContentRows {
    pub fn new(plan: &CRTPlan) -> Self {
        let mut s = Self {
            midpoint: [0; 2],
            next: [0; 2],
            count: 0,
        };
        for i in 0..plan.count {
            let p = plan.prime[i];
            if p != 5 && p != 13 {
                continue;
            }
            let mut e = 0;
            let mut t = plan.power[i];
            while t > 1 {
                e += 1;
                t /= p;
            }
            if e % 2 != 0 || e == 0 {
                continue;
            }
            let mid = p.pow(e / 2);
            s.midpoint[s.count] = mid;
            s.next[s.count] = mid * p;
            s.count += 1;
        }
        s
    }
    pub fn mark(&self, permitted: &mut [u8], n: u32) {
        for i in 0..self.count {
            for t in (self.midpoint[i]..=n).step_by(self.midpoint[i] as usize) {
                if t % self.next[i] != 0 {
                    permitted[t as usize] = 0;
                }
            }
        }
    }
    pub fn accepts(&self, h: u32) -> bool {
        (0..self.count).all(|i| h % self.midpoint[i] != 0 || h % self.next[i] == 0)
    }
}
#[derive(Default)]
pub struct GammaPowerPlan {
    pub entries: Vec<(u32, u32, u32)>,
    pub relevant: u32,
}
impl GammaPowerPlan {
    pub fn new(plan: &CRTPlan) -> Self {
        let mut s = Self {
            entries: Vec::new(),
            relevant: 1,
        };
        for i in 0..plan.count {
            let (p, power) = (plan.prime[i], plan.power[i]);
            if p % 4 != 1 || power == p {
                continue;
            }
            let mut e = 0;
            let mut t = power;
            while t > 1 {
                e += 1;
                t /= p;
            }
            if (p == 5 || p == 13) && e == 2 {
                continue;
            }
            s.entries.push((p, power, e));
            s.relevant *= power;
        }
        s
    }
}
pub fn critical_mask(p: u32) -> u64 {
    let excluded: &[u32] = match p {
        17 => return 0x15554,
        29 => &[0, 1, 4, 8, 22, 26],
        37 => &[0, 1, 2, 4, 8, 12, 17, 21, 26, 30, 34, 36],
        41 => &[
            0, 1, 2, 3, 4, 5, 6, 9, 10, 16, 21, 26, 32, 33, 36, 37, 38, 39, 40,
        ],
        _ => return 0,
    };
    let mut mask = (1u64 << p) - 1;
    for &i in excluded {
        mask &= !(1u64 << i);
    }
    mask
}
pub struct CriticalProjectionPlan {
    entries: Vec<(u32, u32, u64)>,
}
impl CriticalProjectionPlan {
    pub fn new(plan: &CRTPlan) -> Self {
        let mut entries = Vec::new();
        for i in 0..plan.count {
            let p = plan.prime[i];
            let bits = critical_mask(p);
            if bits == 0 || plan.power[i] == p {
                continue;
            }
            let mut e = 0;
            let mut t = plan.power[i];
            while t > 1 {
                e += 1;
                t /= p;
            }
            if e % 2 != 0 {
                continue;
            }
            entries.push((p, p.pow(e / 2), bits));
        }
        Self { entries }
    }
    pub fn accepts(&self, u: u32, h: u32) -> bool {
        for &(p, mid, mask) in &self.entries {
            if h % mid != 0 || h % (mid * p) == 0 {
                continue;
            }
            if u % mid != 0 {
                return false;
            }
            let z =
                (((u / mid) % p) as u64 * inverse_mod((h / mid) % p, p) as u64 % p as u64) as u32;
            if (mask >> z) & 1 == 0 {
                return false;
            }
        }
        true
    }
}
