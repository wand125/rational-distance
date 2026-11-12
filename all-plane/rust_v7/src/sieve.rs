//! Finite-field and lifted projections, fifth-adic marking, critical lifts.
use crate::arithmetic::inverse_mod;

fn normal_projection<const M: usize>() -> (Vec<u8>, Vec<u8>) {
    let m = M as u32;
    let mut roots: Vec<Vec<u32>> = vec![Vec::new(); M];
    let mut qr = vec![0; M];
    for a in 0..m {
        roots[(a * a % m) as usize].push(a);
        qr[(a * a % m) as usize] = 1;
    }
    let mut normal = vec![0; M * M];
    for x in 0..m {
        for y in 0..m {
            let (xx, yy) = ((1 + m - x) % m, (1 + m - y) % m);
            let (aa, bb, cc, dd) = (
                (x * x + y * y) % m,
                (xx * xx + y * y) % m,
                (x * x + yy * yy) % m,
                (xx * xx + yy * yy) % m,
            );
            if qr[dd as usize] == 0 {
                continue;
            }
            for &a in &roots[aa as usize] {
                for &b in &roots[bb as usize] {
                    for &c in &roots[cc as usize] {
                        normal[(((b + m - a) % m) * m + (c + m - a) % m) as usize] = 1;
                    }
                }
            }
        }
    }
    (normal, qr)
}
pub struct PrimeProjection<const P: usize> {
    cube: Vec<u8>,
    any: Vec<u8>,
    qr: Vec<u8>,
}
impl<const P: usize> PrimeProjection<P> {
    pub fn new() -> Self {
        let (normal, qr) = normal_projection::<P>();
        let mut s = Self {
            cube: vec![0; P * P * P],
            any: vec![0; P * P],
            qr,
        };
        let p = P as u32;
        for n in 0..p {
            let inv = if n != 0 { inverse_mod(n, p) } else { 0 };
            for u in 0..p {
                for v in 0..p {
                    let okay = if n != 0 {
                        normal[((u * inv % p) * p + v * inv % p) as usize]
                    } else {
                        u8::from(u * v * (u + p - v) % p == 0)
                    };
                    s.cube[((n * p + u) * p + v) as usize] = okay;
                    s.any[(n * p + u) as usize] |= okay;
                }
            }
        }
        s
    }
    #[inline]
    pub fn accepts(&self, n: u32, u: u32, v: u32) -> bool {
        let p = P as u32;
        self.cube[(((n % p) * p + u % p) * p + v % p) as usize] != 0
    }
    pub fn row_exists(&self, n: u32, u: u32) -> bool {
        let p = P as u32;
        self.any[((n % p) * p + u % p) as usize] != 0
    }
    pub fn verify(&self) {
        for n in 0..P as u32 {
            for u in 0..P as u32 {
                for v in 0..P as u32 {
                    let mut expected = false;
                    let inv = if n != 0 {
                        inverse_mod(2 * n % P as u32, P as u32)
                    } else {
                        0
                    };
                    for a in 0..P as u32 {
                        if expected {
                            break;
                        }
                        if n == 0 {
                            expected =
                                u * (2 * a + u) % P as u32 == 0 && v * (2 * a + v) % P as u32 == 0;
                            continue;
                        }
                        let x = ((n * n) as i32 - (u * u) as i32 - (2 * a * u) as i32)
                            .rem_euclid(P as i32) as u32
                            * inv
                            % P as u32;
                        let y = ((n * n) as i32 - (v * v) as i32 - (2 * a * v) as i32)
                            .rem_euclid(P as i32) as u32
                            * inv
                            % P as u32;
                        let (xx, yy) =
                            ((n + P as u32 - x) % P as u32, (n + P as u32 - y) % P as u32);
                        expected = (x * x + y * y) % P as u32 == a * a % P as u32
                            && self.qr[((xx * xx + yy * yy) % P as u32) as usize] != 0;
                    }
                    assert_eq!(
                        expected,
                        self.accepts(n, u, v),
                        "finite-field projection mismatch mod {P}"
                    );
                }
            }
        }
        eprintln!("PrimeProjection<{P}> validation passed.");
    }
}
pub struct LocalSieve {
    pub p5: PrimeProjection<5>,
    pub p7: PrimeProjection<7>,
    pub p11: PrimeProjection<11>,
    pub p13: PrimeProjection<13>,
    pub p17: PrimeProjection<17>,
    pub p19: PrimeProjection<19>,
    pub p23: PrimeProjection<23>,
    pub p29: PrimeProjection<29>,
    pub p31: PrimeProjection<31>,
    pub p37: PrimeProjection<37>,
    pub p41: PrimeProjection<41>,
    pub p43: PrimeProjection<43>,
}
impl LocalSieve {
    pub fn new() -> Self {
        Self {
            p5: PrimeProjection::new(),
            p7: PrimeProjection::new(),
            p11: PrimeProjection::new(),
            p13: PrimeProjection::new(),
            p17: PrimeProjection::new(),
            p19: PrimeProjection::new(),
            p23: PrimeProjection::new(),
            p29: PrimeProjection::new(),
            p31: PrimeProjection::new(),
            p37: PrimeProjection::new(),
            p41: PrimeProjection::new(),
            p43: PrimeProjection::new(),
        }
    }
    #[inline]
    pub fn without19(&self, n: u32, u: u32, v: u32) -> bool {
        self.p5.accepts(n, u, v)
            && self.p11.accepts(n, u, v)
            && self.p7.accepts(n, u, v)
            && self.p29.accepts(n, u, v)
            && self.p37.accepts(n, u, v)
            && self.p13.accepts(n, u, v)
            && self.p17.accepts(n, u, v)
            && self.p23.accepts(n, u, v)
            && self.p31.accepts(n, u, v)
            && self.p41.accepts(n, u, v)
            && self.p43.accepts(n, u, v)
    }
    pub fn accepts(&self, n: u32, u: u32, v: u32) -> bool {
        self.p19.accepts(n, u, v) && self.without19(n, u, v)
    }
    pub fn rows_possible(&self, n: u32, u: u32) -> bool {
        self.p11.row_exists(n, u) && self.p19.row_exists(n, u)
    }
    pub fn verify(&self) {
        self.p5.verify();
        self.p7.verify();
        self.p11.verify();
        self.p13.verify();
        self.p17.verify();
        self.p19.verify();
        self.p23.verify();
        self.p29.verify();
        self.p31.verify();
        self.p37.verify();
        self.p41.verify();
        self.p43.verify();
        eprintln!("LocalSieve validation passed.");
    }
}
pub struct LiftedProjection<const M: usize, const P: usize> {
    normal: Vec<u8>,
    inverse: Vec<u16>,
    qr: Vec<u8>,
}
impl<const M: usize, const P: usize> LiftedProjection<M, P> {
    pub fn new() -> Self {
        assert!(M % P == 0 && M % 2 == 1);
        let (normal, qr) = normal_projection::<M>();
        let mut inverse = vec![0; M];
        for (a, inv) in inverse.iter_mut().enumerate() {
            if a % P != 0 {
                *inv = inverse_mod(a as u32, M as u32) as u16;
            }
        }
        Self {
            normal,
            inverse,
            qr,
        }
    }
    #[inline]
    pub fn accepts(&self, n: u32, u: u32, v: u32) -> bool {
        let m = M as u32;
        let inv = self.inverse[(n % m) as usize] as u32;
        inv == 0 || self.normal[((u % m * inv % m) * m + v % m * inv % m) as usize] != 0
    }
    pub fn verify(&self) {
        let inv2 = (M as u32 + 1) / 2;
        for u in 0..M as u32 {
            for v in 0..M as u32 {
                let mut exists = false;
                for a in 0..M as u32 {
                    if exists {
                        break;
                    }
                    let x = (1 - (u * u) as i32 - (2 * a * u) as i32).rem_euclid(M as i32) as u32
                        * inv2
                        % M as u32;
                    let y = (1 - (v * v) as i32 - (2 * a * v) as i32).rem_euclid(M as i32) as u32
                        * inv2
                        % M as u32;
                    let (xx, yy) = ((1 + M as u32 - x) % M as u32, (1 + M as u32 - y) % M as u32);
                    exists = (x * x + y * y) % M as u32 == a * a % M as u32
                        && self.qr[((xx * xx + yy * yy) % M as u32) as usize] != 0;
                }
                assert_eq!(
                    exists,
                    self.normal[u as usize * M + v as usize] != 0,
                    "lifted projection differs from radius reconstruction mod {M}"
                );
            }
        }
        eprintln!("LiftedProjection<{M},{P}> validation passed.");
    }
}
pub struct FifthPowerCondition;
impl FifthPowerCondition {
    pub fn accepts_one(n: u32, u: u32) -> bool {
        if n % 5 == 0 || (u % 5 != n % 5 && u % 5 != (5 - n % 5) % 5) {
            return true;
        }
        let mut value = (n - u) as u64 * (n + u) as u64;
        if value == 0 {
            return true;
        }
        let mut e = 0;
        while value % 5 == 0 {
            value /= 5;
            e += 1;
        }
        e % 2 == 0 && (value % 5 == 1 || value % 5 == 4)
    }
    #[allow(dead_code)]
    pub fn accepts(n: u32, u: u32, v: u32) -> bool {
        Self::accepts_one(n, u) && Self::accepts_one(n, v)
    }
    #[allow(dead_code)]
    pub fn mark(n: u32, permitted: &mut [u8]) {
        if n % 5 == 0 {
            return;
        }
        for start in [n % 5, (5 - n % 5) % 5] {
            for u in (start..=n).step_by(5) {
                if !Self::accepts_one(n, u) {
                    permitted[u as usize] = 0;
                }
            }
        }
    }
}
pub struct FifthFastMarker {
    forbidden: [[u8; 44]; 125],
}
impl FifthFastMarker {
    pub fn new() -> Self {
        let mut s = Self {
            forbidden: [[0; 44]; 125],
        };
        for r in 1u32..125 {
            if r % 5 == 0 {
                continue;
            }
            let mut count = 0;
            for u in 0..125 {
                if u % 5 != r % 5 && u % 5 != (5 - r % 5) % 5 {
                    continue;
                }
                let mut w = (r * r + 15625 - u * u) % 125;
                if w == 0 {
                    continue;
                }
                let mut e = 0;
                while w % 5 == 0 {
                    w /= 5;
                    e += 1;
                }
                if e % 2 != 0 || (w % 5 != 1 && w % 5 != 4) {
                    assert!(count < 44);
                    s.forbidden[r as usize][count] = u as u8;
                    count += 1;
                }
            }
            assert_eq!(count, 44, "fifth prefix cardinality");
        }
        s
    }
    pub fn mark(&self, n: u32, permitted: &mut [u8]) {
        if n % 5 == 0 {
            return;
        }
        let r = n % 125;
        for &residue in &self.forbidden[r as usize] {
            for u in (residue as u32..=n).step_by(125) {
                permitted[u as usize] = 0;
            }
        }
        for residue in [r, (125 - r) % 125] {
            for u in (residue..=n).step_by(125) {
                if !FifthPowerCondition::accepts_one(n, u) {
                    permitted[u as usize] = 0;
                }
            }
        }
    }
    pub fn verify(&self) {
        let check = |n| {
            let mut marked = vec![1; n as usize + 1];
            self.mark(n, &mut marked);
            for u in 0..=n {
                assert_eq!(
                    marked[u as usize] != 0,
                    FifthPowerCondition::accepts_one(n, u),
                    "fast fifth mark differs from exact Q_5 definition"
                );
            }
        };
        for n in 1..=1000 {
            check(n);
        }
        for n in [15624, 78124, 999984, 999996, 1000000] {
            check(n);
        }
        eprintln!("FifthFastMarker validation passed.");
    }
}
pub struct PrimeData {
    prime: u32,
    allowed: [u8; 282],
    lifts: Vec<(u32, u32)>,
}
pub struct CriticalExtension {
    tables: Vec<PrimeData>,
}
impl CriticalExtension {
    fn quartic(p: u32, iota: u32, r: u32) -> u32 {
        let t = (r as u64 * (r + p - 1) as u64 % p as u64) as u32;
        let left =
            ((2 * r + p - 1) as u64 % p as u64 * ((2 * t + p - 1) % p) as u64 % p as u64) as u32;
        let right = (2u64 * iota as u64 % p as u64 * t as u64 % p as u64 * (t + 2) as u64
            % p as u64) as u32;
        (left + p - right) % p
    }
    pub fn new() -> Self {
        let mut tables = Vec::new();
        for p in (17u32..=281).step_by(4) {
            let mut prime = true;
            let mut d = 3;
            while d * d <= p {
                if p % d == 0 {
                    prime = false;
                    break;
                }
                d += 2;
            }
            if !prime {
                continue;
            }
            let mut table = PrimeData {
                prime: p,
                allowed: [0; 282],
                lifts: Vec::new(),
            };
            let mut squares = [0u8; 282];
            let mut iota = 0;
            for a in 0..p {
                let sq = a * a % p;
                squares[sq as usize] = 1;
                if sq == p - 1 && iota == 0 {
                    iota = a;
                }
            }
            for r in 2..p {
                let denominator = (r as u64 + p as u64 - 1 + p as u64
                    - iota as u64 * r as u64 % p as u64)
                    % p as u64;
                if denominator != 0 && squares[Self::quartic(p, iota, r) as usize] != 0 {
                    table.allowed[r as usize] = 1;
                    table.allowed[((1 + p - r) % p) as usize] = 1;
                }
            }
            table.lifts = match p {
                41 => vec![
                    (7, 48),
                    (12, 750),
                    (15, 1327),
                    (27, 355),
                    (30, 932),
                    (35, 1634),
                ],
                113 => vec![(9, 10518), (50, 3892), (64, 8878), (105, 2252)],
                137 => vec![(43, 15113), (55, 9234), (83, 9536), (95, 3657)],
                257 => vec![
                    (23, 36260),
                    (62, 44523),
                    (117, 3715),
                    (141, 62335),
                    (196, 21527),
                    (235, 29790),
                ],
                _ => Vec::new(),
            };
            tables.push(table);
        }
        Self { tables }
    }
    // Active per-n state belongs to the worker; the tables are shared.
    pub fn prepare(&self, n: u32, active: &mut Vec<(usize, u32)>) {
        active.clear();
        if n == 0 {
            return;
        }
        for (index, table) in self.tables.iter().enumerate() {
            let p = table.prime;
            if n % (p * p) != 0 {
                continue;
            }
            let mut remaining = n;
            let mut e = 0;
            loop {
                remaining /= p;
                e += 1;
                if remaining % p != 0 {
                    break;
                }
            }
            if e & 1 != 0 {
                continue;
            }
            active.push((index, p.pow(e / 2)));
        }
    }
    pub fn accepts(&self, u: u32, h: u32, active: &[(usize, u32)]) -> bool {
        for &(index, mid) in active {
            let table = &self.tables[index];
            let p = table.prime;
            if h % mid != 0 || (h / mid) % p == 0 {
                continue;
            }
            if u % mid != 0 {
                return false;
            }
            let (big_u, big_h) = (u / mid, h / mid);
            let residue = ((big_u % p) as u64 * inverse_mod(big_h % p, p) as u64 % p as u64) as u32;
            if table.allowed[residue as usize] == 0 {
                return false;
            }
            for &(r, root) in &table.lifts {
                if residue == r {
                    let modulus = p * p;
                    if (big_u % modulus) as u64
                        != root as u64 * (big_h % modulus) as u64 % modulus as u64
                    {
                        return false;
                    }
                    break;
                }
            }
        }
        true
    }
}
pub struct AdditionalSieve {
    p49: LiftedProjection<49, 7>,
    p121: LiftedProjection<121, 11>,
    p169: LiftedProjection<169, 13>,
    p125: LiftedProjection<125, 5>,
    fifth_marker: FifthFastMarker,
    critical: CriticalExtension,
}
impl AdditionalSieve {
    pub fn new() -> Self {
        Self {
            p49: LiftedProjection::new(),
            p121: LiftedProjection::new(),
            p169: LiftedProjection::new(),
            p125: LiftedProjection::new(),
            fifth_marker: FifthFastMarker::new(),
            critical: CriticalExtension::new(),
        }
    }
    pub fn mark(&self, n: u32, permitted: &mut [u8], active: &mut Vec<(usize, u32)>) {
        self.fifth_marker.mark(n, permitted);
        self.critical.prepare(n, active);
    }
    #[inline]
    pub fn accepts(&self, n: u32, u: u32, v: u32, active: &[(usize, u32)]) -> bool {
        self.p125.accepts(n, u, v)
            && self.p169.accepts(n, u, v)
            && self.p121.accepts(n, u, v)
            && self.p49.accepts(n, u, v)
            && self.critical.accepts(u, u - v, active)
    }
    pub fn verify(&self) {
        self.p49.verify();
        self.p121.verify();
        self.p169.verify();
        self.p125.verify();
        self.fifth_marker.verify();
        eprintln!("AdditionalSieve validation passed.");
    }
}
