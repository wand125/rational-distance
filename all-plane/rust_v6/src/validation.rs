//! Independent checks from the frozen C++ sources. No C++ code is executed here.
use crate::arithmetic::*;
use crate::generator::*;
use crate::wide::*;
use crate::{Stats, Tables};
use std::collections::BTreeSet;

fn validation_search(
    n: u32,
    t: &Tables,
    u_first: bool,
    local: bool,
    difference: bool,
) -> (BTreeSet<[u32; 3]>, BTreeSet<[i128; 6]>) {
    let mut permitted = Vec::new();
    permitted_valuations(n, &t.arithmetic.plans[n as usize], &mut permitted);
    let mut parameters = BTreeSet::new();
    let mut points = BTreeSet::new();
    let mut ignored = Stats::default();
    let mut visit = |u: u32, v: u32| {
        if permitted[u as usize] == 0 || permitted[v as usize] == 0 {
            return;
        }
        if difference && !permitted_difference(n, u, v, &t.arithmetic) {
            return;
        }
        if local && !t.local.accepts(n, u, v) {
            return;
        }
        if !permitted_scale(n, u, v, gcd(n, u)) {
            return;
        }
        parameters.insert([2 * n, 2 * u, 2 * v]);
        evaluate_v5(2 * n, 2 * u, 2 * v, true, &mut ignored, |x, y, a| {
            points.insert([2 * n as i128, 2 * u as i128, 2 * v as i128, x, y, a]);
        });
    };
    if u_first {
        for u in 0..=n {
            if permitted[u as usize] == 0 || (local && !t.local.rows_possible(n, u)) {
                continue;
            }
            generate_crt(n, u, &t.arithmetic, |v| visit(u, v), gcd(n, u));
        }
    } else {
        let mut gap = |h| {
            generate_gap_crt(n, h, &t.arithmetic, |u| visit(u, u - h), gcd(n, h));
        };
        if local {
            t.gap.generate(n, gap);
        } else {
            for h in 0..=n {
                gap(h);
            }
        }
    }
    (parameters, points)
}
fn validate_gap(t: &Tables) {
    let mut all_crt_pairs = 0u64;
    for n in 1..=256 {
        let mut actual = Vec::new();
        let mut expected = Vec::new();
        for h in 0..=n {
            generate_gap_crt(n, h, &t.arithmetic, |u| actual.push([u, u - h]), 0);
        }
        actual.sort_unstable();
        for u in 0..=n {
            for v in 0..=u {
                let d = gcd(n / gcd(n, u), n / gcd(n, v));
                if (u - v) % d == 0 {
                    expected.push([u, v]);
                }
            }
        }
        assert_eq!(actual, expected, "gap CRT set/uniqueness mismatch");
        all_crt_pairs += actual.len() as u64;
    }
    let mut parameter_count = 0;
    let mut point_count = 0;
    for n in (12..=1000).step_by(12) {
        let reference = validation_search(n, t, true, true, true).0;
        let gap = validation_search(n, t, false, true, true).0;
        assert_eq!(
            reference, gap,
            "gap filtered parameter set mismatch with u-first"
        );
        parameter_count += gap.len();
        if n <= 500 {
            let reference = validation_search(n, t, true, false, false).1;
            let gap = validation_search(n, t, false, false, false).1;
            assert_eq!(
                reference, gap,
                "gap primitive three-distance positive examples mismatch"
            );
            point_count += gap.len();
        }
    }
    assert!(point_count > 0);
    assert!(
        permitted_scale(10140, 3900, 1235, gcd(10140, 3900)),
        "common-scale regression failure"
    );
    eprintln!("Gap validation passed: {all_crt_pairs} compatible pairs for every n=1..256; {parameter_count} fully filtered parameter tuples at N<=2000 exactly match u-first; {point_count} positive primitive three-distance examples retained.");
}
fn validate_field_wheel(t: &Tables) {
    let mut total = 0;
    for n in 1..=128 {
        for h in 0..=n {
            let mut actual = Vec::new();
            let mut expected = Vec::new();
            generate_field_gap_crt(n, h, &t.arithmetic, &t.wheel, |u| actual.push(u), 0);
            generate_gap_crt(
                n,
                h,
                &t.arithmetic,
                |u| {
                    if t.local.p19.accepts(n, u, u - h) {
                        expected.push(u);
                    }
                },
                0,
            );
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected, "progression wheel mismatch");
            total += actual.len();
        }
    }
    eprintln!("Wheel validation passed: {total} emitted CRT members.");
}
fn validate_v4_generators(t: &Tables) {
    let tables = ArithmeticTables::new(4913);
    let (mut rows, mut classified, mut wheel_members, mut lifted) = (0u64, 0u64, 0u64, 0u64);
    let mut verify_row = |n: u32, h: u32| {
        let plan = GammaPowerPlan::new(&tables.plans[n as usize]);
        let rule = GammaPowerRow::new(h, &plan);
        let gh = gamma_value(plan.relevant, h);
        let (mut expected, mut generated, mut refined) = (Vec::new(), Vec::new(), Vec::new());
        generate_gap_crt(
            n,
            h,
            &tables,
            |u| {
                let exact =
                    gamma_value(plan.relevant, u) == gh && gamma_value(plan.relevant, u - h) == gh;
                assert_eq!(
                    rule.accepts(u),
                    exact,
                    "gamma row differs from direct gcd definition"
                );
                classified += 1;
            },
            0,
        );
        generate_field_gap_crt(
            n,
            h,
            &tables,
            &t.wheel,
            |u| {
                assert_eq!(
                    t.local.accepts(n, u, u - h),
                    t.local.without19(n, u, u - h),
                    "removing redundant19 changed a wheel member"
                );
                if gamma_value(plan.relevant, u) == gh && gamma_value(plan.relevant, u - h) == gh {
                    expected.push(u);
                }
            },
            0,
        );
        let mut direct_stats = GenerationStats::default();
        let mut lift_stats = GenerationStats::default();
        generate_gamma_field_crt(
            n,
            h,
            &tables,
            &t.wheel,
            &rule,
            false,
            |u| generated.push(u),
            &mut direct_stats,
            false,
        );
        generate_gamma_field_crt(
            n,
            h,
            &tables,
            &t.wheel,
            &rule,
            true,
            |u| refined.push(u),
            &mut lift_stats,
            false,
        );
        expected.sort_unstable();
        generated.sort_unstable();
        refined.sort_unstable();
        assert_eq!(expected, generated, "gamma direct CRT/wheel set mismatch");
        assert_eq!(expected, refined, "gamma lifted CRT/wheel set mismatch");
        assert!(
            refined.windows(2).all(|w| w[0] != w[1]),
            "gamma CRT/wheel uniqueness mismatch"
        );
        let interior: Vec<u32> = expected
            .iter()
            .copied()
            .filter(|&u| u > h && u < n)
            .collect();
        let (mut range, mut trimmed_direct, mut trimmed_lift) =
            (Vec::new(), Vec::new(), Vec::new());
        generate_field_gap_crt_range(
            n,
            h,
            &tables,
            &t.wheel,
            |u| {
                if gamma_value(plan.relevant, u) == gh && gamma_value(plan.relevant, u - h) == gh {
                    range.push(u);
                }
            },
            gcd(n, h),
            h + 1,
            n - 1,
        );
        generate_gamma_field_crt(
            n,
            h,
            &tables,
            &t.wheel,
            &rule,
            false,
            |u| trimmed_direct.push(u),
            &mut direct_stats,
            true,
        );
        generate_gamma_field_crt(
            n,
            h,
            &tables,
            &t.wheel,
            &rule,
            true,
            |u| trimmed_lift.push(u),
            &mut lift_stats,
            true,
        );
        range.sort_unstable();
        trimmed_direct.sort_unstable();
        trimmed_lift.sort_unstable();
        assert_eq!(interior, range, "endpoint-restricted range mismatch");
        assert_eq!(
            interior, trimmed_direct,
            "endpoint-restricted direct mismatch"
        );
        assert_eq!(interior, trimmed_lift, "endpoint-restricted lift mismatch");
        rows += 1;
        wheel_members += expected.len() as u64;
        lifted += lift_stats.lifted_progressions;
    };
    for n in 1..=256 {
        let plan = GammaPowerPlan::new(&tables.plans[n as usize]);
        let labels = GammaLabels::new(n, &plan);
        for u in 0..=n {
            if !labels.value.is_empty() {
                assert_eq!(
                    labels.value[u as usize],
                    gamma_value(plan.relevant, u),
                    "gamma label differs from direct gcd definition"
                );
            }
        }
        for h in 0..=n {
            verify_row(n, h);
        }
    }
    for n in [625u32, 1250, 2197, 3375, 4913] {
        let mut gaps = BTreeSet::from([0, 1, n - 1, n]);
        let f = &tables.plans[n as usize];
        for i in 0..f.count {
            let mut q = f.prime[i];
            while q <= f.power[i] {
                gaps.insert(q);
                if q < n {
                    gaps.insert(q + 1);
                }
                if 2 * q <= n {
                    gaps.insert(2 * q);
                }
                if q == f.power[i] {
                    break;
                }
                q *= f.prime[i];
            }
        }
        for h in gaps {
            verify_row(n, h);
        }
    }
    let plan = GammaPowerPlan::new(&tables.plans[125]);
    let row = GammaPowerRow::new(25, &plan);
    assert!(
        row.accepts(30) && row.accepts(50) && !row.accepts(25),
        "high-k gamma branch regression"
    );
    eprintln!("V4 validation passed: {rows} gap rows, {classified} direct gcd comparisons, {wheel_members} retained CRT/wheel members, {lifted} lifted progressions; both high-k branches retained.");
}

// std::mt19937_64, including the reference's deterministic seeds and draw order.
struct Mt19937 {
    state: [u64; 312],
    index: usize,
}
impl Mt19937 {
    fn new(seed: u64) -> Self {
        let mut state = [0u64; 312];
        state[0] = seed;
        for i in 1..312 {
            state[i] = 6364136223846793005u64
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 62))
                .wrapping_add(i as u64);
        }
        Self { state, index: 312 }
    }
    fn next(&mut self) -> u64 {
        if self.index >= 312 {
            for i in 0..312 {
                let x =
                    (self.state[i] & 0xffffffff80000000) | (self.state[(i + 1) % 312] & 0x7fffffff);
                self.state[i] = self.state[(i + 156) % 312]
                    ^ (x >> 1)
                    ^ if x & 1 != 0 { 0xb5026f5aa96619e9 } else { 0 };
            }
            self.index = 0;
        }
        let mut x = self.state[self.index];
        self.index += 1;
        x ^= (x >> 29) & 0x5555555555555555;
        x ^= (x << 17) & 0x71d67fffeda60000;
        x ^= (x << 37) & 0xfff7eee000000000;
        x ^= x >> 43;
        x
    }
}
fn validate_wide_arithmetic() {
    let mut rng = Mt19937::new(20261105);
    let mut checks = 0u64;
    let mut check_root = |r: u128| {
        let s = mul128(r, r);
        assert_eq!(square256(s), Some(r), "wide square of root failed");
        assert_eq!(isqrt256(s), r);
        if r > 0 {
            let below = sub256(s, U256::from(1));
            assert_eq!(isqrt256(below), r - 1, "wide r^2-1 failed");
            if r > 1 {
                assert!(square256(below).is_none());
            }
        }
        let above = add256(s, U256::from(2 * r));
        assert_eq!(isqrt256(above), r, "wide r^2+2r root failed");
        if r > 0 {
            assert!(square256(above).is_none());
        }
        checks += 4;
    };
    for r in 0..5000 {
        check_root(r);
    }
    for bits in 1..=127 {
        for _ in 0..200 {
            let mut r = ((rng.next() as u128) << 64) | rng.next() as u128;
            r &= (1u128 << bits) - 1;
            check_root(r);
            check_root(1u128 << (bits - 1));
            check_root((1u128 << bits) - 1);
        }
    }
    for _ in 0..200000 {
        let (a, b) = (rng.next() as u128, rng.next() as u128);
        assert_eq!(mul128(a, b), U256::from(a * b), "mul128 small product");
        let c = ((rng.next() as u128) << 40) | ((rng.next() as u128) >> 24);
        let d = (rng.next() as u128) >> 16;
        let left = mul128(c, d);
        let high = mul128(c >> 64, d);
        assert_eq!(high.0[3], 0);
        let shifted = U256([0, high.0[0], high.0[1], high.0[2]]);
        let right = add256(mul128(c as u64 as u128, d), shifted);
        assert_eq!(left, right, "mul128 limb decomposition");
        checks += 2;
    }
    for t in 0..20000 {
        let x = ((((rng.next() as u128) << 60) ^ rng.next() as u128) as i128)
            * if t & 1 != 0 { -1 } else { 1 };
        let small = (x >> 70) as i64;
        assert_eq!(
            small.to_string().parse::<i64>().unwrap(),
            small,
            "i128 print"
        );
        assert_eq!(x.to_string().parse::<i128>().unwrap(), x, "i128 full print");
    }
    // Also exercise the top limb and the full unsigned root range.
    for r in [1u128 << 127, u128::MAX - 1, u128::MAX] {
        let s = mul128(r, r);
        assert_eq!(isqrt256(s), r);
        assert_eq!(square256(s), Some(r));
        assert_eq!(isqrt256(sub256(s, U256::from(1))), r - 1);
    }
    eprintln!("V5 wide arithmetic validation passed: {checks} exact checks.");
}

// Independent native-integer square test for the coordinate brute force.
// Unlike the frozen core's floating estimate this uses integer Newton division.
fn small_square(z: u64) -> bool {
    if (0b10011u32 >> (z % 5)) & 1 == 0 || (0b10111u32 >> (z % 7)) & 1 == 0 {
        return false;
    }
    if z == 0 {
        return true;
    }
    let mut r = 1u64 << ((64 - z.leading_zeros() + 1) / 2);
    loop {
        let next = (r + z / r) / 2;
        if next >= r {
            return r * r == z;
        }
        r = next;
    }
}
fn evaluate_points(n: u32, p: u32, q: u32) -> Vec<[i128; 3]> {
    let mut points = Vec::new();
    evaluate_v5(n, p, q, false, &mut Stats::default(), |x, y, a| {
        points.push([x, y, a])
    });
    points
}
fn validate_v5_reconstruction() {
    let mut total = 0;
    let mut seeds = Vec::new();
    for n in 1u32..=12 {
        let mut algebraic = BTreeSet::new();
        let mut brute = BTreeSet::new();
        for p in 0..=n {
            for q in 0..=p {
                for point in evaluate_points(n, p, q) {
                    algebraic.insert([point[0], point[1]]);
                    seeds.push([n, p, q]);
                }
            }
        }
        let n = n as i64;
        let radius = 2 * n * n * n;
        for y in -radius..=n / 2 {
            for x in -radius..=y {
                if !small_square((x * x + y * y) as u64) {
                    continue;
                }
                if !small_square(((n - x) * (n - x) + y * y) as u64) {
                    continue;
                }
                if !small_square((x * x + (n - y) * (n - y)) as u64) {
                    continue;
                }
                brute.insert([x as i128, y as i128]);
            }
        }
        assert_eq!(algebraic, brute, "v5 brute validation mismatch at N={n}");
        total += algebraic.len();
    }
    eprintln!("V5 brute reconstruction validation passed: N=1..12 ({total} points).");
    seeds.push([8400, 6720, 5040]);
    assert_eq!(
        evaluate_points(8400, 6720, 5040),
        vec![[164, 1677, 1685]],
        "linear-case radius failure"
    );
    assert!(
        evaluate_points(24, 24, 0).is_empty(),
        "Q=T=0 validation failure"
    );
    let mut scaled = 0;
    for [n, p, q] in seeds {
        let base = evaluate_points(n, p, q);
        for k in [
            1,
            MAX_DENOMINATOR / n,
            MAX_DENOMINATOR / n - 1,
            MAX_DENOMINATOR / n / 2 + 1,
        ] {
            if k < 1 || k * n > MAX_DENOMINATOR {
                continue;
            }
            let mut got = evaluate_points(k * n, k * p, k * q);
            let mut expected: Vec<[i128; 3]> = base
                .iter()
                .map(|b| [k as i128 * b[0], k as i128 * b[1], k as i128 * b[2]])
                .collect();
            got.sort_unstable();
            expected.sort_unstable();
            assert_eq!(
                got,
                expected,
                "scaled three-distance point lost at N={}",
                k * n
            );
            scaled += 1;
        }
    }
    let n = MAX_DENOMINATOR;
    for p in [0, 1, n / 2, n - 1, n] {
        for q in [0, 1, n / 3, n - 1, n] {
            if q <= p {
                evaluate_points(n, p, q);
            }
        }
    }
    eprintln!("V5 scaled-seed and extreme-parameter validation passed: {scaled} scaled seed families up to N=10^7 (arithmetic only).");
}
pub(super) fn validate(t: &Tables) {
    validate_gap(t);
    validate_field_wheel(t);
    validate_v4_generators(t);
    if let Some(additional) = &t.additional {
        additional.verify();
    }
    t.local.verify();
    validate_wide_arithmetic();
    validate_v5_reconstruction();
}
