//! Exact limb-for-limb port of wide_arith_v5.hpp.
use crate::arithmetic::gcd;
use crate::Stats;
use std::cmp::Ordering;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct U256(pub [u64; 4]);
impl From<u128> for U256 {
    fn from(x: u128) -> Self {
        Self([x as u64, (x >> 64) as u64, 0, 0])
    }
}
impl Ord for U256 {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..4).rev() {
            match self.0[i].cmp(&other.0[i]) {
                Ordering::Equal => {}
                order => return order,
            }
        }
        Ordering::Equal
    }
}
impl PartialOrd for U256 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
pub fn add256(a: U256, b: U256) -> U256 {
    let mut r = U256::default();
    let mut carry = 0u128;
    for i in 0..4 {
        let s = a.0[i] as u128 + b.0[i] as u128 + carry;
        r.0[i] = s as u64;
        carry = s >> 64;
    }
    assert_eq!(carry, 0, "U256 add overflow");
    r
}
pub fn sub256(a: U256, b: U256) -> U256 {
    let mut r = U256::default();
    let mut borrow = 0u64;
    for i in 0..4 {
        let d = (a.0[i] as u128)
            .wrapping_sub(b.0[i] as u128)
            .wrapping_sub(borrow as u128);
        r.0[i] = d as u64;
        borrow = (d >> 64) as u64 & 1;
    }
    assert_eq!(borrow, 0, "U256 subtract underflow");
    r
}
pub fn shr256(a: U256, s: u32) -> U256 {
    let mut r = U256::default();
    for i in 0..4 {
        r.0[i] = a.0[i] >> s;
        if i < 3 {
            r.0[i] |= a.0[i + 1] << (64 - s);
        }
    }
    r
}
pub fn mul128(a: u128, b: u128) -> U256 {
    let (a0, a1, b0, b1) = (a as u64, (a >> 64) as u64, b as u64, (b >> 64) as u64);
    let p00 = a0 as u128 * b0 as u128;
    let p01 = a0 as u128 * b1 as u128;
    let p10 = a1 as u128 * b0 as u128;
    let p11 = a1 as u128 * b1 as u128;
    let mid = (p00 >> 64) + p01 as u64 as u128 + p10 as u64 as u128;
    let high = (mid >> 64) + (p01 >> 64) + (p10 >> 64) + p11 as u64 as u128;
    U256([
        p00 as u64,
        mid as u64,
        high as u64,
        ((high >> 64) + (p11 >> 64)) as u64,
    ])
}
pub fn magnitude128(x: i128) -> u128 {
    if x < 0 {
        0u128.wrapping_sub(x as u128)
    } else {
        x as u128
    }
}
pub fn square_i128(x: i128) -> U256 {
    let m = magnitude128(x);
    mul128(m, m)
}
pub fn mod256(a: U256, m: u32) -> u32 {
    let base = (1u128 << 64) % m as u128;
    let mut r = 0u128;
    for i in (0..4).rev() {
        r = (r * base + a.0[i] as u128 % m as u128) % m as u128;
    }
    r as u32
}
pub fn isqrt256(x: U256) -> u128 {
    if x == U256::default() {
        return 0;
    }
    let mut top = 255usize;
    while (x.0[top / 64] >> (top % 64)) & 1 == 0 {
        top -= 1;
    }
    let b = top & !1;
    let mut bit = U256::default();
    bit.0[b / 64] = 1 << (b % 64);
    let mut rem = x;
    let mut res = U256::default();
    loop {
        let trial = add256(res, bit);
        if trial <= rem {
            rem = sub256(rem, trial);
            res = add256(shr256(res, 1), bit);
        } else {
            res = shr256(res, 1);
        }
        if bit == U256::from(1) {
            break;
        }
        bit = shr256(bit, 2);
    }
    assert_eq!(
        (res.0[2], res.0[3]),
        (0, 0),
        "isqrt256 root exceeds 128 bits"
    );
    (res.0[1] as u128) << 64 | res.0[0] as u128
}
struct SquareResidues {
    q64: [u8; 64],
    q63: [u8; 63],
    q65: [u8; 65],
    q11: [u8; 11],
}
impl SquareResidues {
    fn new() -> Self {
        let mut s = Self {
            q64: [0; 64],
            q63: [0; 63],
            q65: [0; 65],
            q11: [0; 11],
        };
        for a in 0..64 {
            s.q64[a * a % 64] = 1;
        }
        for a in 0..63 {
            s.q63[a * a % 63] = 1;
        }
        for a in 0..65 {
            s.q65[a * a % 65] = 1;
        }
        for a in 0..11 {
            s.q11[a * a % 11] = 1;
        }
        s
    }
}
pub fn square256(z: U256) -> Option<u128> {
    static RESIDUES: OnceLock<SquareResidues> = OnceLock::new();
    let s = RESIDUES.get_or_init(SquareResidues::new);
    if s.q64[(z.0[0] & 63) as usize] == 0 {
        return None;
    }
    let r = mod256(z, 45045) as usize;
    if s.q63[r % 63] == 0 || s.q65[r % 65] == 0 || s.q11[r % 11] == 0 {
        return None;
    }
    let root = isqrt256(z);
    (mul128(root, root) == z).then_some(root)
}

pub const MAX_DENOMINATOR: u32 = 10_000_000;
pub fn evaluate_v5(
    n: u32,
    p: u32,
    q: u32,
    four: bool,
    st: &mut Stats,
    mut accept: impl FnMut(i128, i128, i128),
) {
    assert!(
        n >= 1 && n <= MAX_DENOMINATOR && q <= p && p <= n,
        "evaluate_v5 input outside proved range"
    );
    st.parameter_pairs += 1;
    let (n, p, q) = (n as i128, p as i128, q as i128);
    let (nn, pp, qq) = (n * n, p * p, q * q);
    let (u, v) = (nn - pp, nn - qq);
    let big_q = pp + qq - nn;
    let t = p * u + q * v;
    let r = u * u + v * v;
    let mut roots = [0i128; 2];
    let mut nr = 0;
    if big_q == 0 {
        st.linear_pairs += 1;
        if t == 0 {
            assert_ne!(r, 0, "identically zero quadratic");
            return;
        }
        if r % (4 * t) != 0 {
            return;
        }
        roots[nr] = r / (4 * t);
        nr += 1;
    } else {
        let delta = p - q;
        let w = 2 * nn - delta * delta;
        assert!(u >= 0 && v >= 0 && w > 0, "discriminant factor sign");
        let k = mul128(u as u128 * v as u128, w as u128);
        let Some(sr) = square256(k) else {
            return;
        };
        st.square_discriminants += 1;
        assert_eq!(sr >> 100, 0, "sqrt(K) beyond proved bound");
        let den = 2 * big_q;
        let (plus, minus) = (t + sr as i128, t - sr as i128);
        if plus % den == 0 {
            roots[nr] = plus / den;
            nr += 1;
        }
        if sr != 0 && minus % den == 0 {
            roots[nr] = minus / den;
            nr += 1;
        }
    }
    for &a in &roots[..nr] {
        if a < 0 {
            continue;
        }
        st.integral_a += 1;
        if four && a % 2 == 0 {
            continue;
        }
        let (xn, yn, den) = (u - 2 * a * p, v - 2 * a * q, 2 * n);
        if xn % den != 0 || yn % den != 0 {
            continue;
        }
        let (x, y) = (xn / den, yn / den);
        assert!(x <= y && 2 * y <= n, "symmetry failure");
        assert_eq!(
            add256(square_i128(x), square_i128(y)),
            square_i128(a),
            "three-distance identity failure"
        );
        assert_eq!(
            add256(square_i128(n - x), square_i128(y)),
            square_i128(a + p),
            "three-distance identity failure"
        );
        assert_eq!(
            add256(square_i128(x), square_i128(n - y)),
            square_i128(a + q),
            "three-distance identity failure"
        );
        assert!(a <= 2 * n * n * n, "radius bound failure");
        st.integer_coordinates += 1;
        if four {
            let gx = (magnitude128(x) % n as u128) as u32;
            let gy = (magnitude128(y) % n as u128) as u32;
            if gcd(gcd(n as u32, gx), gy) != 1 {
                continue;
            }
            let (xodd, yodd) = (x % 2 != 0, y % 2 != 0);
            if xodd == yodd {
                continue;
            }
            if (if xodd { y } else { x }) % 4 != 0 {
                continue;
            }
            if (x % 3 == 0) == (y % 3 == 0) {
                continue;
            }
            st.primitive_three_candidates += 1;
        }
        accept(x, y, a);
    }
}
