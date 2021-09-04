use super::trlwe::TRLWE;
use super::util::ops::{pmul, rmadd, vadd, vsub};
use super::util::params::{tlwe, trgsw};
use super::util::{float_to_torus, rotate_ring, zpoly_to_ring, RingLv0, RingLv1, Torus};

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const BG: u32 = trgsw::BG;
const BGBIT: u32 = trgsw::BGBIT;
const LBG: u32 = L as u32 * BGBIT;
const NBIT: usize = trgsw::NBIT;

type Z = i8;
type ZRing = [Z; N];
type Decomposition = ([ZRing; L], [ZRing; L]);
type TRGSWMatrix = [[RingLv1; 2]; 2 * L];
type BootstrappingKey = [TRGSWMatrix; tlwe::N];

fn intpoly_mul_as_torus<const N: usize>(zs: [i8; N], t: Torus) -> [Torus; N] {
    let mut ring = [0; N];
    for i in 0..N {
        let sc = zs[i];
        let base = t.wrapping_mul(sc.abs() as Torus);
        ring[i] = if sc < 0 {
            0u32.wrapping_sub(base)
        } else {
            base
        };
    }
    ring
}

pub struct TRGSW {
    s: RingLv1,
}

impl TRGSW {
    pub fn new(s: RingLv1) -> Self {
        Self { s }
    }

    fn zero_matrix(&self) -> TRGSWMatrix {
        let trlwe = TRLWE::from_secret(self.s.clone());
        let zero_ring = [float_to_torus(0.); N];
        let zero = [trlwe.encrypt_torus(zero_ring); 2 * L];

        let mut matrix: TRGSWMatrix = [[[0; N]; 2]; 2 * L];
        for i in 0..(2 * L) {
            let (a, b) = zero[i];
            matrix[i][0] = a;
            matrix[i][1] = b;
        }
        matrix
    }

    pub fn coefficient_matrix(&self, mu: [i8; N]) -> TRGSWMatrix {
        rmadd(&self._coefficient_matrix(mu), &self.zero_matrix())
    }

    fn _coefficient_matrix(&self, mu: [i8; N]) -> TRGSWMatrix {
        let mut matrix: TRGSWMatrix = [[[0; N]; 2]; 2 * L];
        let zero_ring = [0; N];
        for i in 0..L {
            let w = 2u32.pow(32 - (i as u32 + 1) * (BGBIT as u32));
            matrix[i][0] = intpoly_mul_as_torus(mu, w);
            matrix[i][1] = zero_ring;
            matrix[i + L][0] = zero_ring;
            matrix[i + L][1] = intpoly_mul_as_torus(mu, w);
        }
        matrix
    }

    pub fn coefficient(&self, m: i8) -> TRGSWMatrix {
        let mut mu = [0; N];
        mu[0] = m;
        self.coefficient_matrix(mu)
    }

    pub fn cmux(
        &self,
        c0: (RingLv1, RingLv1),
        c1: (RingLv1, RingLv1),
        flag: bool,
    ) -> (RingLv1, RingLv1) {
        let matrix = self.coefficient(flag as i8);
        cmux(matrix, c0, c1)
    }
}

fn _decomposition(poly: RingLv1) -> [ZRing; L] {
    let mut decomped: [[Z; N]; L] = [[0; N]; L];
    for n in 0..N {
        let mut a = poly[n] >> (32 - LBG);
        let mut cflag = false;
        for l in 0..L {
            let r = a % BG + cflag as u32;

            let s = if r >= (BG / 2) {
                cflag = true;
                (r - (BG / 2)) as Z - (BG / 2) as Z
            } else {
                cflag = false;
                r as Z
            };

            if s < trgsw::SIGN_MIN {
                assert!(s >= trgsw::SIGN_MIN, "{} -> {}", r, s);
            } else if s > trgsw::SIGN_MAX {
                assert!(s <= trgsw::SIGN_MAX, "{} -> {}", r, s);
            }

            decomped[L - l - 1][n] = s;
            a /= BG;
        }
    }
    decomped
}

pub fn cmux(
    matrix: TRGSWMatrix,
    (a0, b0): (RingLv1, RingLv1),
    (a1, b1): (RingLv1, RingLv1),
) -> (RingLv1, RingLv1) {
    let a_true = vsub(&a1, &a0);
    let b_true = vsub(&b1, &b0);
    let (enc_a, enc_b) = external_product(a_true, b_true, matrix);
    let a_ret = vadd(&enc_a, &a0);
    let b_ret = vadd(&enc_b, &b0);
    (a_ret, b_ret)
}

pub fn decomposition(a: RingLv1, b: RingLv1) -> Decomposition {
    let a_decomp = _decomposition(a);
    let b_decomp = _decomposition(b);
    (a_decomp, b_decomp)
}

pub fn external_product(a: RingLv1, b: RingLv1, matrix: TRGSWMatrix) -> (RingLv1, RingLv1) {
    let decomp = decomposition(a, b);
    _external_product(decomp, matrix)
}

pub fn _external_product((a_bar, b_bar): Decomposition, matrix: TRGSWMatrix) -> (RingLv1, RingLv1) {
    let mut a_: RingLv1 = [0; N];
    let mut b_: RingLv1 = [0; N];
    for i in 0..L {
        a_ = vadd(&a_, &pmul(&zpoly_to_ring(a_bar[i]), &matrix[i][0]));
        a_ = vadd(&a_, &pmul(&zpoly_to_ring(b_bar[i]), &matrix[i + L][0]));
        b_ = vadd(&b_, &pmul(&zpoly_to_ring(a_bar[i]), &matrix[i][1]));
        b_ = vadd(&b_, &pmul(&zpoly_to_ring(b_bar[i]), &matrix[i + L][1]));
    }
    (a_, b_)
}

pub fn blind_rotate(
    (a0, b0): (RingLv0, Torus),
    bk: BootstrappingKey,
    (a1, b1): (RingLv1, RingLv1),
) -> (RingLv1, RingLv1) {
    let b_floor = (b0 >> (31 - NBIT)) as usize;
    let offset = 2u32.pow(30 - NBIT as u32);

    let mut a_ret = rotate_ring(a1, 64 - b_floor);
    let mut b_ret = rotate_ring(b1, 64 - b_floor);

    for j in 0..tlwe::N {
        let a_floor = ((a0[j] + offset) >> (31 - NBIT)) as usize;
        let a_ret_rot = rotate_ring(a_ret, a_floor);
        let b_ret_rot = rotate_ring(b_ret, a_floor);
        let cmuxed = cmux(bk[j], (a_ret_rot, b_ret_rot), (a_ret, b_ret));
        a_ret = cmuxed.0;
        b_ret = cmuxed.1;
    }

    (a_ret, b_ret)
}

#[test]
fn test_decomposition() {
    use super::trlwe::TRLWE;
    use super::util::sampling::random_bool_initialization;

    // (1) random BRing, (2) encrypt by TRLWE
    let bs: [bool; N] = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);

    // (3) decomposition a and b
    let a_bar = _decomposition(a);
    let b_bar = _decomposition(b);

    // (4) reconstruct a and b from decomposition
    let mut a_: RingLv1 = [0; N];
    let mut b_: RingLv1 = [0; N];
    for j in 0..N {
        let mut sa: u32 = 0;
        let mut sb: u32 = 0;
        for i in 0..L {
            let w = 2u32.pow(32 - (i as u32 + 1) * (BGBIT as u32));
            sa = sa.wrapping_add((a_bar[i][j] as u32).wrapping_mul(w));
            sb = sb.wrapping_add((b_bar[i][j] as u32).wrapping_mul(w));
        }
        a_[j] = sa;
        b_[j] = sb;
    }

    // (5) decrypt by TRLWE, and assertion!
    assert_eq!(bs, trlwe.decrypt(a_, b_))
}

#[test]
fn test_zero_matrix_add() {
    use super::trlwe::TRLWE;
    use super::util::ops::vadd;
    use super::util::sampling::random_bool_initialization;

    let trlwe = TRLWE::new();
    let trgsw = TRGSW::new(trlwe.get_secret());

    let zm = trgsw.zero_matrix();

    for i in 0..zm.len() {
        let za = zm[i][0];
        let zb = zm[i][1];
        let bs: [bool; N] = random_bool_initialization();
        let (a, b) = trlwe.encrypt(bs);
        let dec_bs = trlwe.decrypt(vadd(&a, &za), vadd(&b, &zb));
        assert_eq!(bs, dec_bs);
    }
}

#[test]
fn test_zero_matrix_multiple() {
    use super::trlwe::TRLWE;
    use super::util::ops::vsub;
    use super::util::sampling::random_bool_initialization;

    let bs: [bool; N] = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);
    let decomp = decomposition(a, b);
    let trgsw = TRGSW::new(trlwe.get_secret());
    let matrix = trgsw.zero_matrix();
    let (a_, b_) = _external_product(decomp, matrix);
    let ring = trlwe.decrypt_torus(a_, b_);
    let offset = [2u32.pow(28); N];

    let m = vsub(&ring, &offset);
    let mut counter = 0;
    for i in 0..N {
        counter += (m[i] <= 2u32.pow(31)) as usize;
    }
    assert!(counter == 0, "counter is {}", counter);
}

#[test]
fn test_external_product() {
    use super::trlwe::TRLWE;
    use super::util::sampling::random_bool_initialization;

    let bs: [bool; N] = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);
    let decomp = decomposition(a, b);

    let trgsw = TRGSW::new(trlwe.get_secret());
    let matrix = trgsw.coefficient(1);

    let (a_, b_) = _external_product(decomp, matrix);
    let dec_bs = trlwe.decrypt(a_, b_);

    let mut counter = 0;
    for i in 0..N {
        counter += (bs[i] != dec_bs[i]) as usize;
    }

    assert!(counter == 0, "counter is {}", counter);
}

#[test]
fn test_cmux() {
    use super::trlwe::TRLWE;
    use super::util::sampling::random_bool_initialization;

    let trlwe = TRLWE::new();

    let bs1: [bool; N] = random_bool_initialization();
    let bs0: [bool; N] = random_bool_initialization();
    let enc0 = trlwe.encrypt(bs0);
    let enc1 = trlwe.encrypt(bs1);

    let trgsw = TRGSW::new(trlwe.get_secret());

    let (a, b) = trgsw.cmux(enc0, enc1, true);
    let dec_bs = trlwe.decrypt(a, b);
    let mut c1 = 0;
    for i in 0..N {
        c1 += (bs1[i] != dec_bs[i]) as usize;
    }

    let (a, b) = trgsw.cmux(enc0, enc1, false);
    let dec_bs = trlwe.decrypt(a, b);
    let mut c2 = 0;
    for i in 0..N {
        c2 += (bs0[i] != dec_bs[i]) as usize;
    }
    assert!(c1 == 0, "counter is {}", c1);
    assert!(c2 == 0, "counter is {}", c2);
}
