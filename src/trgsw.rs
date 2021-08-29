use super::trlwe::TRLWE;
use super::util::ops::{intpoly_mul_as_torus, rmadd};
use super::util::params::trgsw;
use super::util::params::trlwe::Ring;
use super::util::{float_to_torus, zpoly_to_ring};

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const BG: u32 = trgsw::BG;
const BGBIT: u32 = trgsw::BGBIT;
const LBG: u32 = L as u32 * BGBIT;

type Z = trgsw::Z;
type Zpoly = trgsw::Zpoly;
type Decomposition = ([Zpoly; L], [Zpoly; L]);
type TRGSWMatrix = [[Ring; 2]; 2 * L];

pub struct TRGSW {
    s: Ring,
}

impl TRGSW {
    pub fn new(s: Ring) -> Self {
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

    pub fn coefficient(self, m: i8) -> TRGSWMatrix {
        let mut mu = [0; N];
        mu[0] = m;
        self.coefficient_matrix(mu)
    }
}

fn _decomposition(poly: Ring) -> [Zpoly; L] {
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

pub fn decomposition(a: Ring, b: Ring) -> Decomposition {
    let a_decomp = _decomposition(a);
    let b_decomp = _decomposition(b);
    (a_decomp, b_decomp)
}

pub fn external_product((a_bar, b_bar): Decomposition, matrix: TRGSWMatrix) -> (Ring, Ring) {
    use super::util::ops::{pmul, vadd};
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
    for i in 0..L {
        a_ = vadd(&a_, &pmul(&zpoly_to_ring(a_bar[i]), &matrix[i][0]));
        a_ = vadd(&a_, &pmul(&zpoly_to_ring(b_bar[i]), &matrix[i + L][0]));
        b_ = vadd(&b_, &pmul(&zpoly_to_ring(a_bar[i]), &matrix[i][1]));
        b_ = vadd(&b_, &pmul(&zpoly_to_ring(b_bar[i]), &matrix[i + L][1]));
    }
    (a_, b_)
}

#[test]
fn test_decomposition() {
    use super::trlwe::TRLWE;
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    // (1) random BRing, (2) encrypt by TRLWE
    let bs: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);

    // (3) decomposition a and b
    let a_bar = _decomposition(a);
    let b_bar = _decomposition(b);

    // (4) reconstruct a and b from decomposition
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
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
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    let trlwe = TRLWE::new();
    let trgsw = TRGSW::new(trlwe.get_secret());

    let zm = trgsw.zero_matrix();

    for i in 0..zm.len() {
        let za = zm[i][0];
        let zb = zm[i][1];
        let bs: BRing = random_bool_initialization();
        let (a, b) = trlwe.encrypt(bs);
        let dec_bs = trlwe.decrypt(vadd(&a, &za), vadd(&b, &zb));
        assert_eq!(bs, dec_bs);
    }
}

#[test]
fn test_zero_matrix_multiple() {
    use super::trlwe::TRLWE;
    use super::util::ops::vsub;
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    let bs: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);
    let decomp = decomposition(a, b);
    let trgsw = TRGSW::new(trlwe.get_secret());
    let matrix = trgsw.zero_matrix();
    let (a_, b_) = external_product(decomp, matrix);
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
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    let bs: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);
    let decomp = decomposition(a, b);

    let trgsw = TRGSW::new(trlwe.get_secret());
    let matrix = trgsw.coefficient(1);

    let (a_, b_) = external_product(decomp, matrix);
    let dec_bs = trlwe.decrypt(a_, b_);

    let mut counter = 0;
    for i in 0..N {
        counter += (bs[i] != dec_bs[i]) as usize;
    }

    assert!(counter == 0, "counter is {}", counter);
}
