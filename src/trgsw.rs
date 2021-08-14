use super::trlwe::TRLWE;
use super::util::float_to_torus;
use super::util::ops::{intpoly_mul_as_torus, rmadd};
use super::util::params::trgsw;
use super::util::params::trlwe::Ring;

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const BG: u32 = trgsw::BG;
const BGBIT: u32 = trgsw::BGBIT;
const LBG: u32 = L as u32 * BGBIT;

type Z = trgsw::Z;
type Zpoly = [Z; N];
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
        let mut matrix: TRGSWMatrix = [[[0; N]; 2]; 2 * L];
        let zero_ring = [0; N];

        for l in 0..L {
            let w = 1 << (32 - (l + 1) * L);
            matrix[l][0] = intpoly_mul_as_torus(mu, w);
            matrix[l][1] = zero_ring;
            matrix[l + L][0] = zero_ring;
            matrix[l + L][1] = intpoly_mul_as_torus(mu, w);
        }

        rmadd(&matrix, &self.zero_matrix())
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
                println!("{} -> {}", r, s);
                assert!(s >= trgsw::SIGN_MIN);
            } else if s > trgsw::SIGN_MAX {
                println!("{} -> {}", r, s);
                assert!(s <= trgsw::SIGN_MAX);
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
    let mut decomp: Decomposition = ([[0; N]; L], [[0; N]; L]);
    for i in 0..L {
        decomp.0[i] = a_decomp[i];
        decomp.1[i] = b_decomp[i];
    }
    decomp
}

pub fn external_product((a_bar, b_bar): Decomposition, matrix: TRGSWMatrix) -> (Ring, Ring) {
    unimplemented!()
}

#[test]
fn test_decomposition() {
    use super::trlwe::TRLWE;
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    // (1) random BRing, (2) encrypt by TRLWE
    let b_poly: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(b_poly);

    // (3) decomposition a and b
    let a_bar = _decomposition(a);
    let b_bar = _decomposition(b);

    // (4) reconstruct a and b by decomposition
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
    for j in 0..N {
        let mut sa: i64 = 0;
        let mut sb: i64 = 0;
        for i in 0..L {
            sa += (a_bar[i][j] as i64) * (BG as i64).pow((L - i) as u32 - 1);
            sb += (b_bar[i][j] as i64) * (BG as i64).pow((L - i) as u32 - 1);
        }
        a_[j] = (sa.rem_euclid(2i64.pow(LBG)) << (32 - LBG)) as u32;
        b_[j] = (sb.rem_euclid(2i64.pow(LBG)) << (32 - LBG)) as u32;
    }

    // (5) decrypt by TRLWE, and assertion!
    assert_eq!(b_poly, trlwe.decrypt(a_, b_))
}

#[test]
fn test_external_product() {
    use super::trlwe::TRLWE;
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    let b_poly: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(b_poly);
    let decomp = decomposition(a, b);

    let trgsw = TRGSW::new(trlwe.get_secret());
    let matrix = trgsw.coefficient(1);

    let (a_, b_) = external_product(decomp, matrix);
    assert_eq!(b_poly, trlwe.decrypt(a_, b_))
}
