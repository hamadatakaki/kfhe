use super::key::SecretKey;
use super::ops::{pmul, rmadd, vadd};
use super::params::{tlwe, trgsw};
use super::tlwe::{CipherTLWELv0, CipherTLWELv1};
use super::trlwe::{sample_extract_index, CipherTRLWE, TRLWE};
use super::util::{float_to_torus, rotate_ring, zpoly_to_ring, RingLv1, Torus};

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const BG: u32 = trgsw::BG;
const BGBIT: u32 = trgsw::BGBIT;
const LBG: u32 = L as u32 * BGBIT;
const NBIT: usize = trgsw::NBIT;

pub type Z = i8;
pub type TRGSWMatrix = [[RingLv1; 2]; 2 * L];

type ZRing = [Z; N];
type Decomposition = ([ZRing; L], [ZRing; L]);

#[derive(Clone, Debug)]
pub struct BootstrappingKey(pub Vec<TRGSWMatrix>);

impl BootstrappingKey {
    pub fn get(&self, k: usize) -> TRGSWMatrix {
        self.0[k]
    }

    pub fn set(&mut self, k: usize, val: TRGSWMatrix) {
        self.0[k] = val;
    }
}

pub fn uninitialized_trgsw_matrix() -> TRGSWMatrix {
    [[[0; N]; 2]; 2 * L]
}

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

#[derive(Clone, Copy, Debug)]
pub struct TRGSW {
    sk: SecretKey,
}

impl TRGSW {
    pub fn new(sk: SecretKey) -> Self {
        Self { sk }
    }

    fn zero_matrix(&self) -> TRGSWMatrix {
        let trlwe = TRLWE::new(self.sk);
        let zero_ring = [float_to_torus(0.); N];
        let zero = [trlwe.encrypt_torus(zero_ring).describe(); 2 * L];

        let mut matrix: TRGSWMatrix = [[[0; N]; 2]; 2 * L];
        for i in 0..(2 * L) {
            let (a, b) = zero[i];
            matrix[i][0] = a;
            matrix[i][1] = b;
        }
        matrix
    }

    pub fn coefficient_matrix(&self, mu: [i8; N]) -> TRGSWMatrix {
        rmadd(self._coefficient_matrix(mu), self.zero_matrix())
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

    pub fn coefficient_bool(&self, b: bool) -> TRGSWMatrix {
        self.coefficient(b as i8)
    }

    pub fn cmux(&self, flag: bool, c0: CipherTRLWE, c1: CipherTRLWE) -> CipherTRLWE {
        let matrix = self.coefficient(flag as i8);
        cmux(matrix, c0, c1)
    }

    pub fn bootstrapping_key(&self) -> BootstrappingKey {
        let n = tlwe::N;
        let lv0 = self.sk.lv0;
        let mut bk = BootstrappingKey(vec![[[[0; N]; 2]; 2 * L]; tlwe::N]);

        for j in 0..n {
            bk.set(j, self.coefficient(lv0[j] as i8));
        }
        bk
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

pub fn decomposition(c: CipherTRLWE) -> Decomposition {
    let (a, b) = c.describe();
    let a_decomp = _decomposition(a);
    let b_decomp = _decomposition(b);
    (a_decomp, b_decomp)
}

pub fn external_product(matrix: TRGSWMatrix, c: CipherTRLWE) -> CipherTRLWE {
    let decomp = decomposition(c);
    let (a, b) = _external_product(matrix, decomp);
    CipherTRLWE(a, b)
}

fn _external_product(matrix: TRGSWMatrix, (a_bar, b_bar): Decomposition) -> (RingLv1, RingLv1) {
    let mut a_: RingLv1 = [0; N];
    let mut b_: RingLv1 = [0; N];
    for i in 0..L {
        a_ = vadd(a_, pmul(zpoly_to_ring(a_bar[i]), matrix[i][0]));
        a_ = vadd(a_, pmul(zpoly_to_ring(b_bar[i]), matrix[i + L][0]));
        b_ = vadd(b_, pmul(zpoly_to_ring(a_bar[i]), matrix[i][1]));
        b_ = vadd(b_, pmul(zpoly_to_ring(b_bar[i]), matrix[i + L][1]));
    }
    (a_, b_)
}

pub fn cmux(matrix: TRGSWMatrix, c0: CipherTRLWE, c1: CipherTRLWE) -> CipherTRLWE {
    external_product(matrix, c0 - c1) + c1
}

fn rotate_trlwe_cipher(c: CipherTRLWE, k: usize) -> CipherTRLWE {
    let (a, b) = c.describe();
    let a_rot = rotate_ring(a, k);
    let b_rot = rotate_ring(b, k);
    CipherTRLWE(a_rot, b_rot)
}

pub fn blind_rotate(c0: CipherTLWELv0, bk: BootstrappingKey, c1: CipherTRLWE) -> CipherTRLWE {
    let (a0, b0) = c0.describe();
    let b_floor = (b0 >> (31 - NBIT)) as usize;
    let offset = 2u32.pow(30 - NBIT as u32);

    let mut c_ret = rotate_trlwe_cipher(c1, 2 * N - b_floor);
    for j in 0..tlwe::N {
        let a_floor = ((a0[j] + offset) >> (31 - NBIT)) as usize;
        let c_ret_rot = rotate_trlwe_cipher(c_ret, a_floor);
        c_ret = cmux(bk.get(j), c_ret_rot, c_ret);
    }
    c_ret
}

pub fn gate_bootstrapping(c0: CipherTLWELv0, sk: SecretKey) -> CipherTLWELv1 {
    let tv = TRLWE::new(sk).test_vector();
    let bk = TRGSW::new(sk).bootstrapping_key();
    let c = blind_rotate(c0, bk, tv);
    sample_extract_index(c, 0)
}

#[test]
fn test_decomposition() {
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    // (1) random BRing, (2) encrypt by TRLWE
    let sk = SecretKey::new();
    let bs: [bool; N] = random_bool_initialization();
    let trlwe = TRLWE::new(sk);
    let (a, b) = trlwe.encrypt(bs).describe();

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
    assert_eq!(bs, trlwe.decrypt(CipherTRLWE(a_, b_)))
}

#[test]
fn test_zero_matrix_add() {
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    let sk = SecretKey::new();
    let trlwe = TRLWE::new(sk);
    let trgsw = TRGSW::new(sk);

    let zm = trgsw.zero_matrix();

    for i in 0..zm.len() {
        let bs: [bool; N] = random_bool_initialization();
        let c = trlwe.encrypt(bs);
        let z = CipherTRLWE(zm[i][0], zm[i][1]);
        let dec_bs = trlwe.decrypt(c + z);
        assert_eq!(bs, dec_bs);
    }
}

#[test]
fn test_zero_matrix_multiple() {
    use super::ops::vsub;
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    let sk = SecretKey::new();

    let bs: [bool; N] = random_bool_initialization();
    let offset = [2u32.pow(28); N];

    // Encrypt as TRLWE
    let trlwe = TRLWE::new(sk);
    let c = trlwe.encrypt(bs);

    // Get matrix
    let trgsw = TRGSW::new(sk);
    let matrix = trgsw.zero_matrix();

    // Calc external product
    let c_ = external_product(matrix, c);
    let ring = trlwe.decrypt_torus(c_);

    let m = vsub(ring, offset);
    let mut counter = 0;
    for i in 0..N {
        counter += (m[i] <= 2u32.pow(31)) as usize;
    }
    assert!(counter == 0, "counter is {}", counter);
}

#[test]
fn test_external_product() {
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    let sk = SecretKey::new();

    let bs: [bool; N] = random_bool_initialization();

    let trlwe = TRLWE::new(sk);
    let c = trlwe.encrypt(bs);

    let trgsw = TRGSW::new(sk);
    let matrix = trgsw.coefficient(1);

    let c_ = external_product(matrix, c);
    let dec_bs = trlwe.decrypt(c_);

    let mut counter = 0;
    for i in 0..N {
        counter += (bs[i] != dec_bs[i]) as usize;
    }

    assert!(counter == 0, "counter is {}", counter);
}

#[test]
fn test_cmux() {
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    let sk = SecretKey::new();
    let trlwe = TRLWE::new(sk);

    let bs1: [bool; N] = random_bool_initialization();
    let bs0: [bool; N] = random_bool_initialization();
    let c0 = trlwe.encrypt(bs0);
    let c1 = trlwe.encrypt(bs1);

    let trgsw = TRGSW::new(sk);

    let c = trgsw.cmux(true, c0, c1);
    let dec_bs = trlwe.decrypt(c);
    let mut counter1 = 0;
    for i in 0..N {
        counter1 += (bs0[i] != dec_bs[i]) as usize;
    }

    let c = trgsw.cmux(false, c0, c1);
    let dec_bs = trlwe.decrypt(c);
    let mut counter2 = 0;
    for i in 0..N {
        counter2 += (bs1[i] != dec_bs[i]) as usize;
    }
    assert!(
        (counter1 == 0) && (counter2 == 0),
        "counter is ... true: {}, false: {}",
        counter1,
        counter2
    );
}

// FFTじゃないと重すぎて時間がかかるのでやめとく（通ることは確認済み）

// #[test]
// fn test_gate_bootstrapping() {
//     use super::ops::dot;
//     use super::tlwe::TLWE;

//     fn decrypt_as_tlwe_lv1(ext_a: RingLv1, ext_b: Torus, s: RingLv1) -> bool {
//         let m = ext_b.wrapping_sub(dot(ext_a, s)).wrapping_sub(2u32.pow(28));
//         m < 2u32.pow(31)
//     }

//     const LOOP: usize = 2;
//     let bs: [bool; LOOP] = [true, false];

//     for i in 0..LOOP {
//         let pt = bs[i];
//         let sk = SecretKey::new();
//         let tlwe = TLWE::new(sk);
//         let c = tlwe.encrypt(pt);
//         let (a_, b_) = gate_bootstrapping(c, sk).describe();
//         let msg = decrypt_as_tlwe_lv1(a_, b_, sk.lv1);

//         assert_eq!(pt, !msg);
//     }
// }
