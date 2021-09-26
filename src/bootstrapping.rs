use super::key::SecretKey;
use super::params::{tlwe, trgsw};
use super::tlwe::{CipherTLWELv0, CipherTLWELv1};
use super::trgsw::{cmux, TRGSWMatrix, TRGSW};
use super::trlwe::{sample_extract_index, CipherTRLWE, TRLWE};
use super::util::rotate_ring;

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const NBIT: usize = trgsw::NBIT;

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

pub fn bootstrapping_key(sk: SecretKey) -> BootstrappingKey {
    let trgsw = TRGSW::new(sk);
    let lv0 = sk.lv0;
    let mut bk = BootstrappingKey(vec![[[[0; N]; 2]; 2 * L]; tlwe::N]);

    for (j, &item) in lv0.iter().enumerate().take(tlwe::N) {
        bk.set(j, trgsw.coefficient(item as i8));
    }
    bk
}

fn rotate_trlwe_cipher(c: CipherTRLWE, k: usize) -> CipherTRLWE {
    let (a, b) = c.describe();
    let a_rot = rotate_ring(a, k);
    let b_rot = rotate_ring(b, k);
    CipherTRLWE(a_rot, b_rot)
}

pub fn blind_rotate(c0: CipherTLWELv0, c1: CipherTRLWE, bk: BootstrappingKey) -> CipherTRLWE {
    let (a0, b0) = c0.describe();
    let b_floor = (b0 >> (31 - NBIT)) as usize;
    let offset = 2u32.pow(30 - NBIT as u32);

    let mut c_ret = rotate_trlwe_cipher(c1, 2 * N - b_floor);
    for (j, a0) in a0.iter().enumerate().take(tlwe::N) {
        let a_floor = ((a0 + offset) >> (31 - NBIT)) as usize;
        let c_ret_rot = rotate_trlwe_cipher(c_ret, a_floor);
        c_ret = cmux(bk.get(j), c_ret_rot, c_ret);
    }
    c_ret
}

pub fn gate_bootstrapping(c0: CipherTLWELv0, sk: SecretKey) -> CipherTLWELv1 {
    let tv = TRLWE::new(sk).test_vector();
    let bk = bootstrapping_key(sk);
    let c = blind_rotate(c0, tv, bk);
    sample_extract_index(c, 0)
}

#[test]
fn test_trlwe_rotate() {
    use rand_distr::{Distribution, Uniform};

    use super::tlwe::TLWELv1;

    let sk = SecretKey::new();
    let trlwe = TRLWE::new(sk);
    let tlwe1 = TLWELv1::new(sk);
    let tv = trlwe.test_vector();

    // k is random in [0, 2*N)
    let uni = Uniform::new_inclusive(0, 2 * N - 1);
    let mut rng = rand::thread_rng();
    let k = uni.sample(&mut rng);

    let rot = rotate_trlwe_cipher(tv, k);
    for i in 0..N {
        let q = (k + i) / N;
        let r = (k + i) % N;
        let mut rot_c = sample_extract_index(rot, r);
        for _ in 0..q {
            rot_c = -rot_c;
        }
        let rot_b = tlwe1.decrypt(rot_c);
        let c = sample_extract_index(tv, i);
        let b = tlwe1.decrypt(c);
        assert_eq!(b, rot_b);
    }
}

#[test]
fn test_blind_rotate() {
    use super::tlwe::TLWELv1;

    let sk = SecretKey::new();
    let tlwe1 = TLWELv1::new(sk);
    let c0 = CipherTLWELv0::clearly_true();
    let trlwe = TRLWE::new(sk);
    let tv = trlwe.test_vector();
    let bk = bootstrapping_key(sk);

    let rot_c = blind_rotate(c0, tv, bk);
    let c1 = sample_extract_index(rot_c, 0);
    let b1 = tlwe1.decrypt(c1);

    assert_eq!(true, b1);
}

#[test]
fn test_blind_rotate_2() {
    use super::tlwe::TLWELv1;

    let sk = SecretKey::new();
    let tlwe1 = TLWELv1::new(sk);
    let c0 = CipherTLWELv0::clearly_true();
    let trlwe = TRLWE::new(sk);
    let tv = trlwe.test_vector();
    let bk = bootstrapping_key(sk);

    let rot_c = blind_rotate(-c0, tv, bk);
    let c1 = sample_extract_index(rot_c, 0);
    let b1 = tlwe1.decrypt(c1);

    assert_eq!(false, b1);
}

// FFTじゃないと重すぎて時間がかかる

// #[test]
// fn test_gate_bootstrapping() {
//     use super::sampling::random_bool_initialization;
//     use super::tlwe::{TLWELv1, TLWE};

//     let bs: [bool; 16] = random_bool_initialization();
//     let mut count = 0;

//     for b in bs {
//         let sk = SecretKey::new();
//         let tlwe0 = TLWE::new(sk);
//         let tlwe1 = TLWELv1::new(sk);
//         let c = tlwe0.encrypt(b);
//         let c1 = gate_bootstrapping(c, sk);
//         let msg = tlwe1.decrypt(c1);

//         count += (b ^ msg) as usize;
//     }

//     assert!(count == 0, "count: {}", count);
// }
