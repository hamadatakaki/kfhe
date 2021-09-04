use super::key::SecretKey;
use super::ops::{dot, vadd};
use super::params::tlwe;
use super::sampling::{modular_normal_dist, ndim_torus_uniform};
use super::util::{bool_normalization, float_to_torus, RingLv0, Torus};

#[derive(Clone, Copy, Debug)]
pub struct CipherTLWELv0(RingLv0, Torus);

impl CipherTLWELv0 {
    pub fn describe(self) -> (RingLv0, Torus) {
        (self.0, self.1)
    }
}

impl std::ops::Add for CipherTLWELv0 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let (a0, b0) = self.describe();
        let (a1, b1) = rhs.describe();
        let a = vadd(a0, a1);
        let b = b0.wrapping_add(b1);
        CipherTLWELv0(a, b)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TLWE {
    s: RingLv0,
}

impl TLWE {
    pub fn new(sk: SecretKey) -> Self {
        Self { s: sk.lv0 }
    }

    pub fn encrypt_torus(&self, torus: Torus) -> CipherTLWELv0 {
        let s = self.s.clone();
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(a, s).wrapping_add(torus).wrapping_add(e);
        CipherTLWELv0(a, b)
    }

    pub fn encrypt(&self, msg: bool) -> CipherTLWELv0 {
        let m = float_to_torus(bool_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, c: CipherTLWELv0) -> Torus {
        let (a, b) = c.describe();
        let s = self.s.clone();
        b.wrapping_sub(dot(a, s))
    }

    pub fn decrypt(&self, c: CipherTLWELv0) -> bool {
        let m = self.decrypt_torus(c);
        m < 2u32.pow(31)
    }
}

#[test]
fn test_tlwe_enc_and_dec() {
    use super::sampling::random_bool_initialization;

    fn _run_tlwe(sk: SecretKey, msg: bool) -> bool {
        let tlwe = TLWE::new(sk);
        let c = tlwe.encrypt(msg);
        let m = tlwe.decrypt_torus(c);
        m.wrapping_sub(2u32.pow(28)) < 2u32.pow(31)
    }

    const T: usize = 1000;
    let bs: [bool; T] = random_bool_initialization();
    let sk = SecretKey::new();
    for i in 0..T {
        let b = bs[i];
        assert_eq!(b, _run_tlwe(sk, b), "{}", i + 1);
    }
}
