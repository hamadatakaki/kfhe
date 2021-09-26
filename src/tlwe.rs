use super::key::SecretKey;
use super::ops::{dot, vadd, vsub};
use super::params::tlwe;
use super::sampling::{modular_normal_dist, ndim_torus_uniform};
use super::util::{
    bool_normalization, float_to_torus, ring_negative, torus_negative, RingLv0, RingLv1, Torus,
};

#[derive(Clone, Copy, Debug)]
pub struct CipherTLWELv0(pub RingLv0, pub Torus);

impl CipherTLWELv0 {
    pub fn describe(self) -> (RingLv0, Torus) {
        (self.0, self.1)
    }

    pub fn empty() -> Self {
        let a = [0; tlwe::N];
        let b = 0;
        Self(a, b)
    }

    pub fn clearly_true() -> Self {
        let a = [0; tlwe::N];
        let b = float_to_torus(0.125);
        Self(a, b)
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

impl std::ops::Sub for CipherTLWELv0 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let (a0, b0) = self.describe();
        let (a1, b1) = rhs.describe();
        let a = vsub(a0, a1);
        let b = b0.wrapping_sub(b1);
        CipherTLWELv0(a, b)
    }
}

impl std::ops::Neg for CipherTLWELv0 {
    type Output = Self;
    fn neg(self) -> Self {
        let (a0, b0) = self.describe();
        let b = torus_negative(b0);
        let a = ring_negative(a0);
        CipherTLWELv0(a, b)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CipherTLWELv1(pub RingLv1, pub Torus);

impl CipherTLWELv1 {
    pub fn describe(self) -> (RingLv1, Torus) {
        (self.0, self.1)
    }
}

impl std::ops::Neg for CipherTLWELv1 {
    type Output = Self;
    fn neg(self) -> Self {
        let (a0, b0) = self.describe();
        let b = torus_negative(b0);
        let a = ring_negative(a0);
        CipherTLWELv1(a, b)
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
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(a, self.s).wrapping_add(torus).wrapping_add(e);
        CipherTLWELv0(a, b)
    }

    pub fn encrypt(&self, msg: bool) -> CipherTLWELv0 {
        let m = float_to_torus(bool_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, c: CipherTLWELv0) -> Torus {
        let (a, b) = c.describe();
        b.wrapping_sub(dot(a, self.s))
    }

    pub fn decrypt(&self, c: CipherTLWELv0) -> bool {
        let m = self.decrypt_torus(c);
        m < 2u32.pow(31)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TLWELv1 {
    s: RingLv1,
}

impl TLWELv1 {
    pub fn new(sk: SecretKey) -> Self {
        Self { s: sk.lv1 }
    }

    pub fn encrypt_torus(&self, torus: Torus) -> CipherTLWELv1 {
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(a, self.s).wrapping_add(torus).wrapping_add(e);
        CipherTLWELv1(a, b)
    }

    pub fn encrypt(&self, msg: bool) -> CipherTLWELv1 {
        let m = float_to_torus(bool_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, c: CipherTLWELv1) -> Torus {
        let (a, b) = c.describe();
        b.wrapping_sub(dot(a, self.s))
    }

    pub fn decrypt(&self, c: CipherTLWELv1) -> bool {
        let m = self.decrypt_torus(c);
        m < 2u32.pow(31)
    }
}

pub fn tlwe_nand(x: bool, y: bool) -> bool {
    let sk = SecretKey::new();
    let tlwe = TLWE::new(sk);
    let c0 = CipherTLWELv0::clearly_true();
    let c1 = tlwe.encrypt(x);
    let c2 = tlwe.encrypt(y);
    let c = c0 - (c1 + c2);
    tlwe.decrypt(c)
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
    for b in bs {
        assert_eq!(b, _run_tlwe(sk, b));
    }
}
