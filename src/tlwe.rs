use super::util::ops::dot;
use super::util::params::tlwe;
use super::util::sampling::{modular_normal_dist, ndim_bin_uniform, ndim_torus_uniform};
use super::util::Torus;
use super::util::{bool_normalization, float_to_torus};

const M: usize = tlwe::N;

#[derive(Clone, Debug)]
pub struct TLWE {
    s: [Torus; M],
}

impl TLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform();
        Self { s }
    }

    pub fn encrypt_torus(&self, torus: Torus) -> ([Torus; M], Torus) {
        let s = self.s.clone();
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(&a, &s).wrapping_add(torus).wrapping_add(e);
        (a, b)
    }

    pub fn encrypt(&self, msg: bool) -> ([Torus; M], Torus) {
        let m = float_to_torus(bool_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, a: [Torus; M], b: Torus) -> Torus {
        let s = self.s.clone();
        b.wrapping_sub(dot(&a, &s))
    }

    pub fn decrypt(&self, a: [Torus; M], b: Torus) -> bool {
        let m = self.decrypt_torus(a, b);
        m < 2u32.pow(31)
    }
}

#[test]
fn test_tlwe_enc_and_dec() {
    use super::util::sampling::random_bool_initialization;

    fn _run_tlwe(msg: bool) -> bool {
        let tlwe = TLWE::new();
        let (a, b) = tlwe.encrypt(msg);
        let m = tlwe.decrypt_torus(a, b);
        m.wrapping_sub(2u32.pow(28)) < 2u32.pow(31)
    }

    const T: usize = 1000;
    let bs: [bool; T] = random_bool_initialization();
    for i in 0..T {
        let b = bs[i];
        assert_eq!(b, _run_tlwe(b), "{}", i + 1);
    }
}
