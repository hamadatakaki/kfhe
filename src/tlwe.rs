use super::params::tlwe;
use super::sampling::{modular_normal_dist, ndim_bin_uniform, ndim_torus_uniform};
use super::util::ops::dot;
use super::util::{bool_normalization, float_to_torus, Torus};

const N: usize = tlwe::N;

#[derive(Clone, Debug)]
pub struct TLWE {
    s: [Torus; N],
}

impl TLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform();
        Self { s }
    }

    pub fn encrypt(&self, msg: bool) -> ([Torus; N], Torus) {
        let m = float_to_torus(bool_normalization(msg));
        let s = self.s.clone();
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(&a, &s).wrapping_add(m).wrapping_add(e);
        (a, b)
    }

    pub fn decrypt(&self, a: [Torus; N], b: Torus) -> bool {
        let s = self.s.clone();
        let m = b.wrapping_sub(dot(&a, &s)).wrapping_sub(2u32.pow(28));
        m < 2u32.pow(31)
    }
}

#[test]
fn test_tlwe_enc_and_dec() {
    use super::sampling::random_bool_initialization;

    fn _run_tlwe(msg: bool) -> bool {
        let tlwe = TLWE::new();
        let (a, b) = tlwe.encrypt(msg);
        tlwe.decrypt(a, b)
    }

    const T: usize = 1000;
    let bs: [bool; T] = random_bool_initialization();
    for i in 0..T {
        let b = bs[i];
        assert_eq!(b, _run_tlwe(b), "{}", i + 1);
    }
}
