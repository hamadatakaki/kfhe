use super::params::tlwe;
use super::sampling::{modular_normal_dist, ndim_bin_uniform, ndim_torus_uniform};
use super::util::vec::dot;
use super::util::{bool_normalization, float_to_torus, Torus};

#[derive(Clone, Debug)]
pub struct TLWE {
    s: [Torus; tlwe::N],
}

impl TLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform();
        Self { s }
    }

    pub fn encrypt(&self, msg: bool) -> ([Torus; tlwe::N], Torus) {
        let m = float_to_torus(bool_normalization(msg));
        let s = self.s.clone();
        let a = ndim_torus_uniform();
        let e = modular_normal_dist(0., tlwe::ALPHA);
        let b = dot(&a, &s).wrapping_add(m).wrapping_add(e);
        return (a, b);
    }

    pub fn decrypt(&self, a: [Torus; tlwe::N], b: Torus) -> bool {
        let s = self.s.clone();
        let m = b.wrapping_sub(dot(&a, &s)).wrapping_sub(2u32.pow(29)) as i64;
        let m = m % 2i64.pow(32);
        let m = 2i64.pow(31) - m;
        m.signum() < 0
    }
}
