use super::params::N;
use super::sampling::{modular_normal_dist, ndim_bin_uniform, ndim_torus_uniform};
use super::util::torus_to_ring;
use super::util::vec::dot;

#[derive(Clone, Debug)]
pub struct TLWE {
    s: Vec<u32>,
    alpha: f64,
}

impl TLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform(N);
        let alpha = 2f64.powi(-15);
        Self { s, alpha }
    }

    pub fn encrypt(&self, msg: u32) -> (Vec<u32>, u32) {
        let m = torus_to_ring((2. * (msg as f64) - 1.) / 8.);
        let s = self.s.clone();
        let a = ndim_torus_uniform(N);
        let e = modular_normal_dist(0., self.alpha);
        let b = dot(&a, &s).wrapping_add(m).wrapping_add(e);
        return (a, b);
    }

    pub fn decrypt(&self, a: Vec<u32>, b: u32) -> u32 {
        let s = self.s.clone();
        let m = b.wrapping_sub(dot(&a, &s)).wrapping_add(2u32.pow(29)) as i64;
        let m = m % 2i64.pow(32);
        let m = 2i64.pow(31) - m;
        ((1 + m.signum()) / 2) as u32
    }
}
