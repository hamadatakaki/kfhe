use super::params::trlwe;
use super::sampling::{ndim_bin_uniform, ndim_modular_normal_dist, ndim_torus_uniform};
use super::util::ops::{pmul, vadd, vsub};
use super::util::{boolpoly_normalization, fring_to_torus_ring, BRing, Ring, Torus};

pub struct TRLWE {
    s: Ring,
}

impl TRLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform();
        Self { s }
    }

    pub fn encrypt(&self, msg: BRing) -> (Ring, Ring) {
        let m = fring_to_torus_ring(boolpoly_normalization(msg));
        let s = self.s.clone();
        let a: Ring = ndim_torus_uniform();
        let e: Ring = ndim_modular_normal_dist(0., trlwe::ALPHA);
        let b = vadd(&vadd(&pmul(&a, &s), &m), &e);
        (a, b)
    }

    pub fn decrypt(&self, a: Ring, b: Ring) -> BRing {
        let mut bs = [false; trlwe::N];
        let offset = [2u32.pow(29); trlwe::N];
        let s = self.s.clone();
        let m = vsub(&vsub(&b, &pmul(&a, &s)), &offset);
        for i in 0..trlwe::N {
            bs[i] = m[i] > 2u32.pow(31);
        }
        bs
    }
}

pub fn sample_extract_index((a, b): (Ring, Ring), k: usize) -> (Ring, Torus) {
    let n = trlwe::N;
    if k > n - 1 {
        panic!("ArrayIndexOutOfBoundsException")
    }
    let mut ext_a = [0; trlwe::N];
    for i in 0..n {
        if i <= k {
            ext_a[i] = a[k - i];
        } else {
            ext_a[i] = 0u32.wrapping_sub(a[n + k - i]);
        }
    }
    (ext_a, b[k])
}
