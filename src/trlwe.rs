use super::params::trlwe;
use super::sampling::{ndim_bin_uniform, ndim_modular_normal_dist, ndim_torus_uniform};
use super::util::ops::{pmul, vadd, vsub};
use super::util::{boolpoly_normalization, fring_to_torus_ring, BRing, Ring, Torus};

const N: usize = trlwe::N;

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
        let mut bs = [false; N];
        let offset = [2u32.pow(28); N];
        let s = self.s.clone();
        let m = vsub(&vsub(&b, &pmul(&a, &s)), &offset);
        for i in 0..N {
            bs[i] = m[i] <= 2u32.pow(31);
        }
        bs
    }
}

pub fn sample_extract_index((a, b): (Ring, Ring), k: usize) -> (Ring, Torus) {
    let n = N;
    if k > n - 1 {
        panic!("ArrayIndexOutOfBoundsException")
    }
    let mut ext_a = [0; N];
    for i in 0..n {
        if i <= k {
            ext_a[i] = a[k - i];
        } else {
            ext_a[i] = 0u32.wrapping_sub(a[n + k - i]);
        }
    }
    (ext_a, b[k])
}

#[test]
fn test_trlwe_enc_and_dec() {
    use super::sampling::random_bool_initialization;

    fn _run_trlwe(bs: BRing) -> BRing {
        let trlwe = TRLWE::new();
        let (a, b) = trlwe.encrypt(bs);
        trlwe.decrypt(a, b)
    }

    let bs: BRing = random_bool_initialization();
    assert_eq!(bs, _run_trlwe(bs));
}

#[test]
fn test_sample_extract_index() {
    use super::sampling::random_bool_initialization;
    use super::util::ops::dot;
    use super::util::Torus;

    use rand;
    use rand_distr::{Distribution, Uniform};

    fn decrypt_as_tlwe(ext_a: Ring, ext_b: Torus, s: Ring) -> bool {
        let m = ext_b
            .wrapping_sub(dot(&ext_a, &s))
            .wrapping_sub(2u32.pow(28));
        m < 2u32.pow(31)
    }

    let index: usize = Uniform::new(0, N).sample(&mut rand::thread_rng());
    let bs: BRing = random_bool_initialization();

    // Encrypt as TRLWE
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(bs);

    // Sample Extract Index
    let (ext_a, ext_b) = sample_extract_index((a, b), index);

    let msg = decrypt_as_tlwe(ext_a, ext_b, trlwe.s);
    assert_eq!(bs[index], msg);
}
