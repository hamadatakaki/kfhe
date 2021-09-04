use super::util::ops::{pmul, vadd, vsub};
use super::util::params::trlwe;
use super::util::sampling::{ndim_bin_uniform, ndim_modular_normal_dist, ndim_torus_uniform};
use super::util::{boolpoly_normalization, fring_to_torus_ring, RingLv1, Torus};

const N: usize = trlwe::N;

type BRing = [bool; N];

pub struct TRLWE {
    s: RingLv1,
}

impl TRLWE {
    pub fn new() -> Self {
        let s = ndim_bin_uniform();
        Self::from_secret(s)
    }

    pub fn from_secret(s: RingLv1) -> Self {
        Self { s }
    }

    pub fn get_secret(&self) -> RingLv1 {
        self.s.clone()
    }

    pub fn encrypt_torus(&self, msg: RingLv1) -> (RingLv1, RingLv1) {
        let s = self.s.clone();
        let a: RingLv1 = ndim_torus_uniform();
        let e: RingLv1 = ndim_modular_normal_dist(0., trlwe::ALPHA);
        let b = vadd(&vadd(&pmul(&a, &s), &msg), &e);
        (a, b)
    }

    pub fn encrypt(&self, msg: BRing) -> (RingLv1, RingLv1) {
        let m = fring_to_torus_ring(boolpoly_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, a: RingLv1, b: RingLv1) -> RingLv1 {
        vsub(&b, &pmul(&a, &self.get_secret()))
    }

    pub fn decrypt(&self, a: RingLv1, b: RingLv1) -> BRing {
        let m = self.decrypt_torus(a, b);
        let mut bs = [false; N];
        for i in 0..N {
            bs[i] = m[i] <= 2u32.pow(31);
        }
        bs
    }
}

pub fn sample_extract_index(a: RingLv1, b: RingLv1, k: usize) -> (RingLv1, Torus) {
    if k > N - 1 {
        panic!("ArrayIndexOutOfBoundsException")
    }
    let mut ext_a = [0; N];
    for i in 0..N {
        if i <= k {
            ext_a[i] = a[k - i];
        } else {
            ext_a[i] = 0u32.wrapping_sub(a[N + k - i]);
        }
    }
    (ext_a, b[k])
}

#[test]
fn test_trlwe_enc_and_dec() {
    use super::util::sampling::random_bool_initialization;

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
    use super::util::ops::dot;
    use super::util::sampling::random_bool_initialization;

    use rand;
    use rand_distr::{Distribution, Uniform};

    fn decrypt_as_tlwe(ext_a: RingLv1, ext_b: Torus, s: RingLv1) -> bool {
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
    let (ext_a, ext_b) = sample_extract_index(a, b, index);

    let msg = decrypt_as_tlwe(ext_a, ext_b, trlwe.s);
    assert_eq!(bs[index], msg);
}
