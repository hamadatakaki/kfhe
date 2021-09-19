use super::key::SecretKey;
use super::ops::{pmul, vadd, vsub};
use super::params::trlwe;
use super::sampling::{ndim_modular_normal_dist, ndim_torus_uniform};
use super::tlwe::CipherTLWELv1;
use super::util::{boolpoly_normalization, fring_to_torus_ring, RingLv1};

const N: usize = trlwe::N;

type BRing = [bool; N];

#[derive(Clone, Copy, Debug)]
pub struct CipherTRLWE(pub RingLv1, pub RingLv1);

impl CipherTRLWE {
    pub fn describe(self) -> (RingLv1, RingLv1) {
        (self.0, self.1)
    }
}

impl std::ops::Add for CipherTRLWE {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let (a0, b0) = self.describe();
        let (a1, b1) = rhs.describe();
        let a = vadd(a0, a1);
        let b = vadd(b0, b1);
        CipherTRLWE(a, b)
    }
}

impl std::ops::Sub for CipherTRLWE {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let (a0, b0) = self.describe();
        let (a1, b1) = rhs.describe();
        let a = vsub(a0, a1);
        let b = vsub(b0, b1);
        CipherTRLWE(a, b)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TRLWE {
    s: RingLv1,
}

impl TRLWE {
    pub fn new(sk: SecretKey) -> Self {
        Self { s: sk.lv1 }
    }

    pub fn get_secret(&self) -> RingLv1 {
        self.s
    }

    pub fn encrypt_torus(&self, msg: RingLv1) -> CipherTRLWE {
        let a = ndim_torus_uniform();
        let e = ndim_modular_normal_dist(0., trlwe::ALPHA);
        let b = vadd(vadd(pmul(a, self.s), msg), e);
        CipherTRLWE(a, b)
    }

    pub fn encrypt(&self, msg: BRing) -> CipherTRLWE {
        let m = fring_to_torus_ring(boolpoly_normalization(msg));
        self.encrypt_torus(m)
    }

    pub fn decrypt_torus(&self, c: CipherTRLWE) -> RingLv1 {
        let (a, b) = c.describe();
        vsub(b, pmul(a, self.get_secret()))
    }

    pub fn decrypt(&self, c: CipherTRLWE) -> BRing {
        let m = self.decrypt_torus(c);
        let mut bs = [false; N];
        for i in 0..N {
            bs[i] = m[i] <= 2u32.pow(31);
        }
        bs
    }

    pub fn test_vector(&self) -> CipherTRLWE {
        self.encrypt([true; trlwe::N])
    }
}

pub fn sample_extract_index(c: CipherTRLWE, k: usize) -> CipherTLWELv1 {
    if k > N - 1 {
        panic!("ArrayIndexOutOfBoundsException")
    }

    let (a, b) = c.describe();
    let mut ext_a = [0; N];
    for i in 0..N {
        if i <= k {
            ext_a[i] = a[k - i];
        } else {
            ext_a[i] = 0u32.wrapping_sub(a[N + k - i]);
        }
    }
    CipherTLWELv1(ext_a, b[k])
}

#[test]
fn test_trlwe_enc_and_dec() {
    use super::sampling::random_bool_initialization;

    fn _run_trlwe(bs: BRing) -> BRing {
        let sk = SecretKey::new();
        let trlwe = TRLWE::new(sk);
        let c = trlwe.encrypt(bs);
        trlwe.decrypt(c)
    }

    let bs: BRing = random_bool_initialization();
    assert_eq!(bs, _run_trlwe(bs));
}

#[test]
fn test_sample_extract_index() {
    use super::ops::dot;
    use super::sampling::random_bool_initialization;
    use super::util::Torus;

    use rand;
    use rand_distr::{Distribution, Uniform};

    fn decrypt_as_tlwe_lv1(ext_a: RingLv1, ext_b: Torus, s: RingLv1) -> bool {
        let m = ext_b.wrapping_sub(dot(ext_a, s)).wrapping_sub(2u32.pow(28));
        m < 2u32.pow(31)
    }

    let sk = SecretKey::new();

    const LOOP: usize = 64;

    let uni = Uniform::new_inclusive(0, N - 1);
    let mut rng = rand::thread_rng();
    let mut indecies = [0; LOOP];
    indecies[0] = 0;
    for i in 1..LOOP {
        indecies[i] = uni.sample(&mut rng);
    }

    let mut counter = 0;
    loop {
        if counter >= LOOP {
            break;
        }

        let index = indecies[counter];
        let bs: BRing = random_bool_initialization();

        // Encrypt as TRLWE
        let trlwe = TRLWE::new(sk);
        let c = trlwe.encrypt(bs);

        // Sample Extract Index
        let (ext_a, ext_b) = sample_extract_index(c, index).describe();

        let msg = decrypt_as_tlwe_lv1(ext_a, ext_b, trlwe.s);
        assert_eq!(bs[index], msg);

        counter += 1;
    }
}
