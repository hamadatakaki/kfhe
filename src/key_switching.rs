use super::key::SecretKey;
use super::params::{tlwe, trgsw};
use super::tlwe::{CipherTLWELv0, CipherTLWELv1, TLWE};
use super::util::{float_to_torus, RingLv0, Torus};

const N: usize = trgsw::N;
const T: usize = trgsw::T;
const BASEBIT: u32 = trgsw::BASEBIT;

const K: usize = 1 << BASEBIT;

#[derive(Clone, Debug)]
pub struct KeySwitchingKey(pub Vec<CipherTLWELv0>);

impl KeySwitchingKey {
    fn new(sk: SecretKey) -> Self {
        let s = sk.lv1;
        let tlwe = TLWE::new(sk);

        let mut v = vec![CipherTLWELv0::empty(); (K - 1) * T * N];

        for k in 1..K {
            for j in 0..T {
                for i in 0..N {
                    let base = 1 << (j as Torus + 1) * BASEBIT;
                    let msg = (k as f64 * s[i] as f64) / base as f64;
                    let msg = float_to_torus(msg - 0.5);
                    v[i + j * N + (k - 1) * N * T] = tlwe.encrypt_torus(msg);
                }
            }
        }

        assert_eq!(v.len(), (K - 1) * T * N);

        Self(v)
    }

    fn access(&self, i: usize, j: usize, k: usize) -> CipherTLWELv0 {
        assert!(i < N);
        assert!(j < T);
        assert!(0 < k && k < K);

        self.0[i + j * N + (k - 1) * N * T]
    }
}

pub fn identity_key_switching(c: CipherTLWELv1, sk: SecretKey) -> CipherTLWELv0 {
    let (a, b) = c.describe();
    let ks = KeySwitchingKey::new(sk);

    let a0: RingLv0 = [0; tlwe::N];
    let mut c0 = CipherTLWELv0(a0, b);

    let offset: Torus = 1 << (31 - T * BASEBIT as usize);

    for i in 0..N {
        let ai_ = a[i].wrapping_add(offset);
        for j in 0..T {
            let shift = 32 - (j + 1) * BASEBIT as usize;
            let k = (ai_ >> shift) as usize % K;
            if k != 0 {
                c0 = c0 - ks.access(i, j, k);
            }
        }
    }

    c0
}

#[test]
fn test_identity_key_switching() {
    use super::sampling::random_bool_initialization;
    use super::tlwe::{TLWELv1, TLWE};

    let bs: [bool; 16] = random_bool_initialization();
    let mut count = 0;

    for b in bs {
        let sk = SecretKey::new();
        let tlwe0 = TLWE::new(sk);
        let tlwe1 = TLWELv1::new(sk);
        let c1 = tlwe1.encrypt(b);
        let c0 = identity_key_switching(c1, sk);
        let msg = tlwe0.decrypt(c0);

        count += (b != msg) as usize;
    }

    assert!(count == 0, "count: {}", count);
}
