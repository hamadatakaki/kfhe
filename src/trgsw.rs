use super::params::trgsw;
use super::util::ops::{pmul, vadd, vsub};
use super::util::{float_to_torus, BRing, Ring, Torus};

type Zn_X = [i8; trgsw::N];

pub fn decomposition(poly: Ring) -> [Zn_X; trgsw::L] {
    unimplemented!()
}

#[test]
fn test_decomposition() {
    use super::sampling::random_bool_initialization;
    use super::trlwe::TRLWE;

    fn unsigning(s: i8) -> u32 {
        if s < 0 {
            (s + 32) as u32 + 32
        } else {
            s as u32
        }
    }

    // (1) random BRing, (2) encrypt by TRLWE
    let b_poly: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(b_poly);

    // (3) decomposition a and b
    let a_hat = decomposition(a);
    let b_hat = decomposition(b);

    // (4) reconstruct a and b by decomposition
    let mut weight: [Torus; trgsw::L] = [0; trgsw::L];
    for i in 0..trgsw::L {
        weight[i] = float_to_torus((trgsw::BG as f64).powi(-(i as i32 + 1)))
    }
    let mut a_: Ring = [0; trgsw::N];
    let mut b_: Ring = [0; trgsw::N];
    for j in 0..trgsw::N {
        let mut sa: Torus = 0;
        let mut sb: Torus = 0;
        for i in 0..trgsw::L {
            sa = sa.wrapping_add(unsigning(a_hat[j][i]).wrapping_mul(weight[i]));
            sb = sb.wrapping_add(unsigning(b_hat[j][i]).wrapping_mul(weight[i]));
        }
        a_[j] = sa;
        b_[j] = sb;
    }

    // (5) decrypt by TRLWE, and assertion!
    assert_eq!(b_poly, trlwe.decrypt(a_, b_))
}
