use super::util::params::trgsw;
use super::util::params::trlwe::Ring;

const N: usize = trgsw::N;
const L: usize = trgsw::L;
const BG: u32 = trgsw::BG;

type Z = trgsw::Z;
type Zpoly = [Z; N];

pub fn decomposition(poly: Ring) -> [Zpoly; L] {
    let mut decomped: [Zpoly; L] = [[0; N]; L];
    for n in 0..N {
        let mut a = poly[n] / 2u32.pow(32 - L as u32 * trgsw::BGBIT);
        let mut cflag = false;

        for l in 0..L {
            let mut r = a % BG;
            r += cflag as u32;

            // type-cast が乱暴なため注意. BGBIT, Zの型が違うと動かなくなる.
            let s = if r >= (BG / 2) {
                cflag = true;
                (r - (BG / 2)) as Z - (BG / 2) as Z
            } else {
                cflag = false;
                r as Z
            };

            if s < trgsw::SIGN_MIN {
                println!("{} -> {}", r, s);
                assert!(s >= trgsw::SIGN_MIN);
            }
            if s > trgsw::SIGN_MAX {
                println!("{} -> {}", r, s);
                assert!(s <= trgsw::SIGN_MAX);
            }

            decomped[L - l - 1][n] = s;
            a /= BG;
        }
    }
    decomped
}

#[test]
fn test_decomposition() {
    // fail case

    use super::trlwe::TRLWE;
    use super::util::float_to_torus;
    use super::util::params::trlwe::{BRing, Torus};
    use super::util::sampling::random_bool_initialization;

    fn unsigning(s: Z) -> u32 {
        if s < 0 {
            (s + (trgsw::BG / 2) as i8) as u32 + (trgsw::BG / 2)
        } else {
            s as u32
        }
    }

    // // Debug
    // fn debug_array_u32_to_f64(ring: Ring) -> [f64; N] {
    //     let mut arr = [0.; N];
    //     for i in 0..N {
    //         let r = ring[i] as f64;
    //         arr[i] = r / 2f64.powi(32) - 0.5;
    //     }
    //     arr
    // }

    // (1) random BRing, (2) encrypt by TRLWE
    let b_poly: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(b_poly);

    // (3) decomposition a and b
    let a_hat = decomposition(a);
    let b_hat = decomposition(b);

    // (4) reconstruct a and b by decomposition
    let mut weight: [Torus; L] = [0; L];
    for i in 0..L {
        weight[i] = float_to_torus((trgsw::BG as f64).powi(-(i as i32) - 1))
    }
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
    for j in 0..N {
        let mut sa: Torus = 0;
        let mut sb: Torus = 0;
        for i in 0..L {
            sa = sa.wrapping_add(unsigning(a_hat[i][j]).wrapping_mul(weight[i]));
            sb = sb.wrapping_add(unsigning(b_hat[i][j]).wrapping_mul(weight[i]));
        }
        a_[j] = sa;
        b_[j] = sb;
    }

    // // Debug
    // println!(" a: {:?}", debug_array_u32_to_f64(a));
    // println!("Da: {:?}", debug_array_u32_to_f64(a_));
    // println!(" b: {:?}", debug_array_u32_to_f64(b));
    // println!("Db: {:?}", debug_array_u32_to_f64(b_));

    // (5) decrypt by TRLWE, and assertion!
    assert_eq!(b_poly, trlwe.decrypt(a_, b_))
}
