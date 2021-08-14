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
    use super::util::params::trlwe::BRing;
    use super::util::sampling::random_bool_initialization;

    // (1) random BRing, (2) encrypt by TRLWE
    let b_poly: BRing = random_bool_initialization();
    let trlwe = TRLWE::new();
    let (a, b) = trlwe.encrypt(b_poly);

    // (3) decomposition a and b
    let a_hat = decomposition(a);
    let b_hat = decomposition(b);

    // (4) reconstruct a and b by decomposition
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
    for j in 0..N {
        let mut sa: i64 = 0;
        let mut sb: i64 = 0;
        for i in 0..L {
            sa += (a_hat[i][j] as i64) * (BG as i64).pow((L - i) as u32 - 1);
            sb += (b_hat[i][j] as i64) * (BG as i64).pow((L - i) as u32 - 1);
        }
        sa = sa.rem_euclid(2i64.pow(L as u32 * trgsw::BGBIT)) << (32 - L as u32 * trgsw::BGBIT);
        sb = sb.rem_euclid(2i64.pow(L as u32 * trgsw::BGBIT)) << (32 - L as u32 * trgsw::BGBIT);
        a_[j] = sa as u32;
        b_[j] = sb as u32;
    }

    // (5) decrypt by TRLWE, and assertion!
    assert_eq!(b_poly, trlwe.decrypt(a_, b_))
}
