use super::super::util::params::trlwe::Ring;
use super::super::util::params::Torus;
use super::{Decomposition, TRGSWMatrix, L, N};

fn _naive_external_product((a_bar, b_bar): Decomposition, matrix: TRGSWMatrix) -> (Ring, Ring) {
    let mut a: Ring = [0; N];
    let mut b: Ring = [0; N];
    for j in 0..N {
        let m = _extract_const_array_from_trgsw_matrix(matrix, j);
        let (a_prime, b_prime) = _const_external_product((a_bar, b_bar), m);
        for i in 0..N {
            if i + j < N {
                a[i + j] = a[i + j].wrapping_add(a_prime[i]);
                b[i + j] = b[i + j].wrapping_add(b_prime[i]);
            } else {
                a[i + j - N] = a[i + j - N].wrapping_sub(a_prime[i]);
                b[i + j - N] = b[i + j - N].wrapping_sub(b_prime[i]);
            }
        }
    }
    (a, b)
}

fn _extract_const_array_from_trgsw_matrix(matrix: TRGSWMatrix, k: usize) -> [[Torus; 2]; 2 * L] {
    // assert
    assert!(k <= N - 1, "ArrayIndexOutOfBoundsException");

    let mut m = [[0; 2]; 2 * L];
    for l in 0..2 * L {
        m[l][0] = matrix[l][0][k];
        m[l][1] = matrix[l][1][k];
    }
    m
}

fn _const_external_product((a_bar, b_bar): Decomposition, m: [[u32; 2]; 2 * L]) -> (Ring, Ring) {
    let mut a_: Ring = [0; N];
    let mut b_: Ring = [0; N];
    for j in 0..N {
        let mut sa: u32 = 0;
        let mut sb: u32 = 0;
        for i in 0..L {
            sa = sa.wrapping_add(
                m[i][0]
                    .wrapping_mul(a_bar[i][j] as u32)
                    .wrapping_add(m[i + L][0].wrapping_mul(b_bar[i][j] as u32)),
            );
            sb = sb.wrapping_add(
                m[i][1]
                    .wrapping_mul(a_bar[i][j] as u32)
                    .wrapping_add(m[i + L][1].wrapping_mul(b_bar[i][j] as u32)),
            );
        }
        a_[j] = sa;
        b_[j] = sb;
    }
    (a_, b_)
}

#[cfg(test)]
mod tests {
    use super::super::super::trlwe::TRLWE;
    use super::super::super::util::params::trlwe::BRing;
    use super::super::super::util::sampling::random_bool_initialization;
    use super::super::{_decomposition, decomposition, BGBIT, L, N, TRGSW};

    #[test]
    fn test_decomposition_by_const_external_product() {
        let bs: BRing = random_bool_initialization();
        let trlwe = TRLWE::new();
        let (a, b) = trlwe.encrypt(bs);
        let a_bar = _decomposition(a);
        let b_bar = _decomposition(b);
        let mut matrix: [[u32; 2]; 2 * L] = [[0; 2]; 2 * L];
        for i in 0..L {
            let w = 2u32.pow(32 - (i as u32 + 1) * (BGBIT as u32));
            matrix[i][0] = w;
            matrix[i + L][1] = w;
        }
        let (a_, b_) = super::_const_external_product((a_bar, b_bar), matrix);
        assert_eq!(bs, trlwe.decrypt(a_, b_))
    }

    #[test]
    fn test_naive_external_product() {
        let bs: BRing = random_bool_initialization();
        let trlwe = TRLWE::new();
        let (a, b) = trlwe.encrypt(bs);
        let decomp = decomposition(a, b);
        let trgsw = TRGSW::new(trlwe.get_secret());
        let matrix = trgsw.coefficient(1);
        let (a_, b_) = super::_naive_external_product(decomp, matrix);
        let dec_bs = trlwe.decrypt(a_, b_);
        let mut counter = 0;
        for i in 0..N {
            counter += (bs[i] != dec_bs[i]) as usize;
        }
        assert!(counter == 0, "counter is {}", counter);
    }

    #[test]
    fn test_naive_external_product_without_zero_matrix_add() {
        let bs: BRing = random_bool_initialization();
        let trlwe = TRLWE::new();
        let (a, b) = trlwe.encrypt(bs);
        let decomp = decomposition(a, b);
        let trgsw = TRGSW::new(trlwe.get_secret());
        let matrix = trgsw._coefficient_without_zero_matrix_add(1);
        let (a_, b_) = super::_naive_external_product(decomp, matrix);
        assert_eq!(bs, trlwe.decrypt(a_, b_))
    }
}
