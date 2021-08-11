use super::Torus;

pub fn vadd<const N: usize>(v: &[Torus; N], w: &[Torus; N]) -> [Torus; N] {
    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = v[i].wrapping_add(w[i]);
    }
    arr
}

pub fn dot<const N: usize>(v: &[Torus; N], w: &[Torus; N]) -> Torus {
    let mut s: Torus = 0;
    for i in 0..N {
        s = s.wrapping_add(v[i].wrapping_mul(w[i]));
    }
    s
}
