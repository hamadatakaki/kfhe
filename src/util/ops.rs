use super::params::Torus;

pub fn vadd<const N: usize>(v: &[Torus; N], w: &[Torus; N]) -> [Torus; N] {
    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = v[i].wrapping_add(w[i]);
    }
    arr
}

pub fn vsub<const N: usize>(v: &[Torus; N], w: &[Torus; N]) -> [Torus; N] {
    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = v[i].wrapping_sub(w[i]);
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

pub fn pmul<const N: usize>(p: &[Torus; N], q: &[Torus; N]) -> [Torus; N] {
    let mut mul: [Torus; N] = [0; N];
    for i in 0..N {
        for j in 0..N {
            if i + j < N {
                mul[i + j] = mul[i + j].wrapping_add(p[i].wrapping_mul(q[j]));
            } else {
                mul[i + j - N] = mul[i + j - N].wrapping_sub(p[i].wrapping_mul(q[j]));
            }
        }
    }
    mul
}
