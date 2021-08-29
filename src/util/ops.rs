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

pub fn tscale(u: Torus, sc: i8) -> Torus {
    let v = if sc < 0 { (-sc) as Torus } else { sc as Torus };
    let base = u.wrapping_mul(v);
    if sc < 0 {
        0u32.wrapping_sub(base)
    } else {
        base
    }
}

pub fn intpoly_mul_as_torus<const N: usize>(zs: [i8; N], t: Torus) -> [Torus; N] {
    let mut ring = [0; N];
    for i in 0..N {
        ring[i] = tscale(t, zs[i]);
    }
    ring
}

pub fn rmadd<const L: usize, const M: usize, const N: usize>(
    a: &[[[Torus; N]; M]; L],
    b: &[[[Torus; N]; M]; L],
) -> [[[Torus; N]; M]; L] {
    let mut c = [[[0; N]; M]; L];
    for l in 0..L {
        for m in 0..M {
            c[l][m] = vadd(&a[l][m], &b[l][m]);
        }
    }
    c
}

pub fn rdot<const L: usize, const N: usize>(
    v: &[[Torus; N]; L],
    w: &[[Torus; N]; L],
) -> [Torus; N] {
    let mut s: [Torus; N] = [0; N];
    for l in 0..L {
        s = vadd(&s, &pmul(&v[l], &w[l]));
    }
    s
}
