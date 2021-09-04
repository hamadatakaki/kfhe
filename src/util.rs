use super::params;

pub type Torus = u32;
pub type RingLv0 = [Torus; params::tlwe::N];
pub type RingLv1 = [Torus; params::trlwe::N];

pub fn bool_normalization(b: bool) -> f64 {
    (2. * ((b as u8) as f64) - 1.) / 8.
}

pub fn boolpoly_normalization<const N: usize>(bs: [bool; N]) -> [f64; N] {
    let mut xs = [0.; N];
    for i in 0..N {
        xs[i] = bool_normalization(bs[i]);
    }
    xs
}

pub fn float_to_torus(x: f64) -> Torus {
    // [-0.5, 0.5) to Torus(u32)
    //      [-0.5, 0) to [2^31, 2^32)
    //      [0, 0.5)  to [0, 2^31)

    assert!(x >= -0.5);
    assert!(x < 0.5);

    let length_ring = 2f64.powi(32);
    ((if x < 0. { x + 1. } else { x }) * length_ring) as u32
}

pub fn torus_to_float(t: Torus) -> f64 {
    if t >= 2u32.pow(31) {
        (t as f64 / 2f64.powi(32)) - 1.
    } else {
        t as f64 / 2f64.powi(32)
    }
}

pub fn fring_to_torus_ring<const N: usize>(xs: [f64; N]) -> [Torus; N] {
    let mut ring = [0; N];
    for i in 0..xs.len() {
        ring[i] = float_to_torus(xs[i]);
    }
    ring
}

pub fn zpoly_to_ring<const N: usize>(zp: [i8; N]) -> [Torus; N] {
    let mut r = [0; N];
    for i in 0..N {
        r[i] = zp[i] as u32;
    }
    r
}

pub fn rotate_ring<const N: usize>(ring: [Torus; N], k: usize) -> [Torus; N] {
    let mut r = [0; N];
    for i in 0..N {
        r[i] = if i + k < N {
            ring[i + k]
        } else {
            0u32.wrapping_sub(ring[N - k + i])
        };
    }
    r
}

#[test]
fn test_float_to_torus() {
    assert_eq!(float_to_torus(0.), 0);
    assert_eq!(float_to_torus(0.5 - 1. / 2f64.powi(32)), 2u32.pow(31) - 1);
    assert_eq!(float_to_torus(-0.5), 2u32.pow(31));
}

#[test]
fn test_torus_to_float() {
    assert_eq!(torus_to_float(float_to_torus(0.)), 0.);
    assert_eq!(torus_to_float(float_to_torus(-0.5)), -0.5);
}
