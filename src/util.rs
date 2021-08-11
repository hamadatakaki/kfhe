pub mod ops;

use super::params::trlwe;

pub type Torus = u32;
pub type Ring = [Torus; trlwe::N];
pub type FRing = [f64; trlwe::N];
pub type BRing = [bool; trlwe::N];

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
    let length_ring = 2f64.powi(32);
    ((x + 0.5) * length_ring) as u32
}

pub fn fring_to_torus_ring<const N: usize>(xs: [f64; N]) -> [Torus; N] {
    let mut ring = [0; N];
    for i in 0..xs.len() {
        ring[i] = float_to_torus(xs[i]);
    }
    ring
}

#[test]
fn test_float_to_torus() {
    assert_eq!(float_to_torus(0.), 2u32.pow(31));
    assert_eq!(
        float_to_torus(0.5 - 1. / 2f64.powi(32)),
        2 * (2u32.pow(31) - 1) + 1
    );
    assert_eq!(float_to_torus(-0.5), 0);
}
