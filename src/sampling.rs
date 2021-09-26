use super::util::{float_to_torus, Torus};

use rand_distr::{Distribution, Normal, Uniform};

pub fn modular_normal_dist(mu: f64, alpha: f64) -> Torus {
    let normal = Normal::new(mu, alpha).unwrap();
    let sample = normal.sample(&mut rand::thread_rng());
    float_to_torus(sample)
}

pub fn ndim_modular_normal_dist<const N: usize>(mu: f64, alpha: f64) -> [Torus; N] {
    let normal = Normal::new(mu, alpha).unwrap();
    let mut rng = rand::thread_rng();
    let mut ret = [0; N];
    for i in 0..ret.len() {
        ret[i] = float_to_torus(normal.sample(&mut rng));
    }
    ret
}

pub fn ndim_bin_uniform<const N: usize>() -> [Torus; N] {
    let bin_uni = Uniform::new_inclusive(0, 1);
    let mut rng = rand::thread_rng();
    let mut ret = [0; N];
    for i in 0..ret.len() {
        ret[i] = bin_uni.sample(&mut rng);
    }
    ret
}

pub fn ndim_torus_uniform<const N: usize>() -> [Torus; N] {
    let torus_uni = Uniform::new_inclusive(-0.5, 0.5);
    let mut rng = rand::thread_rng();
    let mut ret = [0; N];
    for i in 0..ret.len() {
        ret[i] = float_to_torus(torus_uni.sample(&mut rng));
    }
    ret
}

pub fn random_bool_initialization<const N: usize>() -> [bool; N] {
    let mut ret = [true; N];
    for i in 0..ret.len() {
        ret[i] = rand::random::<bool>();
    }
    ret
}
