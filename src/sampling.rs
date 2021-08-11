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

    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = float_to_torus(normal.sample(&mut rng));
    }
    arr
}

pub fn ndim_bin_uniform<const N: usize>() -> [Torus; N] {
    let bin_uni = Uniform::new_inclusive(0, 1);
    let mut rng = rand::thread_rng();
    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = bin_uni.sample(&mut rng);
    }
    arr
}

pub fn ndim_torus_uniform<const N: usize>() -> [Torus; N] {
    let torus_uni = Uniform::new_inclusive(-0.5, 0.5);
    let mut rng = rand::thread_rng();
    let mut arr = [0; N];
    for i in 0..N {
        arr[i] = float_to_torus(torus_uni.sample(&mut rng));
    }
    arr
}

pub fn random_bool_initialization<const N: usize>() -> [bool; N] {
    let mut arr = [false; N];
    for i in 0..N {
        arr[i] = rand::random::<bool>();
    }
    arr
}
