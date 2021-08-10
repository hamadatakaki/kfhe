use super::util::torus_to_ring;
use rand_distr::{Distribution, Normal, Uniform};

pub fn modular_normal_dist(mu: f64, alpha: f64) -> u32 {
    let normal = Normal::new(mu, alpha).unwrap();
    let sample = normal.sample(&mut rand::thread_rng());
    torus_to_ring((sample + 0.5) % 1. - 0.5)
}

pub fn ndim_bin_uniform(n: usize) -> Vec<u32> {
    let bin_uni = Uniform::new_inclusive(0, 1);
    (0..n)
        .map(|_| bin_uni.sample(&mut rand::thread_rng()))
        .collect()
}

pub fn ndim_torus_uniform(n: usize) -> Vec<u32> {
    let torus_uni = Uniform::new_inclusive(-0.5, 0.5);
    (0..n)
        .map(|_| {
            let r = torus_uni.sample(&mut rand::thread_rng());
            torus_to_ring(r)
        })
        .collect()
}
