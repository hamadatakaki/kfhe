use ndarray::{Array, ArrayBase};
use ndarray_rand::rand_distr::{Distribution, Normal, Uniform};
use ndarray_rand::RandomExt;

const N: usize = 8;

fn sample_modular_normal_dist(mu: f64, alpha: f64) -> f64 {
    let normal = Normal::new(mu, alpha).unwrap();
    let sample = normal.sample(&mut rand::thread_rng());
    (sample + 0.5) % 1. - 0.5
}

#[derive(Clone, Debug)]
struct TLWE {
    s: Vec<f64>,
    alpha: f64,
}

impl TLWE {
    fn new() -> Self {
        let bin_uniform = Uniform::new_inclusive(0u64, 1u64);
        let s = Array::random(N, bin_uniform).map(|&x| x as f64);
        let alpha = 2f64.powi(-15);
        Self {
            s: s.to_vec(),
            alpha,
        }
    }

    fn encrypt(&self, msg: f64) -> (Vec<f64>, f64) {
        let m = (2. * msg - 1.) / 8.;
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let torus_uniform = Uniform::new_inclusive(-0.5f64, 0.5f64);
        let a = Array::random(N, torus_uniform);
        let e = sample_modular_normal_dist(0., self.alpha);

        let b = a.dot(&s) + m + e;
        return (a.to_vec(), b);
    }

    fn decrypt(&self, a_: Vec<f64>, b: f64) -> f64 {
        let a = Array::from_shape_vec(N, a_.clone()).unwrap();
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let m = b - a.dot(&s);
        (1. + m.signum()) / 2.
    }
}

fn run_tlwe(msg: f64) -> f64 {
    let tlwe = TLWE::new();
    let (a, b) = tlwe.encrypt(msg);
    tlwe.decrypt(a, b)
}

fn main() {
    let mut c = 0.;
    let mut m = 0.;
    loop {
        match c {
            10. => m += 1.,
            20. => break,
            _ => (),
        }
        println!("{}: {} -> {}", c, m, run_tlwe(m));
        c += 1.;
    }
}
