use ndarray::Array;
use ndarray_rand::rand_distr::{Distribution, Normal, Uniform};
use ndarray_rand::RandomExt;

const N: usize = 8;

fn torus_to_ring(t: f64) -> i64 {
    let max_ring = 2f64.powi(31);
    (((t + 0.5) * 2. * max_ring) - max_ring) as i64
}

fn sample_modular_normal_dist(mu: f64, alpha: f64) -> i64 {
    let normal = Normal::new(mu, alpha).unwrap();
    let sample = normal.sample(&mut rand::thread_rng());
    torus_to_ring((sample + 0.5) % 1. - 0.5)
}

#[derive(Clone, Debug)]
struct TLWE {
    s: Vec<i64>,
    alpha: f64,
}

impl TLWE {
    fn new() -> Self {
        let bin_uniform = Uniform::new_inclusive(0u64, 1u64);
        let s = Array::random(N, bin_uniform).map(|&x| x as i64);
        let alpha = 2f64.powi(-15);
        Self {
            s: s.to_vec(),
            alpha,
        }
    }

    fn encrypt(&self, msg: i64) -> (Vec<i64>, i64) {
        let m = torus_to_ring((2. * (msg as f64) - 1.) / 8.);
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let torus_uniform = Uniform::new_inclusive(-0.5f64, 0.5f64);
        let a = Array::random(N, torus_uniform).map(|&x| torus_to_ring(x));
        let e = sample_modular_normal_dist(0., self.alpha);

        let b = a.dot(&s) + m + e;
        return (a.to_vec(), b);
    }

    fn decrypt(&self, a_: Vec<i64>, b: i64) -> i64 {
        let a = Array::from_shape_vec(N, a_.clone()).unwrap();
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let m = b - a.dot(&s);
        (1 + m.signum()) / 2
    }
}

fn run_tlwe(msg: i64) -> i64 {
    let tlwe = TLWE::new();
    let (a, b) = tlwe.encrypt(msg);
    tlwe.decrypt(a, b)
}

fn main() {
    let mut c = 0;
    let mut m = 0;
    loop {
        match c {
            10 => m += 1,
            20 => break,
            _ => (),
        }
        println!("{}: {} -> {}", c, m, run_tlwe(m));
        c += 1;
    }
}

#[test]
fn test_torus_to_ring() {
    assert_eq!(torus_to_ring(0.), 0);
    assert_eq!(torus_to_ring(0.5), 2i64.pow(31) - 1);
    assert_eq!(torus_to_ring(-0.5), -2i64.pow(31));
}
