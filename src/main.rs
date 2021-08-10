use ndarray::Array;
use ndarray_rand::rand_distr::{Distribution, Normal, Uniform};
use ndarray_rand::RandomExt;

const N: usize = 8;

fn torus_to_ring(t: f64) -> u64 {
    let length_ring = 2u64.pow(32) as f64;
    let r = ((t + 0.5) * length_ring) as u64;
    r % 2u64.pow(32)
}

fn sample_modular_normal_dist(mu: f64, alpha: f64) -> u64 {
    let normal = Normal::new(mu, alpha).unwrap();
    let sample = normal.sample(&mut rand::thread_rng());
    torus_to_ring((sample + 0.5) % 1. - 0.5)
}

#[derive(Clone, Debug)]
struct TLWE {
    s: Vec<u64>,
    alpha: f64,
}

impl TLWE {
    fn new() -> Self {
        let bin_uniform = Uniform::new_inclusive(0u64, 1u64);
        let s = Array::random(N, bin_uniform).map(|&x| x as u64);
        let alpha = 2f64.powi(-15);
        Self {
            s: s.to_vec(),
            alpha,
        }
    }

    fn encrypt(&self, msg: u64) -> (Vec<u64>, u64) {
        let m = torus_to_ring((2. * (msg as f64) - 1.) / 8.);
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let torus_uniform = Uniform::new_inclusive(-0.5f64, 0.5f64);
        let a = Array::random(N, torus_uniform).map(|&x| torus_to_ring(x));
        let e = sample_modular_normal_dist(0., self.alpha);
        let b = a.dot(&s) + m + e;
        return (a.to_vec(), b);
    }

    fn decrypt(&self, a_: Vec<u64>, b: u64) -> u64 {
        let a = Array::from_shape_vec(N, a_.clone()).unwrap();
        let s = Array::from_shape_vec(N, self.s.clone()).unwrap();
        let m = (b - a.dot(&s) + 2u64.pow(29)) as i64;
        let m = m % 2i64.pow(32);
        let m = 2i64.pow(31) - m;
        ((1 + m.signum()) / 2) as u64
    }
}

fn run_tlwe(msg: u64) -> u64 {
    let tlwe = TLWE::new();
    let (a, b) = tlwe.encrypt(msg);
    tlwe.decrypt(a, b)
}

fn nand(x: u64, y: u64) -> u64 {
    let tlwe = TLWE::new();
    let (x_a, x_b) = tlwe.encrypt(x);
    let (y_a, y_b) = tlwe.encrypt(y);

    let a: Vec<u64> = x_a.iter().zip(y_a.iter()).map(|(x, y)| x + y).collect();
    let b = x_b + y_b;

    // debug
    // let s = Array::from_shape_vec(N, tlwe.s.clone()).unwrap();
    // let a_x = Array::from_shape_vec(N, x_a.clone()).unwrap();
    // let a_y = Array::from_shape_vec(N, y_a.clone()).unwrap();
    // let m_x = x_b - a_x.dot(&s);
    // let m_y = y_b - a_y.dot(&s);
    // let a_z = Array::from_shape_vec(N, a.clone()).unwrap();
    // let m_z = b - a_z.dot(&s);
    // println!("{:?}", (m_x, m_y, m_z, m_x + m_y == m_z));

    tlwe.decrypt(a, b)
}

fn main() {
    println!("{}", nand(1, 1));
    println!("{}", nand(1, 0));
    println!("{}", nand(0, 1));
    println!("{}", nand(0, 0));
}

#[test]
fn test_torus_to_ring() {
    assert_eq!(torus_to_ring(0.), 2u64.pow(31));
    assert_eq!(torus_to_ring(0.5 - 1. / 2f64.powi(32)), 2u64.pow(32) - 1);
    assert_eq!(torus_to_ring(-0.5), 0);
}
