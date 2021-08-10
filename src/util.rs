pub mod vec;

pub fn torus_to_ring(t: f64) -> u32 {
    let length_ring = 2f64.powi(32);
    let r = ((t + 0.5) * length_ring) as u32;
    r
}

#[test]
fn test_torus_to_ring() {
    assert_eq!(torus_to_ring(0.), 2u32.pow(31));
    assert_eq!(
        torus_to_ring(0.5 - 1. / 2f64.powi(32)),
        2 * (2u32.pow(31) - 1) + 1
    );
    assert_eq!(torus_to_ring(-0.5), 0);
}
