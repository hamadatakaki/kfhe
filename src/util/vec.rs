pub fn vadd(v: &Vec<u32>, w: &Vec<u32>) -> Vec<u32> {
    v.iter()
        .zip(w.iter())
        .map(|(&x, &y)| x.wrapping_add(y))
        .collect()
}

pub fn dot(v: &Vec<u32>, w: &Vec<u32>) -> u32 {
    v.iter()
        .zip(w.iter())
        .fold(0, |acc, (&x, &y)| acc.wrapping_add(x.wrapping_mul(y)))
}
