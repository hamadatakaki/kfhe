pub mod tlwe {
    pub const N: usize = 635;
    pub const ALPHA: f64 = 3.0517578125e-05;
}

pub mod trlwe {
    pub const N: usize = 1024;
    pub const ALPHA: f64 = 2.98023223876953125e-08;
}

pub mod trgsw {
    pub const N: usize = super::trlwe::N;
    pub const BGBIT: usize = 6;
    pub const BG: usize = 64;
    pub const L: usize = 3;
}
