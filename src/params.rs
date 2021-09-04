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
    pub const BGBIT: u32 = 6;
    pub const BG: u32 = 64;
    pub const L: usize = 3;
    pub const NBIT: usize = 10;

    pub const SIGN_MIN: i8 = -((BG / 2) as i8);
    pub const SIGN_MAX: i8 = (BG / 2) as i8 - 1;
}