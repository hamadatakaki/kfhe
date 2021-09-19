pub mod tlwe {
    pub const N: usize = 635;
    pub const ALPHA: f64 = 3.051_757_812_5e-5;
}

pub mod trlwe {
    pub const N: usize = 1024;
    pub const ALPHA: f64 = 2.980_232_238_769_531_3e-8;
}

pub mod trgsw {
    pub const N: usize = super::trlwe::N;
    pub const BGBIT: u32 = 6;
    pub const BG: u32 = 64;
    pub const L: usize = 3;
    pub const NBIT: usize = 10;

    pub const T: usize = 8;
    pub const BASEBIT: u32 = 2;

    pub const SIGN_MIN: i8 = -((BG / 2) as i8);
    pub const SIGN_MAX: i8 = (BG / 2) as i8 - 1;
}
