use super::util::sampling::ndim_bin_uniform;
use super::util::{RingLv0, RingLv1};

#[derive(Clone, Copy, Debug)]
pub struct SecretKey {
    pub lv0: RingLv0,
    pub lv1: RingLv1,
}

impl SecretKey {
    pub fn new() -> Self {
        let lv0 = ndim_bin_uniform();
        let lv1 = ndim_bin_uniform();
        Self { lv0, lv1 }
    }
}
