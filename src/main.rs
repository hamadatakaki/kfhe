use kfhe::bootstrapping::gate_bootstrapping;
use kfhe::key::SecretKey;
use kfhe::ops::dot;
use kfhe::tlwe::{CipherTLWELv0, TLWE};
use kfhe::util::{RingLv1, Torus};

fn _nand(x: bool, y: bool) -> bool {
    let sk = SecretKey::new();
    let tlwe: TLWE = TLWE::new(sk);
    let c_x = tlwe.encrypt(x);
    let c_y = tlwe.encrypt(y);
    !tlwe.decrypt(c_x + c_y)
}

fn decrypt_as_tlwe_lv1(ext_a: RingLv1, ext_b: Torus, s: RingLv1) -> bool {
    let m = ext_b.wrapping_sub(dot(ext_a, s)).wrapping_sub(2u32.pow(28));
    m < 2u32.pow(31)
}

fn __nand(x: bool, y: bool) -> bool {
    let sk = SecretKey::new();
    let offset = CipherTLWELv0::clearly_true();

    let tlwe = TLWE::new(sk);
    let c1 = tlwe.encrypt(x);
    let c2 = tlwe.encrypt(y);
    let (a_, b_) = gate_bootstrapping(c1 + c2 - offset, sk).describe();
    decrypt_as_tlwe_lv1(a_, b_, sk.lv1)
}

fn main() {
    for i in 0..10 {
        println!("> loop {}", i);
        println!("{}", __nand(true, true));
        println!("{}", __nand(true, false));
        println!("{}", __nand(false, true));
        println!("{}", __nand(false, false));
    }
}

// test_gate_bootstrapping と同様の理由で通常は実行しない

#[test]
fn test_homnand() {
    use kfhe::key_switching::identity_key_switching;
    use kfhe::sampling::random_bool_initialization;
    use kfhe::tlwe::TLWE;

    let bs: [bool; 16] = random_bool_initialization();
    let mut count = 0;

    for b in bs {
        let sk = SecretKey::new();
        let tlwe = TLWE::new(sk);
        let c = tlwe.encrypt(b);
        let c1 = gate_bootstrapping(c, sk);
        let c0 = identity_key_switching(c1, sk);
        let msg = tlwe.decrypt(c0);

        count += (b != msg) as usize;
    }

    assert!(count == 0, "count: {}", count);
}
