use kfhe::tlwe::tlwe_nand;

fn main() {
    for i in 0..10 {
        println!("> loop {}", i);
        println!("{}", tlwe_nand(true, true));
        println!("{}", tlwe_nand(true, false));
        println!("{}", tlwe_nand(false, true));
        println!("{}", tlwe_nand(false, false));
    }
}

// test_gate_bootstrapping と同様の理由で通常は実行しない

#[test]
fn test_homnand() {
    use kfhe::homnand::homnand;
    use kfhe::key::SecretKey;
    use kfhe::sampling::random_bool_initialization;
    use kfhe::tlwe::TLWE;

    const S: usize = 4;
    let mut count = 0;

    let sk = SecretKey::new();
    let tlwe = TLWE::new(sk);
    let bs1: [bool; S] = random_bool_initialization();
    let bs2: [bool; S] = random_bool_initialization();

    for i in 0..bs1.len() {
        let b1 = bs1[i];
        let b2 = bs2[i];
        let nand = !(b1 && b2);

        let c1 = tlwe.encrypt(b1);
        let c2 = tlwe.encrypt(b2);
        let c = homnand(c1, c2, sk);
        let msg = tlwe.decrypt(c);

        count += (msg ^ nand) as usize;
    }

    assert!(count == 0, "count: {}", count);
}
