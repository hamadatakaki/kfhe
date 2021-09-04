use kfhe::key::SecretKey;
use kfhe::tlwe::TLWE;

fn nand(x: bool, y: bool) -> bool {
    let sk = SecretKey::new();
    let tlwe: TLWE = TLWE::new(sk);
    let c_x = tlwe.encrypt(x);
    let c_y = tlwe.encrypt(y);
    !tlwe.decrypt(c_x + c_y)
}

fn main() {
    println!("{}", nand(true, true));
    println!("{}", nand(true, false));
    println!("{}", nand(false, true));
    println!("{}", nand(false, false));
}
