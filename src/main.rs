use kfhe::tlwe::TLWE;
use kfhe::util::ops::vadd;

fn _run_tlwe(msg: bool) -> bool {
    let tlwe: TLWE = TLWE::new();
    let (a, b) = tlwe.encrypt(msg);
    tlwe.decrypt(a, b)
}

fn nand(x: bool, y: bool) -> bool {
    let tlwe: TLWE = TLWE::new();
    let (x_a, x_b) = tlwe.encrypt(x);
    let (y_a, y_b) = tlwe.encrypt(y);

    let a = vadd(&x_a, &y_a);
    let b = x_b.wrapping_add(y_b);

    !tlwe.decrypt(a, b)
}

fn main() {
    println!("{}", nand(true, true));
    println!("{}", nand(true, false));
    println!("{}", nand(false, true));
    println!("{}", nand(false, false));
}
