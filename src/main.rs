use kfhe::tlwe::TLWE;
use kfhe::util::vec::vadd;

fn _run_tlwe(msg: u32) -> u32 {
    let tlwe = TLWE::new();
    let (a, b) = tlwe.encrypt(msg);
    tlwe.decrypt(a, b)
}

fn nand(x: u32, y: u32) -> u32 {
    let tlwe = TLWE::new();
    let (x_a, x_b) = tlwe.encrypt(x);
    let (y_a, y_b) = tlwe.encrypt(y);

    let a: Vec<u32> = vadd(&x_a, &y_a);
    let b = x_b.wrapping_add(y_b);

    tlwe.decrypt(a, b)
}

fn main() {
    println!("{}", nand(1, 1));
    println!("{}", nand(1, 0));
    println!("{}", nand(0, 1));
    println!("{}", nand(0, 0));
}
