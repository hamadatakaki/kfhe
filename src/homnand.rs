use super::bootstrapping::gate_bootstrapping;
use super::key::SecretKey;
use super::key_switching::identity_key_switching;
use super::tlwe::CipherTLWELv0;

pub fn homnand(x: CipherTLWELv0, y: CipherTLWELv0, sk: SecretKey) -> CipherTLWELv0 {
    let c_true = CipherTLWELv0::clearly_true();
    let c = gate_bootstrapping(c_true - (x + y), sk);
    identity_key_switching(c, sk)
}
