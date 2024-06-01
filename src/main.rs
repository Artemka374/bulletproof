#![allow(non_snake_case)]

use std::iter;

use curve25519_dalek::scalar::Scalar;
use rand::thread_rng;

use types::*;

use crate::util::{exp_iter, inner_product};

mod types;
mod util;

pub struct RangeProof {}

impl RangeProof {
    #[allow(unreachable_code)]
    pub fn new(
        pc_gens: &PedersenGens,
        v: u64,
        v_blinding: Scalar,
        n: usize,
    ) -> (InitialState, BitCommitment) {
        let V = pc_gens.commit(v.into(), v_blinding).compress();

        let rng = &mut thread_rng();

        // Генеруємо приховуючі фактори a_blinding, s_blinding
        let a_blinding = Scalar::random(rng);
        let s_blinding = Scalar::random(rng);

        // Генеруємо випадкові числа s_L, s_R
        let s_L: Vec<Scalar> = (0..n).map(|_| Scalar::random(rng)).collect();
        let s_R: Vec<Scalar> = (0..n).map(|_| Scalar::random(rng)).collect();

        let bit_commitment = BitCommitment { V };

        let next_state = InitialState {
            n,
            v,
            v_blinding,
            pc_gens,
            a_blinding,
            s_blinding,
            s_L,
            s_R,
        };
        (next_state, bit_commitment)
    }

    pub fn generate_proof_data(
        data: InitialState,
        vc: &BitChallenge,
    ) -> (ProofData, PolyCommitment) {
        let n = data.n;

        let mut l_poly = util::VecPoly1::zero(n);
        let mut r_poly = util::VecPoly1::zero(n);

        let mut exp_y = Scalar::one(); // start at y^0 = 1
        let mut exp_2 = Scalar::one(); // start at 2^0 = 1

        // Рахуємо поліноми l, r
        for i in 0..n {
            // Рахуємо і-й біт числа v
            let a_L_i = Scalar::from((data.v >> i) & 1);
            // Рахуємо і-й біт - 1 числа v
            let a_R_i = a_L_i - Scalar::one();

            // За формулою l(x) = (a_L + s_L * x) - z1
            // l_poly.0 - частина не залежна від х
            // l_poly.1 - частина залежна від х
            l_poly.0[i] = a_L_i - vc.z;
            l_poly.1[i] = data.s_L[i];

            // За формулою r(x) = y^n * (a_R + s_R * x) + z^2 * 2^n
            // r_poly.0 - частина не залежна від х
            // r_poly.1 - частина залежна від х
            r_poly.0[i] = exp_y * (a_R_i + vc.z) + vc.z * vc.z * exp_2;
            r_poly.1[i] = exp_y * data.s_R[i];

            exp_y *= vc.y; // y^i -> y^(i+1)
            exp_2 = exp_2 + exp_2; // 2^i -> 2^(i+1)
        }

        // Рахуємо поліном t = l * r
        let t_poly = l_poly.inner_product(&r_poly);

        let rng = &mut thread_rng();

        // Генеруємо зобовʼязання T1, T2
        let t_1_blinding = Scalar::random(rng);
        let t_2_blinding = Scalar::random(rng);
        let T_1 = data.pc_gens.commit(t_poly.1, t_1_blinding);
        let T_2 = data.pc_gens.commit(t_poly.2, t_2_blinding);

        let poly_commitment = PolyCommitment { T_1, T_2 };

        let next_state = ProofData {
            v_blinding: data.v_blinding,
            a_blinding: data.a_blinding,
            s_blinding: data.s_blinding,
            l_poly,
            r_poly,
            t_poly,
            t_1_blinding,
            t_2_blinding,
        };

        (next_state, poly_commitment)
    }

    pub fn get_final_proof(prover: ProofData, bc: &BitChallenge, pc: &PolyChallenge) -> FinalProof {
        // Рахуємо зобовʼязання від t(x)
        let t_blinding_poly = util::Poly2(
            bc.z * bc.z * prover.v_blinding,
            prover.t_1_blinding,
            prover.t_2_blinding,
        );

        // підставляємо х в поліноми
        let t_x = prover.t_poly.eval(pc.x);
        let t_x_blinding = t_blinding_poly.eval(pc.x);
        let l_vec = prover.l_poly.eval(pc.x);
        let r_vec = prover.r_poly.eval(pc.x);

        FinalProof {
            t_x_blinding,
            t_x,
            l_vec,
            r_vec,
        }
    }

    pub fn verify_proof(
        proof_data: VerifierData,
        bit_challenge: BitChallenge,
        poly_challenge: PolyChallenge,
    ) -> anyhow::Result<()> {
        // Перевіряємо скалярний добуток l(x) * r(x) = t(x)
        let mult = inner_product(&proof_data.proof.l_vec, &proof_data.proof.r_vec);

        if mult != proof_data.proof.t_x {
            return Err(anyhow::anyhow!("Polynoms multiplication failed"));
        }

        // Генеруємо y^n, 2^n, 1^n
        let y_pow_n: Vec<Scalar> = exp_iter(bit_challenge.y).take(proof_data.n).collect();
        let two_pow_n: Vec<Scalar> = exp_iter(Scalar::from(2u64)).take(proof_data.n).collect();
        let one_n: Vec<Scalar> = iter::repeat(Scalar::one()).take(proof_data.n).collect();

        // (z - z^2) * mul(1, y^n) - z^3*mul(1, 2^n)
        let delta = (bit_challenge.z - bit_challenge.z * bit_challenge.z)
            * inner_product(&one_n, &y_pow_n)
            - bit_challenge.z
                * bit_challenge.z
                * bit_challenge.z
                * inner_product(&one_n, &two_pow_n);

        // Перевіряємо константність за формулою
        // t_x * B + t_x_blinding * B_blinding = z^2 * V + ((z - z^2) * mul(1, y^n) - z^3*mul(1, 2^n)) * B + x * T_1 + x^2 * T_2

        // t_x * B + t_x_blinding * B_blinding
        let left = proof_data.proof.t_x * proof_data.pedersen_gens.B
            + proof_data.proof.t_x_blinding * proof_data.pedersen_gens.B_blinding;

        // z^2 * V + ((z - z^2) * mul(1, y^n) - z^3*mul(1, 2^n)) * B + x * T_1 + x^2 * T_2
        let right =
            bit_challenge.z * bit_challenge.z * proof_data.bit_commitment.V.decompress().unwrap()
                + proof_data.pedersen_gens.B * delta
                + poly_challenge.x * proof_data.poly_commitment.T_1
                + poly_challenge.x * poly_challenge.x * proof_data.poly_commitment.T_2;

        if left != right {
            return Err(anyhow::anyhow!("Constant proof verification failed"));
        }

        Ok(())
    }
}

pub fn generate_proof(
    v: u64,
    n: usize,
    bit_challenge: BitChallenge,
    poly_challenge: PolyChallenge,
) -> VerifierData {
    let pc_gens = PedersenGens::default();

    let v_blinding = Scalar::random(&mut thread_rng());

    tracing::debug!(
        "Initializing prover with v = {}, n = {}, v_blinding = {:?}",
        v,
        n,
        v_blinding
    );

    let (initial_state, bit_commitment) =
        RangeProof::new(&pc_gens, v, Scalar::random(&mut thread_rng()), n);

    tracing::debug!("Received bit commitment: {:?}", bit_commitment);
    tracing::debug!("Continuing with bit challenge: {:?}", bit_challenge);

    let (product_proof_state, poly_commitment) =
        RangeProof::generate_proof_data(initial_state, &bit_challenge);

    tracing::debug!("Received poly commitment: {:?}", poly_commitment);
    tracing::debug!("Continuing with poly challenge: {:?}", poly_challenge);

    let final_proof =
        RangeProof::get_final_proof(product_proof_state, &bit_challenge, &poly_challenge);

    tracing::debug!("Final proof: {:?}", final_proof);

    VerifierData {
        proof: final_proof,
        pedersen_gens: pc_gens,
        bit_commitment,
        poly_commitment,
        n,
    }
}

fn main() {
    let v = 5;
    let n = 8;

    let bit_challenge = BitChallenge {
        y: Scalar::random(&mut thread_rng()),
        z: Scalar::random(&mut thread_rng()),
    };

    let poly_challenge = PolyChallenge {
        x: Scalar::random(&mut thread_rng()),
    };

    let proof_data = generate_proof(v, n, bit_challenge, poly_challenge);

    match RangeProof::verify_proof(proof_data, bit_challenge, poly_challenge) {
        Ok(_) => println!("Proof is valid"),
        Err(e) => println!("Proof is invalid: {:?}", e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correct_proofs() {
        let rng = &mut thread_rng();

        let bit_challenge_1 = BitChallenge {
            y: Scalar::random(rng),
            z: Scalar::random(rng),
        };
        let bit_challenge_2 = BitChallenge {
            y: Scalar::random(rng),
            z: Scalar::random(rng),
        };

        let poly_challenge_1 = PolyChallenge {
            x: Scalar::random(rng),
        };
        let poly_challenge_2 = PolyChallenge {
            x: Scalar::random(rng),
        };

        let proof_data = generate_proof(5, 3, bit_challenge_1, poly_challenge_1);
        assert!(RangeProof::verify_proof(proof_data, bit_challenge_1, poly_challenge_1).is_ok());

        let proof_data = generate_proof(10, 4, bit_challenge_2, poly_challenge_2);
        assert!(RangeProof::verify_proof(proof_data, bit_challenge_2, poly_challenge_2).is_ok());

        let proof_data = generate_proof(20, 5, bit_challenge_1, poly_challenge_2);
        assert!(RangeProof::verify_proof(proof_data, bit_challenge_1, poly_challenge_2).is_ok());

        let proof_data = generate_proof(40, 6, bit_challenge_2, poly_challenge_1);
        assert!(RangeProof::verify_proof(proof_data, bit_challenge_2, poly_challenge_1).is_ok());
    }

    #[test]
    fn test_incorrect_challenges() {
        let rng = &mut thread_rng();

        let bit_challenge_1 = BitChallenge {
            y: Scalar::random(rng),
            z: Scalar::random(rng),
        };
        let bit_challenge_2 = BitChallenge {
            y: Scalar::random(rng),
            z: Scalar::random(rng),
        };

        let poly_challenge_1 = PolyChallenge {
            x: Scalar::random(rng),
        };
        let poly_challenge_2 = PolyChallenge {
            x: Scalar::random(rng),
        };

        let mut proof_data = generate_proof(5, 3, bit_challenge_1, poly_challenge_1);

        assert!(
            RangeProof::verify_proof(proof_data.clone(), bit_challenge_1, poly_challenge_1).is_ok()
        );

        proof_data.proof.t_x = Scalar::random(rng);

        assert!(RangeProof::verify_proof(proof_data, bit_challenge_1, poly_challenge_1).is_err());

        proof_data = generate_proof(10, 4, bit_challenge_2, poly_challenge_2);
        assert!(
            RangeProof::verify_proof(proof_data.clone(), bit_challenge_2, poly_challenge_2).is_ok()
        );
        assert!(
            RangeProof::verify_proof(proof_data.clone(), bit_challenge_1, poly_challenge_2)
                .is_err()
        );
        assert!(
            RangeProof::verify_proof(proof_data.clone(), bit_challenge_2, poly_challenge_1)
                .is_err()
        );
        assert!(
            RangeProof::verify_proof(proof_data.clone(), bit_challenge_1, poly_challenge_1)
                .is_err()
        );
    }

    #[test]
    fn test_incorrect_proofs() {
        let rng = &mut thread_rng();

        let bit_challenge_1 = BitChallenge {
            y: Scalar::random(rng),
            z: Scalar::random(rng),
        };
        let poly_challenge_1 = PolyChallenge {
            x: Scalar::random(rng),
        };

        let mut proof_data = generate_proof(5, 3, bit_challenge_1, poly_challenge_1);
        proof_data.proof.t_x = Scalar::random(rng);
        assert!(RangeProof::verify_proof(proof_data, bit_challenge_1, poly_challenge_1).is_err());
    }
}
