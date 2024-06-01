#![allow(non_snake_case)]

use curve25519_dalek::constants::{RISTRETTO_BASEPOINT_COMPRESSED, RISTRETTO_BASEPOINT_POINT};
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::MultiscalarMul;
use serde::{Deserialize, Serialize};
use sha3::Sha3_512;

use crate::util;

pub struct InitialState<'a> {
    pub n: usize,
    pub v: u64,
    pub v_blinding: Scalar,
    pub pc_gens: &'a PedersenGens,
    pub a_blinding: Scalar,
    pub s_blinding: Scalar,
    pub s_L: Vec<Scalar>,
    pub s_R: Vec<Scalar>,
}

pub struct ProofData {
    pub l_poly: util::VecPoly1,
    pub r_poly: util::VecPoly1,
    pub t_poly: util::Poly2,
    pub v_blinding: Scalar,
    pub a_blinding: Scalar,
    pub s_blinding: Scalar,
    pub t_1_blinding: Scalar,
    pub t_2_blinding: Scalar,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct FinalProof {
    pub(super) t_x: Scalar,
    pub(super) t_x_blinding: Scalar,
    pub(super) l_vec: Vec<Scalar>,
    pub(super) r_vec: Vec<Scalar>,
}

#[derive(Serialize, Deserialize, Copy, Clone, Debug)]
pub struct PolyCommitment {
    pub(super) T_1: RistrettoPoint,
    pub(super) T_2: RistrettoPoint,
}

#[derive(Serialize, Deserialize, Copy, Clone, Debug)]
pub struct BitCommitment {
    pub(super) V: CompressedRistretto,
}

/// Challenge values derived from all parties' [`BitCommitment`]s.
#[derive(Serialize, Deserialize, Copy, Clone, Debug)]
pub struct BitChallenge {
    pub(super) y: Scalar,
    pub(super) z: Scalar,
}

#[derive(Serialize, Deserialize, Copy, Clone, Debug)]
pub struct PolyChallenge {
    pub(super) x: Scalar,
}

#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct PedersenGens {
    pub B: RistrettoPoint,
    pub B_blinding: RistrettoPoint,
}

impl PedersenGens {
    pub fn commit(&self, value: Scalar, blinding: Scalar) -> RistrettoPoint {
        RistrettoPoint::multiscalar_mul(&[value, blinding], &[self.B, self.B_blinding])
    }
}

impl Default for PedersenGens {
    fn default() -> Self {
        PedersenGens {
            B: RISTRETTO_BASEPOINT_POINT,
            B_blinding: RistrettoPoint::hash_from_bytes::<Sha3_512>(
                RISTRETTO_BASEPOINT_COMPRESSED.as_bytes(),
            ),
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct VerifierData {
    pub proof: FinalProof,
    pub pedersen_gens: PedersenGens,
    pub bit_commitment: BitCommitment,
    pub poly_commitment: PolyCommitment,
    pub n: usize,
}
