mod commitment;

use crate::commitment::*;
use ark_poly::{
    Polynomial,
    univariate::DensePolynomial,
    UVPolynomial,
};
use ark_bls12_381::Fr;
use ark_ff::UniformRand;

use rand::Rng;

fn main() {
    let n = 8;
    let mut rng = rand::thread_rng();
    let basis = create_random_basis(n);
    let p_poly: DensePolynomial<Fr> = DensePolynomial::rand(n-1, &mut rng);

    // public poly is used to create the blinding polynomial
    // ensuring correct vanishing adds 1 to the degree
    // blinding adds 1 to the degree
    // so this must be degree n-3 to give a degree n-1 blinding polynomial
    let public_poly: DensePolynomial<Fr> = DensePolynomial::rand(n-3, &mut rng);
    let rand = Fr::rand(&mut rng);
    let h = random_affine_point();
    let p_comm = commit(&p_poly, &basis, rand, h);

    // verifier chooses this as a challenge
    let x = Fr::rand(&mut rng);

    // prover's blinding factors
    let b0 = Fr::rand(&mut rng);
    let b1 = Fr::rand(&mut rng);

    let s_poly = create_blinding_polynomial(x, b0, b1, &public_poly);
    let s_b = Fr::rand(&mut rng);
    let s_comm = commit(&s_poly, &basis, s_b, h);
    
    let iota = Fr::rand(&mut rng);
    let p_blinded = p_poly + s_poly.mul(iota);
}
