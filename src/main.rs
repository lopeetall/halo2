mod commitment;

use crate::commitment::*;
use ark_poly::{
    Polynomial,
    univariate::DensePolynomial,
    UVPolynomial,
};
use ark_bls12_381::Fr;
use ark_ff::{UniformRand, Zero};

use ark_ff::BigInteger256;
use ark_ec::AffineCurve;
use ark_ec::ProjectiveCurve;
//use rand::Rng;

fn main() {
    let mut rng = rand::thread_rng();

    // setup 
    let n = 8;

    let basis = bad_create_random_basis(n);
    // public poly is used to create the blinding polynomial
    // ensuring correct vanishing adds 1 to the degree
    // blinding adds 1 to the degree
    // so this must be degree n-3 to give a degree n-1 blinding polynomial
    let public_poly: DensePolynomial<Fr> = DensePolynomial::rand(n-3, &mut rng);
    let h = bad_random_affine_point();

    // prover's polynomial
    let p_poly: DensePolynomial<Fr> = DensePolynomial::rand(n-1, &mut rng);

    // prover commits with scalar r
    let r = Fr::rand(&mut rng);
    let p_comm = commit(&p_poly, &basis, r, h);

    // verifier chooses this as a challenge
    let x = Fr::rand(&mut rng);

    // prover evaluates their polynomial at x
    let p_x = p_poly.evaluate(&x);

    // prover chooses blinding factors
    let b0 = Fr::rand(&mut rng);
    let b1 = Fr::rand(&mut rng);

    // prover uses blinding factors and a public polynomial to 
    // create a degree n-1 blinding polynomial s with s(x) == 0
    let s_poly = create_blinding_polynomial(x, b0, b1, &public_poly);
    assert!(s_poly.evaluate(&x)==Fr::zero());

    // prover commits to blinding poly with s_b
    let s_b = Fr::rand(&mut rng);
    let s_comm = commit(&s_poly, &basis, s_b, h);
   
    // prover creates blinded polynomial f_p
    let iota = Fr::rand(&mut rng);
    let f_p_poly = &(p_poly + &s_poly*iota) - &DensePolynomial::from_coefficients_vec(vec![p_x]);
    assert!(f_p_poly.evaluate(&x)==Fr::zero());

    assert!(naive_open(p_comm, s_comm, &f_p_poly, basis, r, h, p_x, s_b, iota))
}
