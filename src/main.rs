mod commitment;

use crate::commitment::*;
use ark_bls12_381::Fr;
use ark_ff::{Field, One, UniformRand, Zero, PrimeField};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_ec::{AffineCurve, ProjectiveCurve};

fn main() {
    let mut rng = rand::thread_rng();

    // setup
    let n = 8;

    let basis = create_random_basis(n);
    // public poly is used to create the blinding polynomial
    // ensuring correct vanishing adds 1 to the degree
    // blinding adds 1 to the degree
    // so this must be degree n-3 to give a degree n-1 blinding polynomial
    let public_poly: DensePolynomial<Fr> = DensePolynomial::rand(n - 3, &mut rng);
    let h = random_affine_point();

    // prover's polynomial
    let p_poly: DensePolynomial<Fr> = DensePolynomial::rand(n - 1, &mut rng);

    // prover commits with scalar r
    let rand = Fr::rand(&mut rng);
    let p_comm = commit(&p_poly, &basis, rand, h);

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
    assert!(s_poly.evaluate(&x) == Fr::zero());

    // prover commits to blinding poly with s_b
    let s_b = Fr::rand(&mut rng);
    let s_comm = commit(&s_poly, &basis, s_b, h);

    // prover creates blinded polynomial f_p
    let iota = Fr::rand(&mut rng);
    let f_p_poly = &(p_poly + &s_poly * iota) - &DensePolynomial::from_coefficients_vec(vec![p_x]);
    assert!(f_p_poly.evaluate(&x) == Fr::zero());

    assert!(naive_open(
        p_comm, s_comm, &f_p_poly, &basis, rand, h, p_x, s_b, iota
    ));

    // random elements for the concentrator
    let z = Fr::rand(&mut rng);
    let big_u = random_affine_point();

    // create list of powers of challenge x
    let mut c = Fr::one();
    let x_powers: Vec<Fr> = (0..n)
        .map(|i| {
            c *= x;
            c
        })
        .collect();

    // create list of challenges to use as u values for concentrator
    let u_challenges: Vec<Fr> = (0..n).map(|_i| Fr::rand(&mut rng)).collect();
    let u_challenges_inv = u_challenges.iter().map(|u| u.inverse().unwrap()).collect::<Vec<Fr>>();


    let (fp_g_compressed, big_l, big_r, l, r, blind) = concentrate(
        &f_p_poly,
        &x_powers,
        &basis,
        &u_challenges,
        big_u,
        h,
        z,
        rand,
        s_b,
        iota,
    );

    // reconstruct <f_p, G>
    
    let fp_g = inner_product_group(&f_p_poly.coeffs().to_vec(), &basis);
    println!("naive fp_g \n{:?}", fp_g.into_affine());
    let fp_g_log = fp_g_compressed.into_projective() + h.mul((blind - rand - iota*s_b).into_repr())-inner_product_group(&u_challenges_inv, &big_l) - inner_product_group(&u_challenges, &big_r);
    println!("log compressed fp_g\n{:?}", fp_g_log.into_affine());

    assert_eq!(fp_g, fp_g_log);
}
