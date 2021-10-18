use ark_poly::{
    Polynomial,
    univariate::DensePolynomial,
    UVPolynomial,
    EvaluationDomain,
    GeneralEvaluationDomain
};
use ark_bls12_381::{G1Affine, G1Projective, Fr};
use ark_ec::{ProjectiveCurve, AffineCurve};
use ark_ff::{UniformRand, PrimeField, One};
use rand::Rng;

// NOT CONSTANT TIME
pub fn random_affine_point() -> G1Affine {
    let mut rng = rand::thread_rng();
    let mut random_point = G1Affine::from_random_bytes(&[rng.gen::<u8>(); 32]);
    match random_point {
        None => random_point = Some(random_affine_point()),
        _ => {},
    }
    random_point.unwrap()
}

// NOT CONSTANT TIME
pub fn create_random_basis(n: usize) -> Vec<G1Affine> {
    let rng = rand::thread_rng();
    (0..n).map(|_i| random_affine_point()).collect()
}

pub fn commit(poly: &DensePolynomial<Fr>, basis: &Vec<G1Affine>, r: Fr, h: G1Affine) -> G1Affine {
    (poly
        .coeffs()
        .iter()
        .zip(basis)
        .map(|(c,b)| 
            b
                .into_projective()
                .mul(c.into_repr())
            )
        .sum::<G1Projective>()
    + h.into_projective().mul(r.into_repr())).into_affine()
}

// the blinding polynomial should have a root at the challenge x
// two blinding factors ought to be enough (but check this)
pub fn create_blinding_polynomial(x: Fr, b0: Fr, b1: Fr, public_poly: &DensePolynomial<Fr>) -> DensePolynomial<Fr> {
    // (X - x)(b_1x + b_0)(pub_poly)
    let vanishing_factor = DensePolynomial::from_coefficients_vec(vec![-x, Fr::one()]);
    let blinding_factor = DensePolynomial::from_coefficients_vec(vec![b0, b1]);
    vanishing_factor.naive_mul(&blinding_factor).naive_mul(&public_poly)
}

pub fn naive_open(comm_p: G1Affine, comm_s: G1Affine, blinded_poly: &DensePolynomial<Fr>, basis: Vec<G1Affine>, rand: Fr, h: G1Affine, x: Fr, s_b: Fr) -> Bool {
    
}
