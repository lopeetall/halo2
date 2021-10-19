use ark_poly::{
    univariate::DensePolynomial,
    UVPolynomial,
};
use ark_bls12_381::{G1Affine, G1Projective, Fr};
use ark_ec::{ProjectiveCurve, AffineCurve};
use ark_ff::{PrimeField, One, Zero};
use rand::Rng;

// NOT UNIFORMLY RANDOM
pub fn bad_random_affine_point() -> G1Affine {
    let mut rng = rand::thread_rng();
    let random_u64_scalar = Fr::from(rng.gen::<u64>());
    G1Affine::prime_subgroup_generator().mul(random_u64_scalar.into_repr()).into_affine()
}

// NOT UNIFORMLY RANDOM
pub fn bad_create_random_basis(n: usize) -> Vec<G1Affine> {
    (0..n).map(|_i| bad_random_affine_point()).collect()
}

pub fn inner_product(coeffs: &[Fr], basis: &Vec<G1Affine>) -> G1Projective {
    coeffs
        .iter()
        .zip(basis)
        .map(|(c,b)|
            b
                .into_projective()
                .mul(c.into_repr())
            )
        .sum::<G1Projective>()          
}

pub fn commit(poly: &DensePolynomial<Fr>, basis: &Vec<G1Affine>, r: Fr, h: G1Affine) -> G1Affine {
    (
        inner_product(poly.coeffs(), basis) + h.into_projective().mul(r.into_repr())
    ).into_affine()
}

// the blinding polynomial should have a root at the challenge x
// two blinding factors ought to be enough (but check this)
pub fn create_blinding_polynomial(x: Fr, b0: Fr, b1: Fr, public_poly: &DensePolynomial<Fr>) -> DensePolynomial<Fr> {
    // (X - x)(b_1x + b_0)(pub_poly)
    let vanishing_factor = DensePolynomial::from_coefficients_vec(vec![-x, Fr::one()]);
    let blinding_factor = DensePolynomial::from_coefficients_vec(vec![b0, b1]);
    vanishing_factor.naive_mul(&blinding_factor).naive_mul(&public_poly)
}

pub fn naive_open(p_comm: G1Affine, s_comm: G1Affine, blinded_poly: &DensePolynomial<Fr>, basis: Vec<G1Affine>, rand: Fr, h: G1Affine, p_x: Fr, s_b: Fr, iota: Fr) -> bool {   
    let left = p_comm.into_projective() - basis[0].mul(p_x.into_repr()) + s_comm.mul(iota);
    let right = commit(&blinded_poly, &basis, rand+iota*s_b, h);

    left == right
}
