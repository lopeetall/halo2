use ark_bls12_381::{Fr, G1Affine, G1Projective};
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use rand::Rng;

pub fn random_affine_point() -> G1Affine {
    let mut rng = rand::thread_rng();
    let random_u64_scalar = Fr::rand(&mut rng);
    G1Affine::prime_subgroup_generator()
        .mul(random_u64_scalar.into_repr())
        .into_affine()
}

pub fn create_random_basis(n: usize) -> Vec<G1Affine> {
    (0..n).map(|_i| random_affine_point()).collect()
}

// inner produce of a vector of field elements and a basis of G1 affine points
pub fn inner_product_group(coeffs: &Vec<Fr>, basis: &[G1Affine]) -> G1Projective {
    coeffs
        .iter()
        .zip(basis)
        .map(|(c, b)| b.into_projective().mul(c.into_repr()))
        .sum::<G1Projective>()
}

// inner product of two vectors of field elements
pub fn inner_product_field(left: &Vec<Fr>, right: &Vec<Fr>) -> Fr {
    left.iter().zip(right).map(|(l, r)| *l * r).sum::<Fr>()
}

// inner product of two vectors of field elements
// uses the inverse of each element in `right`
pub fn inner_product_field_inverse(left: &Vec<Fr>, right: &Vec<Fr>) -> Fr {
    left.iter().zip(right).map(|(l, r)| *l * r).sum::<Fr>()
}

pub fn commit(poly: &DensePolynomial<Fr>, basis: &Vec<G1Affine>, r: Fr, h: G1Affine) -> G1Affine {
    (inner_product_group(&poly.coeffs().to_vec(), basis) + h.into_projective().mul(r.into_repr()))
        .into_affine()
}

// the blinding polynomial should have a root at the challenge x
// two blinding factors ought to be enough (but check this)
pub fn create_blinding_polynomial(
    x: Fr,
    b0: Fr,
    b1: Fr,
    public_poly: &DensePolynomial<Fr>,
) -> DensePolynomial<Fr> {
    // (X - x)(b_1x + b_0)(pub_poly)
    let vanishing_factor = DensePolynomial::from_coefficients_vec(vec![-x, Fr::one()]);
    let blinding_factor = DensePolynomial::from_coefficients_vec(vec![b0, b1]);
    vanishing_factor
        .naive_mul(&blinding_factor)
        .naive_mul(&public_poly)
}

pub fn naive_open(
    p_comm: G1Affine,
    s_comm: G1Affine,
    blinded_poly: &DensePolynomial<Fr>,
    basis: &Vec<G1Affine>,
    rand: Fr,
    h: G1Affine,
    p_x: Fr,
    s_b: Fr,
    iota: Fr,
) -> bool {
    let left = p_comm.into_projective() - basis[0].mul(p_x.into_repr()) + s_comm.mul(iota);
    let right = commit(&blinded_poly, &basis, rand + iota * s_b, h);

    left == right
}

/// Concentrates three vectors into a single element each.
/// The values here are shared between the Prover and Verifier
/// so any party can do this concentration
pub fn concentrate(
    a_poly: &DensePolynomial<Fr>,
    b: &Vec<Fr>,
    g: &Vec<G1Affine>,
    u: &Vec<Fr>,
    big_h: G1Affine,
    big_u: G1Affine,
    z: Fr,
    rand: Fr,
    s_b: Fr,
    iota: Fr,
) -> (G1Affine, Vec<G1Affine>, Vec<G1Affine>, Vec<Fr>, Vec<Fr>, Fr) {
    let mut rng = rand::thread_rng();
    let a = a_poly.coeffs();

    // check that all lengths are equal
    assert!(a.len() == b.len());
    assert!(g.len() == u.len());
    assert!(a.len() == g.len());

    // check that the lengths are all powers of 2
    let log = (a.len() as f64).log2().trunc() as usize;
    assert!(a.len() == 2usize.pow(log as u32));

    // Rather than compressing each vector using the log(n) u challenges
    // we can create a single length n vector out of the u challenges
    // which, if ordered correctly, gives the same result as an inner product.
    //
    // The ordering is inspired by binary counting: As an example, for n = 8
    // and log(n) = 3 we start with small basis [1, u0, u1, u2]. After expanding
    // we create the large basis: [1, u0, u1, u1u0, u2, u2u0, u2u1, u2u1u0]
    //
    // The large basis element at index i is the product of non-trivial
    // small basis elements selected by the binary representation of i,
    // with each set bit representing one of the log(n) non-trivial small
    // basis elements and each 0 representing 1.
    // For example, when i = 6 = 110 the large basis element is u2*u1*1
    //
    // create binary ordered u-basis for easy inner product
    let u_base = vec![Fr::rand(&mut rng); log];
    let mut u = vec![Fr::one()];
    for ui in u_base {
        let u_scaled = u.iter().map(|u| ui * u).collect::<Vec<Fr>>();
        u.extend(u_scaled);
    }

    let a_res = inner_product_field_inverse(&a.to_vec(), &u);
    let b_res = inner_product_field(&b, &u.to_vec());
    let g_res = inner_product_group(&u.to_vec(), g);

    // This portion may be simplifiable in the same way as the a, b, and g folds from above
    // I just haven't thought much about it yet.
    let mut big_l = vec![G1Affine::zero(); log];
    let mut big_r = vec![G1Affine::zero(); log];

    let mut size = a.len();

    let mut l = vec![Fr::zero(); log];
    let mut r = vec![Fr::zero(); log];

    let mut a_curr = a;
    let mut b_curr = b;
    let mut g_curr = g;

    let mut blind = rand + s_b * iota;

    for j in (0..log).rev() {
        size = size / 2;

        l[j] = Fr::rand(&mut rng);
        r[j] = Fr::rand(&mut rng);

        big_l[j] = (inner_product_group(&a_curr[size..].to_vec(), &g_curr[..size].to_vec())
            + big_u
                .mul(z * inner_product_field(&a_curr[size..].to_vec(), &b_curr[..size].to_vec()))
            + big_h.mul(l[j].into_repr()))
        .into_affine();

        big_r[j] = (inner_product_group(&a_curr[..size].to_vec(), &g_curr[size..].to_vec())
            + big_u
                .mul(z * inner_product_field(&a_curr[..size].to_vec(), &b_curr[size..].to_vec()))
            + big_h.mul(r[j].into_repr()))
        .into_affine();

        blind += l[j] * u[j].inverse().unwrap() + r[j] * u[j];
    }

    (g_res.mul(a_res.into_repr()).into_affine(), big_l, big_r, l, r, blind)
}
