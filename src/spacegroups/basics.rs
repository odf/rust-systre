use std::hash::Hash;
use std::ops::{Rem, Mul};
use std::collections::HashSet;

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap, ModZ};
use crate::arithmetic::linear_algebra::{extend_basis, LinearAlgebra};
use crate::arithmetic::matrices::Matrix;
use crate::spacegroups::lattices::rational_lattice_basis;


pub fn generate<T, CS>(gens: &Vec<AffineMap<T, CS>>)
    -> Vec<AffineMap<T, CS>>
    where
        T: Clone + Eq + Hash + Zero + One + Rem<Output=T>,
        CS: Clone + Eq + Hash,
        for <'a> &'a AffineMap<T, CS>:
            Mul<&'a AffineMap<T, CS>, Output=AffineMap<T, CS>>
{
    assert!(!gens.is_empty());

    let id = AffineMap::identity(gens[0].dim());
    let mut seen = HashSet::new();
    seen.insert(id.clone());
    let mut result = vec![id];
    let mut i = 0;

    while i < result.len() {
        let a = result[i].clone();
        i += 1;

        for b in gens {
            let ab = (&a * b).mod_z();
            if !seen.contains(&ab) {
                seen.insert(ab.clone());
                result.push(ab);
            }
        }
    }

    result
}


pub struct PrimitiveSetting<CS, CSP> {
    cell: Vec<Vector<BigRational, CS>>,
    to_primitive: CoordinateMap<BigRational, CS, CSP>,
    ops: Vec<AffineMap<BigRational, CSP>>
}


impl<CS, CSP> PrimitiveSetting<CS, CSP>
    where
        CS: Clone + Eq,
        CSP: Clone + Eq + Hash
{
    fn new(ops_all: &Vec<AffineMap<BigRational, CS>>) -> Self {
        assert!(!ops_all.is_empty());

        let dim = ops_all[0].dim();
        let id = Matrix::identity(dim);
        let s0 = Vector::zero(dim);

        let mut translations = vec![];
        for op in ops_all {
            if op.linear_matrix() == id && op.shift() != s0 {
                translations.push(op.shift().into_iter().collect());
            }
        }
        let cell: Vec<_> = rational_lattice_basis(&translations)
            .into_iter().map(|v| Vector::new(&v)).collect();
        let to_primitive = CoordinateMap::to_basis(&cell);

        let mut seen = HashSet::new();
        let mut ops = vec![];

        for op in ops_all {
            let t = (&to_primitive * op).mod_z();
            if !seen.contains(&t) {
                seen.insert(t.clone());
                ops.push(t);
            }
        }

        PrimitiveSetting { cell, to_primitive, ops }
    }
}


pub fn gram_matrix_configuration_space<CS: Clone>(
    ops: &Vec<AffineMap<BigRational, CS>>
)
    -> Option<Matrix<BigRational>>
{
    let dim = ops[0].dim();

    let mut k = Matrix::zero(dim, dim);
    let mut m = 0;

    for i in 0..dim {
        for j in i..dim {
            k[(i, j)] = m;
            k[(j, i)] = m;
            m += 1;
        }
    }

    let k = k;
    let m = m;

    let mut eqns = vec![];

    for op in ops {
        let s = op.linear_matrix();

        for ia in 0..dim {
            for ib in ia..dim {
                let mut eqn = vec![BigRational::zero(); m];
                for ja in 0..dim {
                    for jb in 0..dim {
                        eqn[k[(ja, jb)]] += &s[(ia, ja)] * &s[(ib, jb)];
                    }
                }
                eqn[k[(ia, ib)]] -= BigRational::one();
                extend_basis(&eqn, &mut eqns);
            }
        }
    }

    if eqns.is_empty() {
        Some(Matrix::identity(m))
    } else {
        let eqns: Vec<Matrix<_>> = eqns.iter().map(|v| Matrix::row(v)).collect();
        Matrix::vstack(&eqns).null_space()
    }
}


pub fn shift_space<CS: Clone>(ops: &Vec<AffineMap<BigRational, CS>>)
    -> Option<Matrix<BigRational>>
{
    let dim = ops[0].dim();
    let id = Matrix::identity(dim);

    let mut eqns = vec![];

    for op in ops {
        for row in (op.linear_matrix() - &id).get_rows() {
            extend_basis(&row, &mut eqns);
        }
    }

    let eqns: Vec<Matrix<_>> = eqns.iter().map(|v| Matrix::row(v)).collect();
    Matrix::vstack(&eqns).null_space()
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::*;

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    #[derive(Clone, Debug, Eq, PartialEq, Hash)]
    struct Standard {}

    #[derive(Clone, Debug, Eq, PartialEq, Hash)]
    struct Primitive {}

    #[test]
    fn test_spacegroup_basics_primitive_settings() {
        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(1)])),
            AffineMap::new(
                &Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
                &(Vector::new(&[r(1), r(1)]) / r(2))
            )
        ];
        let ops = generate(&gens);
        let primitive: PrimitiveSetting<Standard, Primitive>
            = PrimitiveSetting::new(&ops);

        let dim = ops[0].dim();
        let id_m = Matrix::identity(dim);
        let id_op = AffineMap::identity(dim);

        for op in ops {
            if op.linear_matrix() == id_m {
                assert_eq!((&primitive.to_primitive * &op).mod_z(), id_op);
            }
        }

        assert_eq!(primitive.to_primitive.determinant(), r(2));
    }

    #[test]
    fn test_spacegroup_basics_gram_matrix_configuration_space() {
        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(1)])),
            AffineMap::new(
                &Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
                &(Vector::new(&[r(1), r(1)]) / r(2))
            )
        ];
        assert_eq!(
            gram_matrix_configuration_space(&gens),
            Some(Matrix::new(2, &[r(1), r(0), r(0), r(0), r(0), r(1)]))
        );

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(0), r(1), r(-1), r(0)])),
        ];
        assert_eq!(
            gram_matrix_configuration_space(&gens),
            Some(Matrix::new(1, &[r(1), r(0), r(1)]))
        );

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(0), r(1), r(1), r(0)])),
        ];
        assert_eq!(
            gram_matrix_configuration_space(&gens),
            Some(Matrix::new(
                3, &[r(1), r(0), r(1), r(0), r(1), r(0)]
            ).transpose())
        );

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(-1)])),
        ];
        assert_eq!(
            gram_matrix_configuration_space(&gens),
            Some(Matrix::new(
                3, &[r(1), r(0), r(0), r(0), r(1), r(0), r(0), r(0), r(1)]
            ))
        );
    }

    #[test]
    fn test_spacegroup_basics_shift_space() {
        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(1)])),
            AffineMap::new(
                &Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
                &(Vector::new(&[r(1), r(1)]) / r(2))
            )
        ];
        assert_eq!(shift_space(&gens), Some(Matrix::new(1, &[r(0), r(1)])));

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(0), r(1), r(-1), r(0)])),
        ];
        assert_eq!(shift_space(&gens), None);

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(0), r(1), r(1), r(0)])),
        ];
        assert_eq!(shift_space(&gens), Some(Matrix::new(1, &[r(1), r(1)])));

        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(-1)])),
        ];
        assert_eq!(shift_space(&gens), None);
    }
}
