use std::hash::Hash;
use std::ops::{Rem, Mul, Neg};
use std::collections::HashSet;

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::linear_algebra::{Scalar, extend_basis, LinearAlgebra};
use crate::arithmetic::matrices::Matrix;
use crate::spacegroups::lattices::rational_lattice_basis;


pub fn mod_z<T, CS>(op: &AffineMap<T, CS>) -> AffineMap<T, CS>
    where
        T: Clone + One + Rem<Output=T>,
        CS: Clone
{
    let s = op.shift().into_iter().map(|x| x % T::one()).collect();
    AffineMap::new(&op.linear_matrix(), &s)
}


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
            let ab = mod_z(&(&a * b));
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
            let t = mod_z(&(&to_primitive * op));
            if !seen.contains(&t) {
                seen.insert(t.clone());
                ops.push(t);
            }
        }

        PrimitiveSetting { cell, to_primitive, ops }
    }
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
    fn test_spacegroup_primitive_settings() {
        let gens: Vec<AffineMap<_, Standard>> = vec![
            AffineMap::from(Matrix::new(2, &[r(-1), r(0), r(0), r(1)])),
            AffineMap::new(
                &Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
                &(Vector::new(&[r(1), r(1)]) / r(2))
            )
        ];
        let ops = generate(&gens);
        eprintln!("{}", ops.len());
        let primitive: PrimitiveSetting<Standard, Primitive>
            = PrimitiveSetting::new(&ops);

        // TODO test this without relying on arbitrary choices in the code
        assert_eq!(primitive.cell, vec![
            Vector::new(&[r(1), r(-1)]) / r(2),
            Vector::new(&[r(0), r(1)])
        ]);
    }
}
