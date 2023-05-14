use std::hash::Hash;
use std::ops::{Rem, Mul, Neg};
use std::collections::HashSet;

use num_traits::{One, Zero};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::linear_algebra::{Scalar, extend_basis, LinearAlgebra};
use crate::arithmetic::matrices::Matrix;


pub fn mod_z<T, CS>(op: &AffineMap<T, CS>) -> AffineMap<T, CS>
    where
        T: Clone + One + Rem<Output=T>,
        CS: Clone
{
    let s = op.shift().into_iter().map(|x| x % T::one()).collect();
    AffineMap::new(&op.linear_matrix(), &s)
}


pub fn generate<T, CS>(gens: Vec<AffineMap<T, CS>>)
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

        for b in &gens {
            let ab = mod_z(&(&a * b));
            if !seen.contains(&ab) {
                seen.insert(ab.clone());
                result.push(ab);
            }
        }
    }

    result
}


pub struct PrimitiveSetting<T, CS, CSP> {
    cell: Vec<Vector<T, CS>>,
    to_primitive: CoordinateMap<T, CS, CSP>,
    ops: Vec<AffineMap<T, CSP>>
}


impl<T, CS, CSP> PrimitiveSetting<T, CS, CSP>
    where
        T: Clone + Eq + Hash,
        T: Scalar + Neg<Output=T> + Rem<Output=T>,
        Matrix<T>: LinearAlgebra<T>,
        for <'a> &'a Matrix<T>: Mul<&'a Matrix<T>, Output=Matrix<T>>,
        for <'a> &'a AffineMap<T, CS>:
            Mul<&'a AffineMap<T, CS>, Output=AffineMap<T, CS>>,
        for <'a> AffineMap<T, CSP>:
            Mul<&'a AffineMap<T, CSP>, Output=AffineMap<T, CSP>>,
        CS: Clone + Eq,
        CSP: Clone + Eq + Hash
{
    fn new(ops_all: &Vec<AffineMap<T, CS>>) -> Self {
        assert!(!ops_all.is_empty());

        let dim = ops_all[0].dim();
        let id = Matrix::identity(dim);
        let s0 = Vector::zero(dim);

        let mut cell = id.get_rows();

        for op in ops_all {
            if op.linear_matrix() == id && op.shift() != s0 {
                let s: Vec<_> = op.shift().into_iter().collect();
                extend_basis(&s, &mut cell);
            }
        }

        let cell: Vec<_> = cell.into_iter().map(|v| Vector::new(&v)).collect();
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
