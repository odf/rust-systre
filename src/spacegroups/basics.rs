use std::ops::{Rem, Mul, Neg};
use std::collections::HashSet;

use num_traits::{One, Zero};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::linear_algebra::{Scalar, extend_basis};
use crate::arithmetic::matrices::Matrix;


pub fn mod_z<T, CS>(op: AffineMap<T, CS>) -> AffineMap<T, CS>
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
        T: Clone + Eq + std::hash::Hash + Zero + One + Rem<Output=T>,
        CS: Clone + Eq + std::hash::Hash,
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
            let ab = mod_z(&a * b);
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
    ops: Vec<AffineMap<T, CS>>
}


impl<T, CS, CSP> PrimitiveSetting<T, CS, CSP>
    where
        T: Clone + Eq + Scalar + Neg<Output=T>,
        CS: Clone + Eq
{
    fn new(ops: &Vec<AffineMap<T, CS>>) -> Self {
        assert!(!ops.is_empty());

        let dim = ops[0].dim();
        let id = Matrix::identity(dim);
        let s0 = Vector::zero(dim);

        let mut cell = id.get_rows();

        for op in ops {
            if op.linear_matrix() == id && op.shift() != s0 {
                let s: Vec<_> = op.shift().into_iter().collect();
                extend_basis(&s, &mut cell);
            }
        }
        todo!()
    }
}