use std::ops::{SubAssign, AddAssign, Mul, Neg, MulAssign};

use num_rational::BigRational;
use num_traits::{Signed, Zero, One};

use crate::arithmetic::matrices::Matrix;
use crate::arithmetic::linear_algebra::LinearAlgebra;


enum CrystalSystem2d {
    Oblique,
    Rectangular,
    Square,
    Hexagonal,
}

enum CrystalSystem3d {
    Cubic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}


fn normalized<T>(v: &[T]) -> Vec<T>
    where
        T: Clone + Signed,
        for <'a> &'a T: Neg<Output=T>
{
    if let Some(x) = v.iter().find(|x| !x.is_zero()) {
        if x.is_negative() {
            return v.iter().map(|x| -x).collect();
        }
    }
    v.to_vec()
}


fn dot<T>(v: &[T], w: &[T]) -> T
    where
        T: Zero,
        for <'a> &'a T: Mul<&'a T, Output=T>
{
    assert_eq!(v.len(), w.len());

    let mut s = T::zero();
    for i in 0..v.len() {
        s = s + &v[i] * &w[i];
    }

    s
}


fn are_orthogonal<T>(v: &[T], w: &[T]) -> bool
    where
        T: PartialEq + Zero,
        for <'a> &'a T: Mul<&'a T, Output=T>
{
    dot(v, w) == T::zero()
}


fn are_collinear<T>(v: &[T], w: &[T]) -> bool
    where
        T: PartialEq + Zero + Signed + Mul<T, Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T>
{
    dot(v, w) * dot(v, w) == dot(v, v) * dot(w, w)
}


fn matrix_order<T>(matrix: &Matrix<T>, max: usize) -> usize
    where
        T: Clone + Zero + One + PartialEq,
        for <'a> Matrix<T>: MulAssign<&'a Matrix<T>>
{
    let (nrows, ncols) = matrix.shape();
    assert_eq!(nrows, ncols);

    let identity = Matrix::identity(nrows);
    let mut a = identity.clone();

    for i in 1..=max {
        a *= matrix;
        if a == identity {
            return i;
        }
    }
    
    0
}


fn preserves_orientation<T>(matrix: &Matrix<T>) -> bool
    where
        T: Signed,
        Matrix<T>: LinearAlgebra<T>
{
    !matrix.determinant().is_negative()
}


fn operator_axis<T>(matrix: &Matrix<T>) -> Option<Vec<T>>
    where
        T: Clone + Signed,
        for <'a> &'a T: Neg<Output=T>,
        Matrix<T>: LinearAlgebra<T> + AddAssign<Matrix<T>> + SubAssign<Matrix<T>>
{
    let (dim, ncols) = matrix.shape();
    assert_eq!(dim, ncols);

    let r = if dim % 2 == 0 && !preserves_orientation(matrix) {
        matrix + Matrix::identity(dim)
    }
    else {
        matrix - Matrix::identity(dim)
    };

    if let Some(z) = r.null_space() {
        if z.ncols == 1 {
            return Some(normalized(&z.get_col(0)));
        }
    }

    None
}


struct OperatorDetails<T> {
    matrix: Matrix<T>,
    dimension: usize,
    axis: Option<Vec<T>>,
    order: usize,
    is_orientation_preserving: bool,
    is_clockwise: bool,
}


impl<T> OperatorDetails<T>
    where
        T: Clone + Zero + One + Signed + PartialEq,
        for <'a> &'a T: Neg<Output=T>,
        Matrix<T>: LinearAlgebra<T> + AddAssign<Matrix<T>> + SubAssign<Matrix<T>>,
        for <'a> Matrix<T>: MulAssign<&'a Matrix<T>>
{
    fn from(matrix: Matrix<T>) -> Self {
        let (nrows, ncols) = matrix.shape();
        assert_eq!(nrows, ncols);

        let dimension = nrows;
        let is_orientation_preserving = preserves_orientation(&matrix);
        let order = matrix_order(&matrix, 6);
        let axis = operator_axis(&matrix);

        todo!()
    }
}
