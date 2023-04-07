use std::ops::{SubAssign, AddAssign, Mul};

use num_rational::BigRational;
use num_traits::{Signed, Zero};

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


fn matrix_order(matrix: &Matrix<BigRational>, max: usize) -> usize {
    let (nrows, ncols) = matrix.shape();
    assert_eq!(nrows, ncols);

    let identity = Matrix::<BigRational>::identity(nrows);
    let mut a = identity.clone();

    for i in 1..=max {
        a = a * matrix;
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
            return Some(z.get_col(0));
        }
    }

    None
}
