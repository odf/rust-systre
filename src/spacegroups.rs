use num_rational::BigRational;
use num_traits::Signed;

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


fn preserves_orientation<T>(matrix: &Matrix<T>)
    -> bool
    where
        T: Signed,
        Matrix<T>: LinearAlgebra<T>
{
    !matrix.determinant().is_negative()
}
