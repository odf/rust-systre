use std::ops::{SubAssign, AddAssign, Mul, Neg, MulAssign};

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
        T: PartialEq + Zero + Signed,
        for <'a> &'a T: Mul<&'a T, Output=T>
{
    dot(v, w) * dot(v, w) == dot(v, v) * dot(w, w)
}


fn matrix_order<T>(matrix: &Matrix<T>, max: usize) -> usize
    where
        T: Clone + Zero + One + PartialEq,
        for <'a> T: MulAssign<&'a T> + AddAssign<&'a T>
{
    let (nrows, ncols) = matrix.shape();
    assert_eq!(nrows, ncols);

    let identity = Matrix::identity(nrows);
    let mut a = identity.clone();

    for i in 1..=max {
        a = &a * matrix;
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
        Matrix<T>: LinearAlgebra<T>,
        Matrix<T>: AddAssign<Matrix<T>> + SubAssign<Matrix<T>>,
{
    let (dim, ncols) = matrix.shape();
    assert_eq!(dim, ncols);

    let r = if dim % 2 != 0 && !preserves_orientation(matrix) {
        matrix + Matrix::identity(dim)
    } else {
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
        for <'a> &'a T: Neg<Output=T> + Mul<&'a T, Output=T>,
        for <'a> T: MulAssign<&'a T> + AddAssign<&'a T>,
        Matrix<T>: LinearAlgebra<T>,
        Matrix<T>: AddAssign<Matrix<T>> + SubAssign<Matrix<T>>
{
    fn from(matrix: Matrix<T>) -> Self {
        let (nrows, ncols) = matrix.shape();
        assert_eq!(nrows, ncols);
        assert!(nrows == 2 || nrows == 3);

        let dimension = nrows;
        let is_orientation_preserving = preserves_orientation(&matrix);
        let order = matrix_order(&matrix, 6);
        let axis = operator_axis(&matrix);

        let is_clockwise = if dimension == 2 {
            if !is_orientation_preserving {
                false
            } else if order == 0 || order > 2 {
                !(matrix[(0, 1)]).is_negative()
            } else {
                true
            }
        } else if (order == 0 || order > 2) && axis.is_some() {
            let axis = axis.clone().unwrap();
            let u0 = [T::one(), T::zero(), T::zero()];
            let u1 = [T::zero(), T::one(), T::zero()];
            let v = if are_collinear(&axis, &u0) { u1 } else { u0 };
            let v = Matrix::col(&v);
            let w = &matrix * &v;
            let a = Matrix::hstack(&[Matrix::col(&axis), v, w]);

            preserves_orientation(&a)
        } else {
            true
        };

        OperatorDetails {
            matrix,
            dimension,
            axis,
            order,
            is_orientation_preserving,
            is_clockwise
        }
    }
}


fn crystal_system_and_basis_2d<T>(ops: &[Matrix<T>])
    -> (CrystalSystem2d, Matrix<T>)
    where
        T: Clone + Zero + One + Signed + PartialEq,
        for <'a> &'a T: Neg<Output=T> + Mul<&'a T, Output=T>,
        for <'a> T: MulAssign<&'a T> + AddAssign<&'a T>,
        Matrix<T>: LinearAlgebra<T>,
        Matrix<T>: AddAssign<Matrix<T>> + SubAssign<Matrix<T>>
{
    let ops_with_details: Vec<_> = ops.iter()
        .map(|op| OperatorDetails::from(op.clone()))
        .collect();
    let mirrors: Vec<_> = ops_with_details.iter()
        .filter(|opd| !opd.is_orientation_preserving)
        .collect();
    let n = ops_with_details.iter().map(|opd| opd.order).max().unwrap();
    let m = if n == 6 { 3 } else { n };

    let crystal_system = match n {
        3 | 6 => CrystalSystem2d::Hexagonal,
        4 => CrystalSystem2d::Square,
        _ => if mirrors.len() > 0 {
            CrystalSystem2d::Rectangular
        } else {
            CrystalSystem2d::Oblique
        }
    };

    let x = if mirrors.len() > 0 {
        Matrix::col(&mirrors[0].axis.clone().unwrap())
    } else {
        Matrix::col(&[T::one(), T::zero()])
    };

    let y = if m >= 3 {
        let s = ops_with_details.iter().find(|s|
            s.is_orientation_preserving && s.is_clockwise && s.order == m
        ).unwrap();
        &s.matrix * &x
    } else if mirrors.len() > 1 {
        Matrix::col(&mirrors[1].axis.clone().unwrap())
    } else {
        let t = if x[(0, 0)] == T::zero() {
            Matrix::col(&[T::one(), T::zero()])
        } else {
            Matrix::col(&[T::zero(), T::one()])
        };
        if mirrors.len() > 0 {
            &t - &mirrors[0].matrix * &t
        } else {
            t
        }
    };
    todo!()
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::arithmetic::matrices::Matrix;

    use super::*;

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    #[test]
    fn test_spacegroups_operator_details_2d() {
        let m = Matrix::new(2, &[r(1), r(0), r(0), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 1);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(2, &[r(1), r(0), r(0), r(-1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(0)]));
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);

        let m = Matrix::new(2, &[r(0), r(1), r(-1), r(0)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 4);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(2, &[r(0), r(-1), r(1), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 6);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, false);

        let m = Matrix::new(2, &[r(-1), r(0), r(1), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(&opd.axis, &Some(vec![r(0), r(1)]));
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);
    }

    #[test]
    fn test_spacegroups_operator_details_3d() {
        let m = Matrix::new(
            3, &[r(1), r(0), r(0), r(0), r(1), r(0), r(0), r(0), r(1)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 1);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(-1), r(0), r(0), r(0), r(-1), r(0), r(0), r(0), r(-1)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(0), r(0), r(1), r(1), r(0), r(0), r(0), r(1), r(0)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(1), r(1)]));
        assert_eq!(opd.order, 3);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(0), r(0), r(-1), r(-1), r(0), r(0), r(0), r(-1), r(0)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(1), r(1)]));
        assert_eq!(opd.order, 6);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);
    }
}
