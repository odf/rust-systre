use std::ops::{Add, SubAssign, Neg, Mul, Div, MulAssign};

use num_rational::{BigRational};
use num_traits::{Zero, One, Inv};

use super::matrices::Matrix;

pub trait Field: Sized + Add + Zero + Neg + Mul + One + Inv {}

impl Field for BigRational {}
impl Field for f64 {}
// tag other types as fields, e.g. prime residual classes


pub fn extend_basis<T>(v: &[T], bs: &mut Vec<Matrix<T>>)
    where
        T: Field + Clone,
        for <'a> &'a T: Div<&'a T, Output=T>,
        Matrix<T>: Neg<Output=Matrix<T>> + SubAssign,
        for <'a> &'a Matrix<T>: Mul<T, Output=Matrix<T>>
{
    let pivot_column = |m: &Matrix<T>| {
        (0..m.ncols).position(|i| !m[(0, i)].is_zero())
    };

    let mut v = Matrix::row(&v);

    for i in 0..bs.len() {
        let b = &bs[i];
        assert!(b.shape() == v.shape());

        if let Some(col) = pivot_column(&v) {
            let col_b = pivot_column(b).unwrap();

            if col < col_b {
                bs.insert(i, if (bs.len() - i) % 2 > 0 { -v } else { v });
                return;
            } else if col == col_b {
                v -= b * (&v[(0, col)] / &b[(0, col)]);
            }
        } else {
            break;
        }
    }

    if pivot_column(&v).is_some() {
        bs.push(v);
    }
}


impl<T> Matrix<T>
    where
        T: Field + Clone + for <'a> MulAssign<&'a T>,
        for <'a> &'a T: Div<&'a T, Output=T>,
        Matrix<T>: Neg<Output=Matrix<T>> + SubAssign,
        for <'a> &'a Matrix<T>: Mul<T, Output=Matrix<T>>
{
    pub fn rank(&self) -> usize {
        let (nrows, _) = self.shape();

        let mut basis = vec![];
        for i in 0..nrows {
            extend_basis(&self.get_row(i), &mut basis);
        }
        basis.len()
    }

    pub fn determinant(&self) -> T {
        let (nrows, ncols) = self.shape();
        assert!(nrows == ncols, "must be a square matrix");

        let mut basis = vec![];
        for i in 0..nrows {
            extend_basis(&self.get_row(i), &mut basis);
        }

        if basis.len() < nrows {
            T::zero()
        } else {
            let mut det = T::one();
            for i in 0..nrows {
                det *= &basis[i][(0, i)];
            }
            det
        }
    }
}


#[test]
fn test_matrix_rank_and_det() {
    let r = |n: i64| BigRational::from(num_bigint::BigInt::from(n));

    let a = Matrix::new(3, &[
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    ]);
    assert_eq!(a.rank(), 3);
    assert_eq!(a.determinant(), 1.0);

    let a = Matrix::new(3, &[
        1.0, 0.3, 0.7,
        0.0, 2.0, 1.2,
        0.0, 0.0, 0.25,
    ]);
    assert_eq!(a.rank(), 3);
    assert_eq!(a.determinant(), 0.5);

    let a = Matrix::new(3, &[
        r(1), r(0), r(0),
        r(0), r(1), r(0),
        r(0), r(0), r(1),
    ]);
    assert_eq!(a.rank(), 3);
    assert_eq!(a.determinant(), r(1));

    let a = Matrix::new(3, &[
        r(1), r(2), r(3),
        r(2), r(4), r(6),
        r(3), r(6), r(9),
    ]);
    assert_eq!(a.rank(), 1);
    assert_eq!(a.determinant(), r(0));

    let a = Matrix::new(3, &[
        r(1), r(2), r(3),
        r(4), r(5), r(6),
        r(7), r(8), r(9),
    ]);
    assert_eq!(a.rank(), 2);
    assert_eq!(a.determinant(), r(0));

    let a = Matrix::new(3, &[
        r(1), r(2), r(3),
        r(4), r(5), r(6),
        r(7), r(8), r(8),
    ]);
    assert_eq!(a.rank(), 3);
    assert_eq!(a.determinant(), r(3));
}
