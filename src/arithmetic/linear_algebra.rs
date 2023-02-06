use std::ops::{Add, SubAssign, Neg, Mul, Div};

use num_rational::{BigRational};
use num_traits::{Zero, One, Inv};

use super::matrices::Matrix;

pub trait Field: Sized + Add + Zero + Neg + Mul + One + Inv {}

impl Field for BigRational {}
// tag other types as fields, e.g. prime residual classes


pub fn extend_basis<T>(v: &[T], bs: &mut Vec<Matrix<T>>)
    where
        T: Field + Clone,
        for <'a> &'a T: Div<&'a T, Output=T>,
        Matrix<T>: Eq + Neg<Output=Matrix<T>> + SubAssign,
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

    if !pivot_column(&v).is_none() {
        bs.push(v);
    }
}
