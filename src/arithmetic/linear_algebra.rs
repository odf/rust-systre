use std::ops::{Add, SubAssign, Neg, Mul, MulAssign, Div};

use num_rational::{BigRational};
use num_traits::{Zero, One, Inv};

use super::matrices::Matrix;

pub trait Field: Sized + Add + Zero + Neg + Mul + One + Inv {}

impl Field for BigRational {}
// tag other types as fields, e.g. prime residual classes


pub fn extend_basis<T>(v: &[T], bs: &mut Vec<Matrix<T>>)
    where T: Clone + Eq + Field + Neg<Output=T> + Div<Output=T> +
        for <'a> MulAssign<&'a T> +
        for <'a> SubAssign<&'a T>
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
                if (bs.len() - i) % 2 > 0 {
                    v = v * -T::one();
                }
                bs.insert(i, v);
                return;
            } else if col == col_b {
                v -= b * (v[(0, col)].clone() / b[(0, col)].clone());
            }
        } else {
            break;
        }
    }

    if v != Matrix::zero(1, v.ncols) {
        bs.push(v);
    }
}
