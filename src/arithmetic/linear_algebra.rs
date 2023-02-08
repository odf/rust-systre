use std::ops::{Add, SubAssign, Neg, Mul, Div, MulAssign, DivAssign};

use num_rational::{BigRational};
use num_traits::{Zero, One, Inv};

use super::matrices::Matrix;

pub trait Field: Sized + Add + Zero + Neg + Mul + One + Inv {}

impl Field for BigRational {}
impl Field for f64 {}
// tag other types as fields, e.g. prime residual classes


pub fn extend_basis<T>(v: &[T], bs: &mut Vec<Vec<T>>)
    where
        T: Field + Clone + SubAssign + for <'a> Div<&'a T, Output=T>,
        for <'a> &'a T: Neg<Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T> ,
{
    let pivot_column = |v: &Vec<T>| { v.iter().position(|x| !x.is_zero()) };

    let mut v = v.to_vec();

    for i in 0..bs.len() {
        let b = &bs[i];
        assert!(b.len() == v.len());

        if let Some(col) = pivot_column(&v) {
            let col_b = pivot_column(b).unwrap();

            if col < col_b {
                if (bs.len() - i) % 2 > 0 {
                    v = v.iter().map(|x| -x).collect();
                }
                bs.insert(i, v);
                return;
            } else if col == col_b {
                let f = std::mem::replace(&mut v[col], T::zero()) / &b[col];
                for k in (col + 1)..v.len() {
                    v[k] -= &b[k] * &f;
                }
            }
        } else {
            break;
        }
    }

    if pivot_column(&v).is_some() {
        bs.push(v);
    }
}

pub fn reduced_basis<T>(m: &Matrix<T>) -> Matrix<T>
    where
        T: Field + Clone + SubAssign + for <'a> Div<&'a T, Output=T>,
        T: for <'a> DivAssign<&'a T>,
        for <'a> &'a T: Neg<Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T>,
{
    let (nrows, _) = m.shape();

    let mut basis = vec![];
    for i in 0..nrows {
        extend_basis(&m.get_row(i), &mut basis);
    }

    let mut col = 0;
    for row in 0..basis.len() {
        while basis[row][col].is_zero() {
            col += 1;
        }

        let f = std::mem::replace(&mut basis[row][col], T::one());
        for k in (col + 1)..basis[row].len() {
            basis[row][k] /= &f;
        }

        for i in 0..row {
            let f = std::mem::replace(&mut basis[i][col], T::zero());
            for k in (col + 1)..basis[i].len() {
                let a = &basis[row][k] * &f;
                basis[i][k] -= a;
            }
        }
    }

    Matrix::vstack(&basis.iter().map(|v| Matrix::row(&v)).collect::<Vec<_>>())
}


#[test]
fn test_reduced_basis() {
    fn r(m: &Matrix<i64>) -> Matrix<BigRational> {
        m.iter()
            .map(|&n| BigRational::from(num_bigint::BigInt::from(n)))
            .collect::<Matrix<_>>()
            .reshape(m.ncols)
    }

    assert_eq!(
        reduced_basis(&r(&Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 8]))),
        r(&Matrix::new(3, &[1, 0, 0, 0, 1, 0, 0, 0, 1]))
    );

    assert_eq!(
        reduced_basis(&r(&Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 9]))),
        r(&Matrix::new(3, &[1, 0, -1, 0, 1, 2]))
    );
}


pub fn solve<T>(lft: &Matrix<T>, rgt: &Matrix<T>) -> Option<Matrix<T>>
    where
        T: Field + Clone + SubAssign + for <'a> Div<&'a T, Output=T>,
        T: for <'a> DivAssign<&'a T> + for <'a> MulAssign<&'a T>,
        for <'a> &'a T: Neg<Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T>,
        Matrix<T>: Mul<Output=Matrix<T>>
{
    let (rowslft, colslft) = lft.shape();
    let (rowsrgt, colsrgt) = rgt.shape();
    assert!(rowslft == rowsrgt);

    let t = reduced_basis(&Matrix::hstack(&[lft.clone(), rgt.clone()]));
    if t.nrows == 0 {
        return Some(Matrix::zero(colslft, colsrgt));
    }

    let lft = t.submatrix(0..t.nrows, 0..colslft);
    let rgt = t.submatrix(0..t.nrows, colslft..t.ncols);
    if lft.rank() < lft.nrows {
        return None;
    }

    let t = reduced_basis(&Matrix::hstack(&[
        lft.transpose(),
        Matrix::identity(lft.ncols)
    ]));
    let u = t.submatrix(0..t.nrows, lft.nrows..t.ncols).transpose();

    let m = rgt.ncols;
    let x = Matrix::vstack(&[rgt, Matrix::zero(lft.ncols - lft.nrows, m)]);
    Some(u * x)
}


#[test]
fn test_solve() {
    fn r(m: &Matrix<i64>) -> Matrix<BigRational> {
        m.iter()
            .map(|&n| BigRational::from(num_bigint::BigInt::from(n)))
            .collect::<Matrix<_>>()
            .reshape(m.ncols)
    }

    let a = r(&Matrix::new(2, &[1, 2, 3, 4]));
    let x = r(&Matrix::new(2, &[1, 0, 1, -3]));
    let b = &a * &x;
    assert_eq!(b, r(&Matrix::new(2, &[3, -6, 7, -12])));
    let s = solve(&a, &b).unwrap();
    assert_eq!(s, x);

    let id = Matrix::identity(2);
    let a_inv = solve(&a, &id).unwrap();
    assert_eq!(a * a_inv, id);

    let a = r(&Matrix::new(2, &[1, 2, 3, 6]));
    let x = r(&Matrix::new(1, &[1, 1]));
    let b = &a * &x;
    assert_eq!(b, r(&Matrix::new(1, &[3, 9])));
    let s = solve(&a, &b).unwrap();
    assert_eq!(a * s, b);

    let a = Matrix::new(2, &[1.0, 2.0, 3.0, 4.0]);
    let x = Matrix::new(2, &[1.0, 0.0, 1.0, -3.0]);
    let b = &a * &x;
    assert_eq!(b, Matrix::new(2, &[3.0, -6.0, 7.0, -12.0]));
    let s = solve(&a, &b).unwrap();
    assert_eq!(s, x);

    let id = Matrix::identity(2);
    let a_inv = solve(&a, &id).unwrap();
    assert_eq!(a * a_inv, id);
}


impl<T> Matrix<T>
    where
        T: Field + Clone + SubAssign + for <'a> Div<&'a T, Output=T>,
        T: for <'a> DivAssign<&'a T> + for <'a> MulAssign<&'a T>,
        for <'a> &'a T: Neg<Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T>,
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
                det *= &basis[i][i];
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
        r(7), r(8), r(8),
        r(4), r(5), r(6),
        r(1), r(2), r(3),
    ]);
    assert_eq!(a.rank(), 3);
    assert_eq!(a.determinant(), r(-3));
}
