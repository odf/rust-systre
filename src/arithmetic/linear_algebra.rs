use std::ops::{Add, Sub, Mul, Div, Neg, SubAssign};
use std::ops::{AddAssign, MulAssign};

use num_rational::{BigRational};
use num_traits::{Zero, One};

use super::matrices::Matrix;


pub fn gcdx<T>(a: T, b: T) -> (T, T, T, T, T) // TODO return a struct?
    where T:
        Copy + Zero + One +
        Div<Output=T> + Sub<Output=T> + Mul<Output=T>
{
    let (mut a, mut a_next) = (a, b);
    let (mut r, mut r_next) = (T::one(), T::zero());
    let (mut s, mut s_next) = (T::zero(), T::one());

    while !a_next.is_zero() {
        let q = a / a_next;
        (a, a_next) = (a_next, a - q * a_next);
        (r, r_next) = (r_next, r - q * r_next);
        (s, s_next) = (s_next, s - q * s_next);
    }

    (a, r, s, r_next, s_next)
}


pub trait Scalar: Sized + Add + Zero + Neg + Mul + One {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>);
    fn normalize_column(col: usize, v: &mut Vec<Self>);
    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>);
    fn solve_row(a: &Vec<Self>, x: &Vec<Vec<Self>>, b: &Vec<Self>)
        -> Option<Vec<Self>>;
}

impl Scalar for BigRational {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>) {
        let f = std::mem::replace(&mut v[col], Self::zero()) / &b[col];
        for k in (col + 1)..v.len() {
            v[k] -= &b[k] * &f;
        }
    }

    fn normalize_column(col: usize, v: &mut Vec<Self>) {
        let f = std::mem::replace(&mut v[col], Self::one());
        for k in (col + 1)..v.len() {
            v[k] /= &f;
        }
    }

    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>) {
        let f = std::mem::replace(&mut v[col], Self::zero());
        for k in (col + 1)..v.len() {
            let a = &b[k] * &f;
            v[k] -= a;
        }
    }

    fn solve_row(_a: &Vec<Self>, _x: &Vec<Vec<Self>>, b: &Vec<Self>)
        -> Option<Vec<Self>>
    {
        Some(b.clone())
    }
}

impl Scalar for f64 {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>) {
        let f = v[col] / b[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }

    fn normalize_column(col: usize, v: &mut Vec<Self>) {
        let f = v[col];
        v[col] = 1.0;

        for k in (col + 1)..v.len() {
            v[k] /= f;
        }
    }

    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>) {
        let f = v[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }

    fn solve_row(_a: &Vec<Self>, _x: &Vec<Vec<Self>>, b: &Vec<Self>)
        -> Option<Vec<Self>>
    {
        Some(b.clone())
    }
}

impl Scalar for i32 {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>) {
        clear_column_primitive_int(col, v, b);
    }

    fn normalize_column(col: usize, v: &mut Vec<Self>) {
        normalize_column_primitive_int(col, v);
    }

    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>) {
        reduce_column_primitive_int(col, v, b);
    }

    fn solve_row(a: &Vec<Self>, x: &Vec<Vec<Self>>, b: &Vec<Self>)
        -> Option<Vec<Self>>
    {
        solve_row_primitive_int(a, x, b)
    }
}

impl Scalar for i64 {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>) {
        clear_column_primitive_int(col, v, b);
    }

    fn normalize_column(col: usize, v: &mut Vec<Self>) {
        normalize_column_primitive_int(col, v);
    }

    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>) {
        reduce_column_primitive_int(col, v, b);
    }

    fn solve_row(a: &Vec<Self>, x: &Vec<Vec<Self>>, b: &Vec<Self>)
        -> Option<Vec<Self>>
    {
        solve_row_primitive_int(a, x, b)
    }
}


fn clear_column_primitive_int<T>(col: usize, v: &mut Vec<T>, b: &mut Vec<T>)
    where T: num_traits::PrimInt
{
    let (_, r, s, t, u) = gcdx(b[col], v[col]);
    let det = r * u - s * t;

    for k in col..v.len() {
        let tmp = det * (b[k] * r + v[k] * s);
        v[k] = b[k] * t + v[k] * u;
        b[k] = tmp;
    }
}

fn normalize_column_primitive_int<T>(col: usize, v: &mut Vec<T>)
    where T: num_traits::PrimInt + Neg<Output=T>
{    if v[col] < T::zero() {
        for k in col..v.len() {
            v[k] = -v[k];
        }
    }
}

fn reduce_column_primitive_int<T>(col: usize, v: &mut Vec<T>, b: &Vec<T>)
    where T: num_traits::PrimInt + SubAssign
{
    let f = v[col] / b[col] - (
        if v[col] < T::zero() { T::one() } else { T::zero() }
    );

    if f != T::zero() {
        for k in col..v.len() {
            v[k] -= b[k] * f;
        }
    }
}

fn solve_row_primitive_int<T>(a: &Vec<T>, x: &Vec<Vec<T>>, b: &Vec<T>)
    -> Option<Vec<T>>
    where T: num_traits::PrimInt + SubAssign
{
    let k = x.len();
    let mut result = vec![T::zero(); b.len()];

    for col in 0..b.len() {
        let mut t = b[col];
        for i in 0..k {
            t -= a[i] * x[i][col];
        }
        if t % a[k] == T::zero() {
            result[col] = t / a[k];
        } else {
            return None;
        }
    }
    Some(result)
}


pub fn extend_basis<T>(v: &[T], bs: &mut Vec<Vec<T>>)
    where T: Scalar + Clone + Neg<Output=T>,
{
    let pivot_column = |v: &Vec<T>| { v.iter().position(|x| !x.is_zero()) };

    let mut v = v.to_vec();

    for i in 0..bs.len() {
        let mut b = &mut bs[i];
        assert!(b.len() == v.len());

        if let Some(col) = pivot_column(&v) {
            let col_b = pivot_column(b).unwrap();

            if col < col_b {
                if (bs.len() - i) % 2 > 0 {
                    v = v.iter().map(|x| -x.clone()).collect();
                }
                bs.insert(i, v);
                return;
            } else if col == col_b {
                Scalar::clear_column(col, &mut v, &mut b);
            }
        } else {
            break;
        }
    }

    if pivot_column(&v).is_some() {
        bs.push(v);
    }
}


pub trait LinearAlgebra<T> where Self: Sized {
    fn rank(&self) -> usize;
    fn determinant(&self) -> T;
    fn reduced_basis(&self) -> Self;
    fn solve(&self, rgt: &Self) -> Option<Self>;
    fn null_space(&self) -> Option<Self>;
    fn inverse(&self) -> Option<Self>;
}


impl<T> LinearAlgebra<T> for Matrix<T>
    where
        T: Scalar + Clone + Sub<Output=T> + Neg<Output=T>,
        for <'a> T: Mul<&'a T, Output=T> + AddAssign<&'a T> + MulAssign<&'a T>,
        for <'a> &'a T: Mul<&'a T, Output=T>,
{
    fn rank(&self) -> usize {
        let (nrows, _) = self.shape();

        let mut basis = vec![];
        for i in 0..nrows {
            extend_basis(&self.get_row(i), &mut basis);
        }
        basis.len()
    }

    fn determinant(&self) -> T {
        let (nrows, ncols) = self.shape();
        assert!(nrows == ncols, "must be a square matrix");

        match self.nrows {
            0 => T::one(),
            1 => self[(0, 0)].clone(),
            2 => {
                &self[(0, 0)] * &self[(1, 1)] - &self[(0, 1)] * &self[(1, 0)]
            },
            3 => {
                &self[(0, 0)] * &self[(1, 1)] * &self[(2, 2)] +
                &self[(0, 1)] * &self[(1, 2)] * &self[(2, 0)] +
                &self[(0, 2)] * &self[(1, 0)] * &self[(2, 1)] -
                &self[(0, 2)] * &self[(1, 1)] * &self[(2, 0)] -
                &self[(0, 0)] * &self[(1, 2)] * &self[(2, 1)] -
                &self[(0, 1)] * &self[(1, 0)] * &self[(2, 2)]
            },
            _ => {
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
    }

    fn reduced_basis(&self) -> Matrix<T> {
        let (nrows, _) = self.shape();

        let mut basis = vec![];
        for i in 0..nrows {
            extend_basis(&self.get_row(i), &mut basis);
        }

        let mut col = 0;
        for row in 0..basis.len() {
            while basis[row][col].is_zero() {
                col += 1;
            }

            Scalar::normalize_column(col, &mut basis[row]);

            let b = basis[row].clone();
            for i in 0..row {
                Scalar::reduce_column(col, &mut basis[i], &b);
            }
        }

        Matrix::vstack(
            &basis.iter().map(|v| Matrix::row(&v)).collect::<Vec<_>>()
        )
    }

    fn solve(&self, rgt: &Matrix<T>) -> Option<Matrix<T>> {
        let (rowslft, colslft) = self.shape();
        let (rowsrgt, colsrgt) = rgt.shape();
        assert!(rowslft == rowsrgt);

        let t = Matrix::hstack(&[self.clone(), rgt.clone()]).reduced_basis();
        if t.nrows == 0 {
            return Some(Matrix::zero(colslft, colsrgt));
        }

        let lft = t.submatrix(0..t.nrows, 0..colslft);
        let rgt = t.submatrix(0..t.nrows, colslft..t.ncols);
        if lft.rank() < lft.nrows {
            return None;
        }

        let t = Matrix::hstack(&[lft.transpose(), Matrix::identity(lft.ncols)])
            .reduced_basis();
        let b = t.submatrix(0..t.nrows, 0..lft.nrows).transpose();
        let u = t.submatrix(0..t.nrows, lft.nrows..t.ncols).transpose();

        let mut x = vec![];
        for i in 0..rgt.nrows {
            match Scalar::solve_row(&b.get_row(i), &x, &rgt.get_row(i)) {
                Some(r) => x.push(r),
                None => return None,
            }
        }
        for _i in rgt.nrows..u.ncols {
            x.push(vec![T::zero(); rgt.ncols]);
        }

        let x = Matrix::vstack(
            &x.iter().map(|v| Matrix::row(&v)).collect::<Vec<_>>()
        );
        Some(u * x)
    }

    fn null_space(&self) -> Option<Matrix<T>> {
        let m = self.transpose();
        let t = Matrix::hstack(&[m.clone(), Matrix::identity(m.nrows)])
            .reduced_basis();
        let lft = t.submatrix(0..t.nrows, 0..m.ncols);
        let rgt = t.submatrix(0..t.nrows, m.ncols..t.ncols);
        let k = lft.rank();

        if k < rgt.nrows {
            Some(rgt.submatrix(k..rgt.nrows, 0..rgt.ncols).transpose())
        } else {
            None
        }
    }

    fn inverse(&self) -> Option<Matrix<T>> {
        let (nrows, ncols) = self.shape();
        assert_eq!(nrows, ncols);
        self.solve(&Matrix::identity(nrows))
    }
}


mod test {
    use super::*;
    use num_bigint::BigInt;

    fn _r(m: &Matrix<i64>) -> Matrix<BigRational> {
        m.iter()
            .map(|&n| BigRational::from(BigInt::from(n)))
            .collect::<Matrix<_>>()
            .reshape(m.ncols)
    }

    #[test]
    fn test_extend_basis_i64() {
        let mut basis = vec![];

        extend_basis(&vec![5, 8, 13], &mut basis);
        assert_eq!(basis, vec![vec![5, 8, 13]]);

        extend_basis(&vec![7, 12, 19], &mut basis);
        assert_eq!(basis, vec![vec![1, 0, 1], vec![0, 4, 4]]);

        extend_basis(&vec![3, 5, 9], &mut basis);
        assert_eq!(basis, vec![vec![1, 0, 1], vec![0, -1, -2], vec![0, 0, -4]]);
    }

    #[test]
    fn test_reduced_basis() {
        assert_eq!(
            _r(&Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 8])).reduced_basis(),
            _r(&Matrix::new(3, &[1, 0, 0, 0, 1, 0, 0, 0, 1]))
        );

        assert_eq!(
            _r(&Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 9])).reduced_basis(),
            _r(&Matrix::new(3, &[1, 0, -1, 0, 1, 2]))
        );

        assert_eq!(
            Matrix::new(3, &[1, 4, 7, 2, 5, 8, 3, 6, 8]).reduced_basis(),
            Matrix::new(3, &[1, 1, 0, 0, 3, 0, 0, 0, 1])
        );
    }

    #[test]
    fn test_determinant() {
        for d in 0..=4 {
            let t = (2 as i32).pow(d as u32);

            let a = Matrix::<f64>::identity(d);
            assert_eq!(a.determinant(), 1.0);
            assert_eq!((&a + &a).determinant(), t as f64);

            let a = Matrix::<BigRational>::identity(d);
            assert_eq!(a.determinant(), BigRational::one());
            assert_eq!(
                (&a + &a).determinant(),
                BigRational::from(BigInt::from(t))
            );
        }

        assert_eq!(
            Matrix::new(3, &[1.0, 0.3, 0.7, 0.0, 2.0, 1.2, 0.0, 0.0, 0.25])
                .determinant(),
            0.5
        );
    }

    #[test]
    fn test_determinant_i64() {
        let a = Matrix::new(3, &[
            5, 8, 13,
            7, 12, 19,
            3, 5, 9
        ]);
        assert_eq!(a.determinant(), 4);

        let a = Matrix::new(3, &[
            5, 8, 13,
            7, 12, 19,
            3, 5, 8
        ]);
        assert_eq!(a.determinant(), 0);
    }

    #[test]
    fn test_matrix_rank() {
        let a = Matrix::new(3, &[
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        ]);
        assert_eq!(a.rank(), 3);

        let a = Matrix::new(3, &[
            1.0, 0.3, 0.7,
            0.0, 2.0, 1.2,
            0.0, 0.0, 0.25,
        ]);
        assert_eq!(a.rank(), 3);

        let a = _r(&Matrix::new(3, &[1, 0, 0, 0, 1, 0, 0, 0, 1]));
        assert_eq!(a.rank(), 3);

        let a = _r(&Matrix::new(3, &[1, 2, 3, 2, 4, 6, 3, 6, 9]));
        assert_eq!(a.rank(), 1);
        assert_eq!(a.determinant(), BigRational::zero());

        let a = _r(&Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 9]));
        assert_eq!(a.rank(), 2);
        assert_eq!(a.determinant(), BigRational::zero());

        let a = _r(&Matrix::new(3, &[7, 8, 8, 4, 5, 6, 1, 2, 3]));
        assert_eq!(a.rank(), 3);
        assert_eq!(a.determinant(), BigRational::from(BigInt::from(-3)));
    }

    #[test]
    fn test_linalg_solve() {
        let a = _r(&Matrix::new(2, &[1, 2, 3, 4]));
        let x = _r(&Matrix::new(2, &[1, 0, 1, -3]));
        let b = &a * &x;
        assert_eq!(b, _r(&Matrix::new(2, &[3, -6, 7, -12])));
        let s = a.solve(&b).unwrap();
        assert_eq!(s, x);

        let id = Matrix::identity(2);
        let a_inv = a.solve(&id).unwrap();
        assert_eq!(a * a_inv, id);

        let a = _r(&Matrix::new(2, &[1, 2, 3, 6]));
        let x = _r(&Matrix::new(1, &[1, 1]));
        let b = &a * &x;
        assert_eq!(b, _r(&Matrix::new(1, &[3, 9])));
        let s = a.solve(&b).unwrap();
        assert_eq!(&a * s, b);

        let a = Matrix::new(2, &[1.0, 2.0, 3.0, 4.0]);
        let x = Matrix::new(2, &[1.0, 0.0, 1.0, -3.0]);
        let b = &a * &x;
        assert_eq!(b, Matrix::new(2, &[3.0, -6.0, 7.0, -12.0]));
        let s = a.solve(&b).unwrap();
        assert_eq!(s, x);

        let id = Matrix::identity(2);
        let a_inv = a.solve(&id).unwrap();
        assert_eq!(a * a_inv, id);
    }

    #[test]
    fn test_linalg_nullspace() {
        let a = _r(&Matrix::new(2, &[1, 2, 3, 4]));
        let n = a.null_space();
        assert_eq!(n, None);

        let a = _r(&Matrix::new(2, &[1, 2, 3, 6]));
        let n = a.null_space().unwrap();
        assert_eq!(&a * &n, Matrix::zero(2, 1));

        let a = Matrix::new(2, &[1.0, 2.0, 3.0, 4.0]);
        let n = a.null_space();
        assert_eq!(n, None);

        let a = Matrix::new(3, &[0.0, 1.0, 0.0]);
        let n = a.null_space().unwrap();
        assert_eq!(&n, &Matrix::new(2, &[1.0, 0.0, 0.0, 0.0, 0.0, 1.0]));
    }

    #[test]
    fn test_inverse() {
        let a = _r(&Matrix::new(2, &[1, 2, 3, 4]));
        let b = a.inverse().unwrap();
        assert_eq!(a * b, Matrix::identity(2));

        let a = _r(&Matrix::new(2, &[1, 2, 3, 6]));
        assert_eq!(a.inverse(), None);
    }

    #[test]
    fn test_inverse_i64() {
        let a = &Matrix::new(2, &[1, 2, 3, 5]);
        let b = a.inverse().unwrap();
        assert_eq!(a * b, Matrix::identity(2));

        let a = &Matrix::new(2, &[1, 2, 3, 4]);
        assert_eq!(a.inverse(), None);
    }
}
