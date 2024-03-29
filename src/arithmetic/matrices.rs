use std::fmt::Display;
use std::ops::{Add, AddAssign, Sub, SubAssign, Neg};
use std::ops::{Div, DivAssign, Mul, MulAssign};
use std::ops::{Index, IndexMut, Range};
use std::slice::Iter;
use std::vec::IntoIter;
use num_traits::{One, Zero};


#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct Matrix<T> {
    pub(crate) nrows: usize,
    pub(crate) ncols: usize,
    pub(crate) data: Vec<T>,
}

impl<T> Matrix<T> {
    pub fn shape(&self) -> (usize, usize) {
        (self.nrows, self.ncols)
    }

    pub fn iter(&self) -> Iter<T> {
        self.data.iter()
    }

    pub fn into_iter(self) -> IntoIter<T> {
        self.data.into_iter()
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &T {
        assert!(i * self.ncols < self.data.len());
        assert!(j < self.ncols);

        &self.data[i * self.ncols + j]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        assert!(i * self.ncols < self.data.len());
        assert!(j < self.ncols);

        &mut self.data[i * self.ncols + j]
    }
}

impl<T: Display> Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (nrows, ncols) = self.shape();

        for i in 0..nrows {
            for j in 0..ncols {
                write!(f, " {}", self[(i, j)])?;
            }
            write!(f, "\n")?;
        }

        Ok(())
    }
}

impl<T: Clone> FromIterator<T> for Matrix<T> {
    fn from_iter<I: IntoIterator<Item=T>>(iter: I) -> Self {
        Matrix::new(1, &Vec::from_iter(iter))
    }
}

impl<T: Clone> Matrix<T> {
    pub fn new(ncols: usize, data: &[T]) -> Matrix<T> {
        let nrows = if data.len() == 0 { 0 } else { data.len() / ncols };
        assert_eq!(data.len(), nrows * ncols);

        Matrix { nrows, ncols, data: data.to_vec() }
    }

    pub fn row(data: &[T]) -> Matrix<T> {
        Self::new(data.len(), data)
    }

    pub fn col(data: &[T]) -> Matrix<T> {
        Self::new(1, data)
    }

    pub fn reshape(self, ncols: usize) -> Matrix<T> {
        Matrix::new(ncols, &self.data)
    }

    pub fn hstack(parts: &[Matrix<T>]) -> Matrix<T> {
        if parts.len() == 0 {
            Matrix::new(0, &[])
        } else {
            let nrows = parts[0].nrows;
            assert!(parts.iter().all(|m| m.nrows == nrows));
            let ncols: usize = parts.iter().map(|m| m.ncols).sum();
            let mut data = Vec::with_capacity(nrows * ncols);

            for irow in 0..nrows {
                for p in parts {
                    for icol in 0..p.ncols {
                        data.push(p[(irow, icol)].clone())
                    }
                }
            }
            Matrix::new(ncols, &data)
        }
    }

    pub fn vstack(parts: &[Matrix<T>]) -> Matrix<T> {
        if parts.len() == 0 {
            Matrix::new(0, &[])
        } else {
            let ncols = parts[0].ncols;
            assert!(parts.iter().all(|m| m.ncols == ncols));
            let nrows: usize = parts.iter().map(|m| m.nrows).sum();
            let mut data = Vec::with_capacity(nrows * ncols);

            for p in parts {
                for irow in 0..p.nrows {
                    for icol in 0..p.ncols {
                        data.push(p[(irow, icol)].clone())
                    }
                }
            }
            Matrix::new(ncols, &data)
        }
    }

    pub fn submatrix(&self, rows: Range<usize>, cols: Range<usize>) -> Self {
        assert!(rows.end <= self.nrows);
        assert!(cols.end <= self.ncols);

        let nrows = rows.end - rows.start;
        let ncols = cols.end - cols.start;
        let mut data = Vec::with_capacity(nrows * ncols);

        for i in 0..nrows {
            for j in 0..ncols {
                data.push(self[(i + rows.start, j + cols.start)].clone());
            }
        }

        Matrix::new(ncols, &data)
    }

    pub fn get_row(&self, i: usize) -> Vec<T> {
        let m = self.ncols;
        assert!(i * m < self.data.len());

        self.data[(i * m)..(i * m + m)].to_vec()
    }

    pub fn get_rows(&self) -> Vec<Vec<T>> {
        (0..self.nrows).map(|i| self.get_row(i)).collect()
    }

    pub fn get_col(&self, j: usize) -> Vec<T> {
        let (nrows, ncols) = self.shape();
        assert!(j < ncols);

        (0..nrows).map(|i| self[(i, j)].clone()).collect()
    }

    pub fn transpose(&self) -> Matrix<T> {
        let (nrows, ncols) = self.shape();
        let data: Vec<_> = (0..ncols)
            .flat_map(|j| (0..nrows).map(move |i| self[(i, j)].clone()))
            .collect();

        Matrix::new(nrows, &data)
    }
}

impl<T: Clone + Zero> Matrix<T> {
    pub fn zero(nrows: usize, ncols: usize) -> Matrix<T> {
        Self::new(ncols, &vec![T::zero(); nrows * ncols])
    }
}

impl<T: Clone + Zero + One> Matrix<T> {
    pub fn unit_vector(n: usize, k: usize) -> Matrix<T> {
        assert!(k < n);
        let mut row = Matrix::zero(1, n);
        row[(0, k)] = T::one();
        row
    }

    pub fn identity(n: usize) -> Matrix<T> {
        let mut id = Matrix::zero(n, n);
        for i in 0..n {
            id[(i, i)] = T::one();
        }
        id
    }

    pub fn is_identity(&self) -> bool
        where T: PartialEq
    {
        let (nrows, ncols) = self.shape();
        if nrows != ncols {
            return false
        }
        for i in 0..nrows {
            for j in 0..ncols {
                if j == i {
                    if !self[(i, j)].is_one() {
                        return false;
                    }
                } else {
                    if !self[(i, j)].is_zero() {
                        return false;
                    }
                }
            }
        }
        true
    }
}


#[test]
pub(crate) fn test_matrix_basics() {
    let a = Matrix::new(3, &[1, 2, 3, 4, 5, 6]);

    assert_eq!(a.shape(), (2, 3));
    assert_eq!(a[(0, 1)], 2);
    assert_eq!(a[(1, 2)], 6);

    assert_eq!(a.get_row(1), vec![4, 5, 6]);
    assert_eq!(a.get_col(1), vec![2, 5]);
    assert_eq!(a.transpose(), Matrix::new(2, &[1, 4, 2, 5, 3, 6]));

    assert_eq!(
        a.iter().cloned().collect::<Vec<_>>(),
        vec![1, 2, 3, 4, 5, 6]
    );
    assert_eq!(a.iter().cloned().collect::<Matrix<_>>().reshape(3), a);
    assert_eq!(
        a.iter().map(|i| i * i).collect::<Matrix<_>>().reshape(3),
        Matrix::new(3, &[1, 4, 9, 16, 25, 36]),
    );

    let mut a = a;
    a[(1, 1)] = -5;
    assert_eq!(a, Matrix::new(3, &[1, 2, 3, 4, -5, 6]));

    assert_eq!(
        Matrix::unit_vector(5, 2),
        Matrix::new(5, &[0, 0, 1, 0, 0])
    );
    assert_eq!(
        Matrix::identity(3),
        Matrix::new(3, &[1, 0, 0, 0, 1, 0, 0, 0, 1])
    );
    assert!(Matrix::<f32>::identity(6).is_identity());
    assert!(!Matrix::new(5, &[0, 0, 1, 0, 0]).is_identity());
    assert!(!Matrix::new(3, &[1, 0, 0, 0, 2, 0, 0, 0, 1]).is_identity());
    assert!(!Matrix::new(3, &[1, 0, 1, 0, 1, 0, 0, 0, 1]).is_identity());

    let t: Matrix<i64> = Matrix::new(3, &[]);
    assert_eq!(t.shape(), (0, 3));

    assert_eq!(
        Matrix::hstack(&[
            Matrix::col(&[1, 2]),
            Matrix::new(2, &[3, 4, 5, 6]),
            Matrix::col(&[7, 8]),
        ]),
        Matrix::new(4, &[1, 3, 4, 7, 2, 5, 6, 8])
    );

    assert_eq!(
        Matrix::vstack(&[
            Matrix::row(&[1, 2]),
            Matrix::new(2, &[3, 4, 5, 6]),
            Matrix::row(&[7, 8]),
        ]),
        Matrix::new(2, &[1, 2, 3, 4, 5, 6, 7, 8])
    );

    assert_eq!(
        Matrix::new(3, &[1, 2, 3, 4, 5, 6, 7, 8, 9]).submatrix(1..3, 0..2),
        Matrix::new(2, &[4, 5, 7, 8]),
    )
}


impl<T> Neg for &Matrix<T>
    where T: Clone + Neg<Output=T>
{
    type Output = Matrix<T>;

    fn neg(self) -> Self::Output {
        let data: Vec<_> = self.data.iter().map(|x| -x.clone()).collect();
        Matrix::new(self.ncols, &data)
    }
}

impl<T> Neg for Matrix<T>
    where T: Clone + Neg<Output=T>
{
    type Output = Matrix<T>;

    fn neg(self) -> Self::Output {
        -(&self)
    }
}


#[test]
pub(crate) fn test_matrix_neg() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!(-&a, Matrix::new(2, &[-1, -2, -3, -4]));
    assert_eq!(-a, Matrix::new(2, &[-1, -2, -3, -4]));
}


impl<T> AddAssign<&Matrix<T>> for Matrix<T>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] += &rhs.data[i];
        }
    }
}

impl<T> AddAssign<Matrix<T>> for Matrix<T>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] += &rhs.data[i];
        }
    }
}

impl<T> AddAssign<&T> for Matrix<T>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &T) {
        for i in 0..self.data.len() {
            self.data[i] += &rhs;
        }
    }
}

impl<T> AddAssign<T> for Matrix<T>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: T) {
        for i in 0..self.data.len() {
            self.data[i] += &rhs;
        }
    }
}

impl<S, T> Add<S> for &Matrix<T>
    where T: Clone, Matrix<T>: AddAssign<S>
{
    type Output = Matrix<T>;

    fn add(self, rhs: S) -> Matrix<T> {
        let mut tmp = (*self).clone();
        tmp += rhs;
        tmp
    }
}

impl<S, T> Add<S> for Matrix<T>
    where T: Clone, Matrix<T>: AddAssign<S>
{
    type Output = Matrix<T>;

    fn add(self, rhs: S) -> Matrix<T> {
        &self + rhs
    }
}

#[test]
pub(crate) fn test_matrix_add() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!(&a + 2, Matrix::new(2, &[3, 4, 5, 6]));
    assert_eq!(&a + &5, Matrix::new(2, &[6, 7, 8, 9]));
    assert_eq!(a.clone() + &a, Matrix::new(2, &[2, 4, 6, 8]));
    assert_eq!(&a + a.clone() + &1, Matrix::new(2, &[3, 5, 7, 9]));
    assert_eq!(a.clone() + a.clone() + -2, Matrix::new(2, &[0, 2, 4, 6]));
}

impl<T> SubAssign<&Matrix<T>> for Matrix<T>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] -= &rhs.data[i];
        }
    }
}

impl<T> SubAssign<Matrix<T>> for Matrix<T>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] -= &rhs.data[i];
        }
    }
}

impl<T> SubAssign<&T> for Matrix<T>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &T) {
        for i in 0..self.data.len() {
            self.data[i] -= &rhs;
        }
    }
}

impl<T> SubAssign<T> for Matrix<T>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..self.data.len() {
            self.data[i] -= &rhs;
        }
    }
}

impl<S, T> Sub<S> for &Matrix<T>
    where T: Clone, Matrix<T>: SubAssign<S>
{
    type Output = Matrix<T>;

    fn sub(self, rhs: S) -> Matrix<T> {
        let mut tmp = (*self).clone();
        tmp -= rhs;
        tmp
    }
}

impl<S, T> Sub<S> for Matrix<T>
    where T: Clone, Matrix<T>: SubAssign<S>
{
    type Output = Matrix<T>;

    fn sub(self, rhs: S) -> Matrix<T> {
        &self - rhs
    }
}

#[test]
pub(crate) fn test_matrix_sub() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!(&a - 2, Matrix::new(2, &[-1, 0, 1, 2]));
    assert_eq!(&a - &5, Matrix::new(2, &[-4, -3, -2, -1]));
    assert_eq!(a.clone() - &a, Matrix::zero(2, 2));
    assert_eq!(&a - a.clone() - &1, Matrix::new(2, &[-1, -1, -1, -1]));
    assert_eq!(a.clone() - a.clone() - -2, Matrix::new(2, &[2, 2, 2, 2]));
}

impl<T> MulAssign<&T> for Matrix<T>
    where T: for <'a> MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &T) {
        for i in 0..self.data.len() {
            self.data[i] *= &rhs;
        }
    }
}

impl<T> MulAssign<T> for Matrix<T>
    where T: for <'a> MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.data.len() {
            self.data[i] *= &rhs;
        }
    }
}

impl<S, T> Mul<S> for &Matrix<T>
    where T: Clone, Matrix<T>: MulAssign<S>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: S) -> Matrix<T> {
        let mut tmp = (*self).clone();
        tmp *= rhs;
        tmp
    }
}

impl<S, T> Mul<S> for Matrix<T>
    where T: Clone, Matrix<T>: MulAssign<S>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: S) -> Matrix<T> {
        &self * rhs
    }
}

impl<T> Mul<&Matrix<T>> for &Matrix<T>
    where T:
        Clone + Zero +
        for <'a> MulAssign<&'a T> +
        for <'a> AddAssign<&'a T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: &Matrix<T>) -> Matrix<T> {
        let (rowslft, colslft) = self.shape();
        let (rowsrgt, colsrgt) = rhs.shape();
        assert!(colslft == rowsrgt);

        let mut result: Matrix<T> = Matrix::zero(rowslft, colsrgt);

        for i in 0..rowslft {
            for j in 0..colsrgt {
                let mut s = T::zero();
                for k in 0..colslft {
                    let mut t = self[(i, k)].clone();
                    t *= &rhs[(k, j)];
                    s += &t;
                }
                result[(i, j)] = s;
            }
        }

        result
    }
}

impl<T> Mul<&Matrix<T>> for Matrix<T>
    where T:
        Clone + Zero +
        for <'a> MulAssign<&'a T> +
        for <'a> AddAssign<&'a T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: &Matrix<T>) -> Matrix<T> {
        &self * rhs
    }
}

impl<T> Mul<Matrix<T>> for &Matrix<T>
    where T:
        Clone + Zero +
        for <'a> MulAssign<&'a T> +
        for <'a> AddAssign<&'a T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: Matrix<T>) -> Matrix<T> {
        self * &rhs
    }
}

impl<T> Mul<Matrix<T>> for Matrix<T>
    where T:
        Clone + Zero +
        for <'a> MulAssign<&'a T> +
        for <'a> AddAssign<&'a T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: Matrix<T>) -> Matrix<T> {
        &self * &rhs
    }
}

#[test]
pub(crate) fn test_matrix_mul() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!(&a * 2 * &3, Matrix::new(2, &[6, 12, 18, 24]));
    assert_eq!(a.clone() * &3 * 5, Matrix::new(2, &[15, 30, 45, 60]));

    assert_eq!(&a * &a, Matrix::new(2, &[7, 10, 15, 22]));
    assert_eq!(&a * a.clone(), Matrix::new(2, &[7, 10, 15, 22]));
    assert_eq!(a.clone() * &a, Matrix::new(2, &[7, 10, 15, 22]));
    assert_eq!(a.clone() * a.clone(), Matrix::new(2, &[7, 10, 15, 22]));

    let a = Matrix::row(&[1, 2, 3]);
    let b = Matrix::new(2, &[1, 2, 3, 4, 5, 6]);
    assert_eq!(a * b, Matrix::row(&[22, 28]));
}

impl<T> DivAssign<&T> for Matrix<T>
    where T: for <'a> DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: &T) {
        for i in 0..self.data.len() {
            self.data[i] /= &rhs;
        }
    }
}

impl<T> DivAssign<T> for Matrix<T>
    where T: for <'a> DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.data.len() {
            self.data[i] /= &rhs;
        }
    }
}

impl<S, T> Div<S> for &Matrix<T>
    where T: Clone, Matrix<T>: DivAssign<S>
{
    type Output = Matrix<T>;

    fn div(self, rhs: S) -> Matrix<T> {
        let mut tmp = (*self).clone();
        tmp /= rhs;
        tmp
    }
}

impl<S, T> Div<S> for Matrix<T>
    where T: Clone, Matrix<T>: DivAssign<S>
{
    type Output = Matrix<T>;

    fn div(self, rhs: S) -> Matrix<T> {
        &self / rhs
    }
}

#[test]
pub(crate) fn test_matrix_div() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!((&a + 3) / 2, Matrix::new(2, &[2, 2, 3, 3]));
    assert_eq!(a / &3, Matrix::new(2, &[0, 0, 1, 1]));

    assert_eq!(
        Matrix::row(&[1.75, 2.25, 0.75]) / 0.25,
        Matrix::row(&[7.0, 9.0, 3.0])
    );
}
