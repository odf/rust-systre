use std::ops::{Add, Sub, Mul, Div};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use num_traits::{Zero, One};


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Matrix<T> {
    ncols: usize,
    data: Vec<T>,
}

impl<T: Clone> Matrix<T> {
    pub fn new(ncols: usize, data: &[T]) -> Matrix<T> {
        assert!(data.len() % ncols == 0);
        Matrix { ncols, data: data.to_vec() }
    }

    pub fn row(data: &[T]) -> Matrix<T> {
        Matrix { ncols: data.len(), data: data.to_vec() }
    }

    pub fn col(data: &[T]) -> Matrix<T> {
        Matrix { ncols: 1, data: data.to_vec() }
    }
}

impl<T: Clone + Zero> Matrix<T> {
    pub fn zero(nrows: usize, ncols: usize) -> Matrix<T> {
        let data = vec![T::zero(); nrows * ncols];
        Matrix { ncols, data }
    }

}

impl<T: Clone + Zero + One> Matrix<T> {
    pub fn unit_vector(n: usize, k: usize) -> Matrix<T> {
        assert!(k < n);

        let mut data = vec![T::zero(); n];
        data[k] = T::one();

        Matrix { ncols: n, data }
    }

    pub fn identity(n: usize) -> Matrix<T> {
        let mut data = vec![T::zero(); n * n];
        for i in 0..n {
            data[i * (n + 1)] = T::one();
        }

        Matrix { ncols: n, data }
    }
}

impl<S, T> AddAssign<Matrix<T>> for Matrix<S>
    where S: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] += &rhs.data[i];
        }
    }
}

impl<S, T> Add<Matrix<T>> for Matrix<S>
    where S: Clone + for <'a> AddAssign<&'a T>
{
    type Output = Matrix<S>;

    fn add(self, rhs: Matrix<T>) -> Matrix<S> {
        let mut tmp = self.clone();
        tmp += rhs;
        tmp
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

impl<T> Add<T> for Matrix<T>
    where T: Clone + for <'a> AddAssign<&'a T>
{
    type Output = Matrix<T>;

    fn add(self, rhs: T) -> Matrix<T> {
        let mut tmp = self.clone();
        tmp += rhs;
        tmp
    }
}

impl<S, T> SubAssign<Matrix<T>> for Matrix<S>
    where S: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: Matrix<T>) {
        assert!(self.data.len() == rhs.data.len());
        assert!(self.ncols == rhs.ncols);

        for i in 0..self.data.len() {
            self.data[i] -= &rhs.data[i];
        }
    }
}

impl<S, T> Sub<Matrix<T>> for Matrix<S>
    where S: Clone + for <'a> SubAssign<&'a T>
{
    type Output = Matrix<S>;

    fn sub(self, rhs: Matrix<T>) -> Matrix<S> {
        let mut tmp = self.clone();
        tmp -= rhs;
        tmp
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

impl<T> Sub<T> for Matrix<T>
    where T: Clone + for <'a> SubAssign<&'a T>
{
    type Output = Matrix<T>;

    fn sub(self, rhs: T) -> Matrix<T> {
        let mut tmp = self.clone();
        tmp -= rhs;
        tmp
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

impl<T> Mul<T> for Matrix<T>
    where T: Clone + for <'a> MulAssign<&'a T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: T) -> Matrix<T> {
        let mut tmp = self.clone();
        tmp *= rhs;
        tmp
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

impl<T> Div<T> for Matrix<T>
    where T: Clone + for <'a> DivAssign<&'a T>
{
    type Output = Matrix<T>;

    fn div(self, rhs: T) -> Matrix<T> {
        let mut tmp = self.clone();
        tmp /= rhs;
        tmp
    }
}


#[test]
fn test_matrix_ops() {
    let m = Matrix::new(2, &[1, 2, 3, 4]);

    let m = m + 2;
    assert_eq!(m, Matrix::new(2, &[3, 4, 5, 6]));

    let m = m - 3;
    assert_eq!(m, Matrix::new(2, &[0, 1, 2, 3]));

    let m = m.clone() + m;
    assert_eq!(m, Matrix::new(2, &[0, 2, 4, 6]));

    let m = m / 3;
    assert_eq!(m, Matrix::new(2, &[0, 0, 1, 2]));

    let m = m * 5;
    assert_eq!(m, Matrix::new(2, &[0, 0, 5, 10]));

    let m = m - Matrix::identity(2) * 5;
    assert_eq!(m, Matrix::new(2, &[-5, 0, 5, 5]));

    assert_eq!(
        Matrix::unit_vector(5, 2),
        Matrix::row(&[0, 0, 1, 0, 0])
    );
}
