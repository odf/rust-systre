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

impl<'a, S, T: 'a> Add<S> for &'a Matrix<T>
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
fn test_matrix_add() {
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

impl<'a, S, T: 'a> Sub<S> for &'a Matrix<T>
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
fn test_matrix_sub() {
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

impl<'a, S, T: 'a> Mul<S> for &'a Matrix<T>
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
        let rowslft = self.data.len() / self.ncols;
        let colslft = self.ncols;
        let rowsrgt = rhs.data.len() / rhs.ncols;
        let colsrgt = rhs.ncols;

        assert!(colslft == rowsrgt);

        let mut result: Matrix<T> = Matrix::zero(rowslft, colsrgt);

        for i in 0..rowslft {
            for j in 0..colsrgt {
                let mut s = T::zero();
                for k in 0..colslft {
                    let mut t = self.data[i * colslft + k].clone();
                    t *= &rhs.data[k * colsrgt + j];
                    s += &t;
                }
                result.data[i * colsrgt + j] = s;
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
fn test_matrix_mul() {
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

impl<'a, S, T: 'a> Div<S> for &'a Matrix<T>
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
fn test_matrix_div() {
    let a = Matrix::new(2, &[1, 2, 3, 4]);

    assert_eq!((&a + 3) / 2, Matrix::new(2, &[2, 2, 3, 3]));
    assert_eq!(a / &3, Matrix::new(2, &[0, 0, 1, 1]));

    assert_eq!(
        Matrix::row(&[1.75, 2.25, 0.75]) / 0.25,
        Matrix::row(&[7.0, 9.0, 3.0])
    );
}
