use std::marker::PhantomData;
use std::ops::{Neg, Add, AddAssign, Sub, SubAssign};
use std::ops::{Mul, MulAssign, Div, DivAssign};
use std::ops::{Index, IndexMut};
use std::slice::Iter;

use num_traits::{Zero, One};

use super::matrices::Matrix;


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Vector<T, CS> {
    coords: Matrix<T>,
    phantom: PhantomData<CS>,
}

impl<T, CS> Vector<T, CS> {
    pub fn dim(&self) -> usize {
        self.coords.shape().0
    }

    pub fn iter(&self) -> Iter<T> {
        self.coords.iter()
    }
}

impl<T, CS> Index<usize> for Vector<T, CS> {
    type Output = T;

    fn index(&self, i: usize) -> &T {
        &self.coords[(i, 0)]
    }
}

impl<T, CS> IndexMut<usize> for Vector<T, CS> {
    fn index_mut(&mut self, i: usize) -> &mut T {
        &mut self.coords[(i, 0)]
    }
}

impl<T: Clone, CS> FromIterator<T> for Vector<T, CS> {
    fn from_iter<I: IntoIterator<Item=T>>(iter: I) -> Self {
        Vector::new(&Vec::from_iter(iter))
    }
}

impl<T: Clone, CS> From<Matrix<T>> for Vector<T, CS> {
    fn from(coords: Matrix<T>) -> Self {
        Vector { coords, phantom: PhantomData::default() }
    }
}

impl<T: Clone, CS> Vector<T, CS> {
    pub fn new(coords: &[T]) -> Vector<T, CS> {
        Vector {
            coords: Matrix::col(coords),
            phantom: PhantomData::default()
        }
    }
}

impl<T: Clone + Zero, CS> Vector<T, CS> {
    pub fn zero(dim: usize) -> Vector<T, CS> {
        Vector::new(&vec![T::zero(); dim])
    }
}

impl<T: Clone + Zero + One, CS> Vector<T, CS> {
    pub fn unit(n: usize, k: usize) -> Vector<T, CS> {
        assert!(k < n);
        let mut v = Vector::zero(n);
        v[k] = T::one();
        v
    }
}


impl<T, CS> Neg for &Vector<T, CS>
    where T: Clone + Neg<Output=T>
{
    type Output = Vector<T, CS>;

    fn neg(self) -> Self::Output {
        Vector::from(-(&self.coords))
    }
}

impl<T, CS> Neg for Vector<T, CS>
    where T: Clone + Neg<Output=T>
{
    type Output = Vector<T, CS>;

    fn neg(self) -> Self::Output {
        -(&self)
    }
}


/// # Examples
/// 
/// ```
/// use rust_systre::arithmetic::geometry::*;
/// 
/// #[derive(Debug, PartialEq)]
/// struct World {}
/// struct Local {}
/// 
/// let mut u = Vector::<_, World>::new(&[1, 2, 3]);
/// u += Vector::<_, World>::new(&[3, 2, 1]);
/// assert_eq!(u, Vector::new(&[4, 4, 4]));
/// ```
/// 
/// ```compile_fail
/// use rust_systre::arithmetic::geometry::*;
/// 
/// #[derive(Debug, PartialEq)]
/// struct World {}
/// struct Local {}
/// 
/// let mut u = Vector::<_, World>::new(&[1, 2, 3]);
/// u += Vector::<_, Local>::new(&[3, 2, 1]);
/// assert_eq!(u, Vector::new(&[4, 4, 4]));
/// ```

impl<T,  CS> AddAssign<Vector<T,  CS>> for Vector<T,  CS>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords += &rhs.coords
    }
}

impl<T,  CS> AddAssign<&Vector<T,  CS>> for Vector<T,  CS>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords += &rhs.coords
    }
}

impl<'a, S, T: 'a, CS> Add<S> for &'a Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: AddAssign<S>
{
    type Output = Vector<T, CS>;

    fn add(self, rhs: S) -> Vector<T, CS> {
        let mut tmp = (*self).clone();
        tmp += rhs;
        tmp
    }
}

impl<S, T, CS> Add<S> for Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: AddAssign<S>
{
    type Output = Vector<T, CS>;

    fn add(self, rhs: S) -> Vector<T, CS> {
        &self + rhs
    }
}

impl<T,  CS> SubAssign<Vector<T,  CS>> for Vector<T,  CS>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords -= &rhs.coords
    }
}

impl<T,  CS> SubAssign<&Vector<T,  CS>> for Vector<T,  CS>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords -= &rhs.coords
    }
}

impl<'a, S, T: 'a, CS> Sub<S> for &'a Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: SubAssign<S>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: S) -> Vector<T, CS> {
        let mut tmp = (*self).clone();
        tmp -= rhs;
        tmp
    }
}

impl<S, T, CS> Sub<S> for Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: SubAssign<S>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: S) -> Vector<T, CS> {
        &self - rhs
    }
}

impl<T, CS> MulAssign<&T> for Vector<T, CS>
    where T: for <'a> MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: &T) {
        self.coords *= rhs;
    }
}

impl<T, CS> MulAssign<T> for Vector<T, CS>
    where T: for <'a> MulAssign<&'a T>
{
    fn mul_assign(&mut self, rhs: T) {
        self.coords *= rhs;
    }
}

impl<'a, S, T: 'a, CS> Mul<S> for &'a Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: MulAssign<S>
{
    type Output = Vector<T, CS>;

    fn mul(self, rhs: S) -> Vector<T, CS> {
        let mut tmp = (*self).clone();
        tmp *= rhs;
        tmp
    }
}

impl<S, T, CS> Mul<S> for Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: MulAssign<S>
{
    type Output = Vector<T, CS>;

    fn mul(self, rhs: S) -> Vector<T, CS> {
        &self * rhs
    }
}

pub struct Point<T, CS> {
    coords: Matrix<T>,
    phantom: PhantomData<CS>,
}


pub struct ScalarProduct<T, CS> {
    coeffs: Matrix<T>,
    phantom: PhantomData<CS>,
}


pub struct AffineMap<T, CS> {
    linear_coeffs: Matrix<T>,
    shift: Vector<T, CS>,
}


pub struct CoordinateMap<T, CSIn, CSOut> {
    linear_coeffs: Matrix<T>,
    origin_shift: Vector<T, CSOut>,
    phantom: PhantomData<CSIn>,
}


#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug, PartialEq)]
    struct World {}

    #[test]
    fn test_vector_basics() {

        let v: Vector<_, World> = Vector::new(&[1, 2, 3]);
        assert_eq!(v.dim(), 3);
        assert_eq!(v[0], 1);
        assert_eq!(v[2], 3);

        assert_eq!(v.iter().cloned().collect::<Vec<_>>(), vec![1, 2, 3]);
        assert_eq!(v.iter().map(|i| i * i).collect::<Vec<_>>(), vec![1, 4, 9]);
        assert_eq!(v.iter().cloned().collect::<Vector<_, World>>(), v);

        let mut v = v;
        v[1] = 4;
        assert_eq!(v, Vector::new(&[1, 4, 3]));

        assert_eq!(Vector::<i32, World>::zero(0), Vector::new(&[]));
        assert_eq!(Vector::<i32, World>::zero(5), Vector::new(&[0, 0, 0, 0, 0]));
        assert_eq!(Vector::<i32, World>::unit(3, 1), Vector::new(&[0, 1, 0]));
    }

    #[test]
    #[should_panic]
    fn test_bad_vector_unit() {
        Vector::<i32, World>::unit(3, 5);
    }

    #[test]
    fn test_vector_negation() {
        let v: Vector<_, World> = Vector::new(&[1, 2, 3]);

        assert_eq!(-&v, Vector::new(&[-1, -2, -3]));
        assert_eq!(-v, Vector::new(&[-1, -2, -3]));
    }

    #[test]
    fn test_vector_addition() {
        let mut v: Vector<_, World> = Vector::new(&[1, 2, 3]);

        v += Vector::new(&[-1, -1, -1]);
        assert_eq!(v, Vector::new(&[0, 1, 2]));

        v += &(v.clone());
        assert_eq!(v, Vector::new(&[0, 2, 4]));

        assert_eq!(v.clone() + v.clone(), Vector::new(&[0, 4, 8]));
        assert_eq!(&v.clone() + v.clone(), Vector::new(&[0, 4, 8]));
        assert_eq!(v.clone() + &v.clone(), Vector::new(&[0, 4, 8]));
        assert_eq!(&v.clone() + &v.clone(), Vector::new(&[0, 4, 8]));
    }

    #[test]
    #[should_panic]
    fn test_bad_vector_addition() {
        let mut v: Vector<_, World> = Vector::new(&[1, 2, 3]);
        v += Vector::new(&[-1, -1, -1, -1]);
    }

    #[test]
    fn test_vector_subtraction() {
        let mut v: Vector<_, World> = Vector::new(&[1, 2, 3]);

        v -= Vector::new(&[1, 1, 1]);
        assert_eq!(v, Vector::new(&[0, 1, 2]));

        v -= &(Vector::new(&[-1, -1, -2]));
        assert_eq!(v, Vector::new(&[1, 2, 4]));

        assert_eq!(v.clone() - Vector::unit(3, 1), Vector::new(&[1, 1, 4]));
        assert_eq!(&v - Vector::unit(3, 0), Vector::new(&[0, 2, 4]));
        assert_eq!(v.clone() - &Vector::unit(3, 2), Vector::new(&[1, 2, 3]));
        assert_eq!(&v - &Vector::unit(3, 0), Vector::new(&[0, 2, 4]));
    }

    #[test]
    #[should_panic]
    fn test_bad_vector_subtraction() {
        let mut v: Vector<_, World> = Vector::new(&[1, 2, 3]);
        v -= Vector::new(&[-1, -1, -1, -1]);
    }

    #[test]
    fn test_vector_multiplication() {
        let mut v: Vector<_, World> = Vector::new(&[1, 2, 3]);

        v *= 2;
        assert_eq!(v, Vector::new(&[2, 4, 6]));

        v *= &3;
        assert_eq!(v, Vector::new(&[6, 12, 18]));

        assert_eq!(v.clone() * 2, Vector::new(&[12, 24, 36]));
        assert_eq!(&v * 3, Vector::new(&[18, 36, 54]));
        assert_eq!(v.clone() * &4, Vector::new(&[24, 48, 72]));
        assert_eq!(&v * &5, Vector::new(&[30, 60, 90]));
    }
}
