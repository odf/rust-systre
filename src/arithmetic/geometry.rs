use std::marker::PhantomData;
use std::ops::{Neg, Add, AddAssign, Sub, SubAssign};
use std::ops::{Mul, MulAssign, Div, DivAssign};
use std::ops::{Index, IndexMut};
use std::slice::Iter;

use num_traits::{Zero, One};

use super::matrices::Matrix;
use super::linear_algebra::Field;


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

impl<S, T, CS> Add<S> for &Vector<T, CS>
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

impl<S, T, CS> Sub<S> for &Vector<T, CS>
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

impl<S, T, CS> Mul<S> for &Vector<T, CS>
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

impl<T, CS> DivAssign<&T> for Vector<T, CS>
    where T: for <'a> DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: &T) {
        self.coords /= rhs;
    }
}

impl<T, CS> DivAssign<T> for Vector<T, CS>
    where T: for <'a> DivAssign<&'a T>
{
    fn div_assign(&mut self, rhs: T) {
        self.coords /= rhs;
    }
}

impl<S, T, CS> Div<S> for &Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: DivAssign<S>
{
    type Output = Vector<T, CS>;

    fn div(self, rhs: S) -> Vector<T, CS> {
        let mut tmp = (*self).clone();
        tmp /= rhs;
        tmp
    }
}

impl<S, T, CS> Div<S> for Vector<T, CS>
    where T: Clone, CS: Clone, Vector<T, CS>: DivAssign<S>
{
    type Output = Vector<T, CS>;

    fn div(self, rhs: S) -> Vector<T, CS> {
        &self / rhs
    }
}


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Point<T, CS> {
    coords: Matrix<T>,
    phantom: PhantomData<CS>,
}

impl<T, CS> Point<T, CS> {
    pub fn dim(&self) -> usize {
        self.coords.shape().0
    }

    pub fn iter(&self) -> Iter<T> {
        self.coords.iter()
    }
}

impl<T, CS> Index<usize> for Point<T, CS> {
    type Output = T;

    fn index(&self, i: usize) -> &T {
        &self.coords[(i, 0)]
    }
}

impl<T, CS> IndexMut<usize> for Point<T, CS> {
    fn index_mut(&mut self, i: usize) -> &mut T {
        &mut self.coords[(i, 0)]
    }
}

impl<T: Clone, CS> FromIterator<T> for Point<T, CS> {
    fn from_iter<I: IntoIterator<Item=T>>(iter: I) -> Self {
        Point::new(&Vec::from_iter(iter))
    }
}

impl<T: Clone, CS> From<Matrix<T>> for Point<T, CS> {
    fn from(coords: Matrix<T>) -> Self {
        Point { coords, phantom: PhantomData::default() }
    }
}

impl<T: Clone, CS> Point<T, CS> {
    pub fn new(coords: &[T]) -> Point<T, CS> {
        Point {
            coords: Matrix::col(coords),
            phantom: PhantomData::default()
        }
    }
}

impl<T: Clone + Zero, CS> Point<T, CS> {
    pub fn origin(dim: usize) -> Point<T, CS> {
        Point::new(&vec![T::zero(); dim])
    }
}

impl<T,  CS> AddAssign<Vector<T,  CS>> for Point<T,  CS>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords += &rhs.coords
    }
}

impl<T,  CS> AddAssign<&Vector<T,  CS>> for Point<T,  CS>
    where T: for <'a> AddAssign<&'a T>
{
    fn add_assign(&mut self, rhs: &Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords += &rhs.coords
    }
}

impl<S, T, CS> Add<S> for &Point<T, CS>
    where T: Clone, CS: Clone, Point<T, CS>: AddAssign<S>
{
    type Output = Point<T, CS>;

    fn add(self, rhs: S) -> Point<T, CS> {
        let mut tmp = (*self).clone();
        tmp += rhs;
        tmp
    }
}

impl<S, T, CS> Add<S> for Point<T, CS>
    where T: Clone, CS: Clone, Point<T, CS>: AddAssign<S>
{
    type Output = Point<T, CS>;

    fn add(self, rhs: S) -> Point<T, CS> {
        &self + rhs
    }
}

impl<T,  CS> SubAssign<Vector<T,  CS>> for Point<T,  CS>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords -= &rhs.coords
    }
}

impl<T,  CS> SubAssign<&Vector<T,  CS>> for Point<T,  CS>
    where T: for <'a> SubAssign<&'a T>
{
    fn sub_assign(&mut self, rhs: &Vector<T,  CS>) {
        assert_eq!(self.dim(), rhs.dim());
        self.coords -= &rhs.coords
    }
}

impl<'a, S, T: 'a, CS> Sub<S> for &'a Point<T, CS>
    where T: Clone, CS: Clone, Point<T, CS>: SubAssign<S>
{
    type Output = Point<T, CS>;

    fn sub(self, rhs: S) -> Point<T, CS> {
        let mut tmp = (*self).clone();
        tmp -= rhs;
        tmp
    }
}

impl<S, T, CS> Sub<S> for Point<T, CS>
    where T: Clone, CS: Clone, Point<T, CS>: SubAssign<S>
{
    type Output = Point<T, CS>;

    fn sub(self, rhs: S) -> Point<T, CS> {
        &self - rhs
    }
}

impl <T, CS> Sub<Point<T, CS>> for &Point<T, CS>
    where T: Clone, Matrix<T>: SubAssign<Matrix<T>>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: Point<T, CS>) -> Self::Output {
        let mut coords = self.coords.clone();
        coords -= rhs.coords;
        Vector::new(&coords.data)
    }
}

impl <T, CS> Sub<&Point<T, CS>> for &Point<T, CS>
    where T: Clone, for <'a> Matrix<T>: SubAssign<&'a Matrix<T>>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: &Point<T, CS>) -> Self::Output {
        let mut coords = self.coords.clone();
        coords -= &rhs.coords;
        Vector::new(&coords.data)
    }
}

impl <T, CS> Sub<Point<T, CS>> for Point<T, CS>
    where T: Clone, Matrix<T>: SubAssign<Matrix<T>>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: Point<T, CS>) -> Self::Output {
        &self - rhs
    }
}

impl <T, CS> Sub<&Point<T, CS>> for Point<T, CS>
    where T: Clone, for <'a> Matrix<T>: SubAssign<&'a Matrix<T>>
{
    type Output = Vector<T, CS>;

    fn sub(self, rhs: &Point<T, CS>) -> Self::Output {
        &self - rhs
    }
}


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ScalarProduct<T, CS> {
    coeffs: Matrix<T>,
    phantom: PhantomData<CS>,
}

impl<T, CS> ScalarProduct<T, CS> {
    pub fn dim(&self) -> usize {
        self.coeffs.shape().0
    }
}

impl<T, CS> Index<(usize, usize)> for ScalarProduct<T, CS> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &T {
        &self.coeffs[(i, j)]
    }
}

impl<T, CS> From<Matrix<T>> for ScalarProduct<T, CS>
    where
        T: std::fmt::Debug + PartialEq + PartialOrd,
        T: Field + Clone + Sub<Output=T> + SubAssign,
        for <'a> T: Div<&'a T, Output=T> + Mul<&'a T, Output=T>,
        for <'a> T: AddAssign<&'a T> + DivAssign<&'a T> + MulAssign<&'a T>,
        for <'a> &'a T: Neg<Output=T> + Mul<&'a T, Output=T>,
{
    fn from(coeffs: Matrix<T>) -> Self {
        assert_eq!(&coeffs, &coeffs.transpose());
        for i in 1..=coeffs.shape().0 {
            assert!(coeffs.submatrix(0..i, 0..i).determinant() > T::zero());
        }

        ScalarProduct { coeffs, phantom: PhantomData::default() }
    }
}

impl<T, CS> ScalarProduct<T, CS> where T: Clone + Zero + One {
    pub fn default(dim: usize) -> ScalarProduct<T, CS> {
        ScalarProduct {
            coeffs: Matrix::identity(dim),
            phantom: PhantomData::default()
        }
    }
}

impl<T, CS> ScalarProduct<T, CS>
    where
        T: Clone + Zero + One + AddAssign<T>,
        for <'a> T: Mul<&'a T, Output=T>,
        for <'a> &'a T: Mul<&'a T, Output=T>
{
    pub fn apply(&self, u: &Vector<T, CS>, v: &Vector<T, CS>) -> T {
        assert_eq!(u.dim(), self.dim());
        assert_eq!(v.dim(), self.dim());

        let mut t = T::zero();
        for i in 0..u.dim() {
            for j in 0..u.dim() {
                t += &u[i] * &self[(i, j)] * &v[j];
            }
        }
        t
    }
}


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AffineMap<T, CS> {
    linear_coeffs: Matrix<T>,
    shift: Vector<T, CS>,
}

impl<T, CS> AffineMap<T, CS> {
    pub fn dim(&self) -> usize {
        self.shift.dim()
    }
}

impl<T: Clone, CS: Clone> AffineMap<T, CS> {
    pub fn new(linear_coeffs: &Matrix<T>, shift: &Vector<T, CS>) -> Self
    {
        let n = shift.dim();
        assert_eq!(linear_coeffs.shape().0, n);
        assert_eq!(linear_coeffs.shape().1, n);

        AffineMap {
            linear_coeffs: linear_coeffs.clone(),
            shift: shift.clone()
        }
    }

    pub fn identity(dim: usize) -> Self
        where T: Zero + One
    {
        AffineMap::new(&Matrix::identity(dim), &Vector::zero(dim))
    }
}

impl<T: Clone, CS: Clone> Mul<&AffineMap<T, CS>> for &AffineMap<T, CS>
    where
        for <'a> &'a Matrix<T>: Mul<&'a Matrix<T>, Output=Matrix<T>>,
        for <'a> Matrix<T>: Add<&'a Matrix<T>, Output=Matrix<T>>
{
    type Output = AffineMap<T, CS>;

    fn mul(self, rhs: &AffineMap<T, CS>) -> Self::Output {
        assert_eq!(self.dim(), rhs.dim());
        let linear_coeffs = &self.linear_coeffs * &rhs.linear_coeffs;
        let shift = Vector::from(
            &self.linear_coeffs * &rhs.shift.coords + &self.shift.coords
        );
        AffineMap::new(&linear_coeffs, &shift)
    }
}

impl<T: Clone, CS> Mul<&Point<T, CS>> for &AffineMap<T, CS>
    where
        for <'a> &'a Matrix<T>: Mul<&'a Matrix<T>, Output=Matrix<T>>,
        for <'a> Matrix<T>: Add<&'a Matrix<T>, Output=Matrix<T>>
{
    type Output = Point<T, CS>;

    fn mul(self, rhs: &Point<T, CS>) -> Self::Output {
        Point::from(&self.linear_coeffs * &rhs.coords + &self.shift.coords)
    }
}

impl<T: Clone, CS> Mul<&Vector<T, CS>> for &AffineMap<T, CS>
    where for <'a> &'a Matrix<T>: Mul<&'a Matrix<T>, Output=Matrix<T>>
{
    type Output = Vector<T, CS>;

    fn mul(self, rhs: &Vector<T, CS>) -> Self::Output {
        Vector::from(&self.linear_coeffs * &rhs.coords)
    }
}

impl<S, T, CS> Mul<S> for &AffineMap<T, CS>
    where for <'a> &'a AffineMap<T, CS>: Mul<&'a S, Output=S>
{
    type Output = S;

    fn mul(self, rhs: S) -> Self::Output {
        self * &rhs
    }
}

impl<S, T, U, CS> Mul<S> for AffineMap<T, CS>
    where for <'a> &'a AffineMap<T, CS>: Mul<S, Output=U>
{
    type Output = U;

    fn mul(self, rhs: S) -> Self::Output {
        &self * rhs
    }
}


pub struct CoordinateMap<T, CSIn, CSOut> {
    linear_coeffs: Matrix<T>,
    origin_shift: Vector<T, CSOut>,
    phantom: PhantomData<CSIn>,
}


/// # Examples
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
/// ```
#[cfg(doctest)]
struct CannotAddVectorsWithDifferentCoordinateSystems();

/// # Examples
/// 
/// ```compile_fail
/// use rust_systre::arithmetic::geometry::*;
/// 
/// #[derive(Debug, PartialEq)]
/// struct World {}
/// 
/// let mut u = Point::<_, World>::new(&[1, 2, 3]);
/// u += Point::<_, World>::new(&[3, 2, 1]);
/// ```
#[cfg(doctest)]
struct CannotAddTwoPoints();


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
    fn test_point_basics() {

        let p: Point<_, World> = Point::new(&[1, 2, 3]);
        assert_eq!(p.dim(), 3);
        assert_eq!(p[0], 1);
        assert_eq!(p[2], 3);

        assert_eq!(p.iter().cloned().collect::<Vec<_>>(), vec![1, 2, 3]);
        assert_eq!(p.iter().map(|i| i * i).collect::<Vec<_>>(), vec![1, 4, 9]);
        assert_eq!(p.iter().cloned().collect::<Point<_, World>>(), p);

        let mut p = p;
        p[1] = 4;
        assert_eq!(p, Point::new(&[1, 4, 3]));

        assert_eq!(Point::<i32, World>::origin(0), Point::new(&[]));
        assert_eq!(Point::<i32, World>::origin(4), Point::new(&[0, 0, 0, 0]));
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
    fn test_point_vector_addition() {
        let mut p: Point<_, World> = Point::new(&[1, 2, 3]);
        let v: Vector<_, World> = Vector::new(&[1, 0, 1]);

        p += Vector::unit(3, 1);
        assert_eq!(p, Point::new(&[1, 3, 3]));

        p += &v;
        assert_eq!(p, Point::new(&[2, 3, 4]));

        assert_eq!(p.clone() + v.clone(), Point::new(&[3, 3, 5]));
        assert_eq!(&p.clone() + v.clone(), Point::new(&[3, 3, 5]));
        assert_eq!(p.clone() + &v.clone(), Point::new(&[3, 3, 5]));
        assert_eq!(&p.clone() + &v.clone(), Point::new(&[3, 3, 5]));
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
    fn test_point_vector_subtraction() {
        let mut p: Point<_, World> = Point::new(&[1, 2, 3]);

        p -= Vector::new(&[1, 1, 1]);
        assert_eq!(p, Point::new(&[0, 1, 2]));

        p -= &(Vector::new(&[-1, -1, -2]));
        assert_eq!(p, Point::new(&[1, 2, 4]));

        assert_eq!(p.clone() - Vector::unit(3, 1), Point::new(&[1, 1, 4]));
        assert_eq!(&p - Vector::unit(3, 0), Point::new(&[0, 2, 4]));
        assert_eq!(p.clone() - &Vector::unit(3, 2), Point::new(&[1, 2, 3]));
        assert_eq!(&p - &Vector::unit(3, 0), Point::new(&[0, 2, 4]));
    }

    #[test]
    fn test_point_point_subtraction() {
        let p: Point<_, World> = Point::new(&[1, 2, 3]);
        let q: Point<_, World> = Point::new(&[1, 0, 1]);

        assert_eq!(p.clone() - q.clone(), Vector::new(&[0, 2, 2]));
        assert_eq!(p.clone() - &q.clone(), Vector::new(&[0, 2, 2]));
        assert_eq!(&p.clone() - q.clone(), Vector::new(&[0, 2, 2]));
        assert_eq!(&p.clone() - &q.clone(), Vector::new(&[0, 2, 2]));
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

    #[test]
    fn test_vector_division() {
        let mut v: Vector<_, World> = Vector::new(&[73, 79, 83]);

        v /= 2;
        assert_eq!(v, Vector::new(&[36, 39, 41]));

        v /= &3;
        assert_eq!(v, Vector::new(&[12, 13, 13]));

        assert_eq!(v.clone() / 2, Vector::new(&[6, 6, 6]));
        assert_eq!(&v / 3, Vector::new(&[4, 4, 4]));
        assert_eq!(v.clone() / &4, Vector::new(&[3, 3, 3]));
        assert_eq!(&v / &5, Vector::new(&[2, 2, 2]));
    }

    #[test]
    fn test_scalar_product() {
        let dot: ScalarProduct<f64, World> = ScalarProduct::default(3);
        let v = Vector::new(&[1.0, 2.0, 3.0]);
        let w = Vector::new(&[-1.0, 2.0, -1.0]);

        assert_eq!(dot.apply(&v, &v), 14.0);
        assert_eq!(dot.apply(&v, &w), 0.0);

        let dot: ScalarProduct<_, World> = ScalarProduct::from(
            Matrix::new(2, &[1.0, -0.5, -0.5, 1.0])
        );
        let v = Vector::new(&[1.0, 0.0]);
        let w = Vector::new(&[1.0, 1.0]);
        assert_eq!(dot.apply(&v, &w), 0.5);
    }

    #[test]
    #[should_panic]
    fn test_bad_scalar_product_not_symmetric() {
        let _dot: ScalarProduct<_, World> = ScalarProduct::from(
            Matrix::new(2, &[1.0, -0.5, -0.51, 1.0])
        );
    }

    #[test]
    #[should_panic]
    fn test_bad_scalar_product_not_positive_definite() {
        let _dot: ScalarProduct<_, World> = ScalarProduct::from(
            Matrix::new(2, &[1.0, -0.5, -0.5, 0.2])
        );
    }

    #[test]
    fn test_affine_map() {
        let a: AffineMap<_, World> = AffineMap::new(
            &Matrix::new(2, &[0, -1, 1, 0]), &Vector::new(&[1, 0])
        );
        let b: AffineMap<_, World> = AffineMap::new(
            &Matrix::new(2, &[0, 1, -1, 0]), &Vector::new(&[0, 1])
        );
        assert_eq!(&a * &b, AffineMap::identity(2));
        assert_eq!(&a * b.clone(), AffineMap::identity(2));
        assert_eq!(a.clone() * &b, AffineMap::identity(2));
        assert_eq!(a.clone() * b.clone(), AffineMap::identity(2));
    }
}
