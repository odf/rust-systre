use std::ops::{SubAssign, AddAssign, MulAssign};

use num_traits::Zero;

use crate::arithmetic::geometry::{AffineMap, CoordinateMap};
use crate::arithmetic::lattices::reduced_lattice_basis;
use crate::arithmetic::matrices::Matrix;
use crate::arithmetic::linear_algebra::{LinearAlgebra, extend_basis};


#[derive(Debug, PartialEq)]
enum CrystalSystem2d {
    Oblique,
    Rectangular,
    Square,
    Hexagonal,
}

enum CrystalSystem3d {
    Cubic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}


enum Centering2d {
    Primitive,
    Centered,
}


enum Centering3d {
    Primitive,
    FaceCentered,
    BodyCentered,
    Rhombohedral,
    AFaceCentered,
    BFaceCentered,
    CFaceCentered,
}


pub trait Coord:
    crate::arithmetic::lattices::Coord
    + for <'a> AddAssign<&'a Self>
    + for <'a> SubAssign<&'a Self>
    + for <'a> MulAssign<&'a Self>
{
    fn epsilon() -> Self;
}


pub trait CoordPtr<T>:
    crate::arithmetic::lattices::CoordPtr<T>
{
}


struct OperatorDetails<T> {
    matrix: Matrix<T>,
    dimension: usize,
    axis: Option<Vec<T>>,
    order: usize,
    is_orientation_preserving: bool,
    is_clockwise: bool,
}

impl<T> OperatorDetails<T>
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    fn from(matrix: Matrix<T>) -> Self {
        let (nrows, ncols) = matrix.shape();
        assert_eq!(nrows, ncols);
        assert!(nrows >= 1 && nrows <= 3);

        let dimension = nrows;
        let is_orientation_preserving = !is_left_handed(&matrix);
        let order = matrix_order(&matrix, 6);
        let axis = operator_axis(&matrix);

        let is_clockwise = if dimension == 2 {
            if !is_orientation_preserving {
                false
            } else if order == 0 || order > 2 {
                !(matrix[(0, 1)]).is_negative()
            } else {
                true
            }
        } else if (order == 0 || order > 2) && axis.is_some() {
            let axis = axis.clone().unwrap();
            let u0 = [T::one(), T::zero(), T::zero()];
            let u1 = [T::zero(), T::one(), T::zero()];
            let v = if are_collinear(&axis, &u0) { u1 } else { u0 };
            let v = Matrix::col(&v);
            let w = &matrix * &v;
            let a = Matrix::hstack(&[Matrix::col(&axis), v, w]);

            !is_left_handed(&a)
        } else {
            true
        };

        OperatorDetails {
            matrix,
            dimension,
            axis,
            order,
            is_orientation_preserving,
            is_clockwise
        }
    }
}


fn positive_direction<T>(v: &[T]) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    if let Some(x) = v.iter().find(|x| !x.is_zero()) {
        if x.is_negative() {
            return v.iter().map(|x| -x).collect();
        }
    }
    v.to_vec()
}


fn dot<T>(v: &[T], w: &[T]) -> T
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    assert_eq!(v.len(), w.len());

    let mut s = T::zero();
    for i in 0..v.len() {
        s = s + &v[i] * &w[i];
    }

    s
}


fn are_orthogonal<T>(v: &[T], w: &[T]) -> bool
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    dot(v, w) == T::zero()
}


fn are_collinear<T>(v: &[T], w: &[T]) -> bool
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    dot(v, w) * dot(v, w) == dot(v, v) * dot(w, w)
}


fn matrix_order<T>(matrix: &Matrix<T>, max: usize) -> usize
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    let (nrows, ncols) = matrix.shape();
    assert_eq!(nrows, ncols);

    let identity = Matrix::identity(nrows);
    let mut a = identity.clone();

    for i in 1..=max {
        a = &a * matrix;
        if a == identity {
            return i;
        }
    }
    
    0
}


fn is_left_handed<T>(matrix: &Matrix<T>) -> bool
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    matrix.determinant().is_negative()
}


fn operator_axis<T>(matrix: &Matrix<T>) -> Option<Vec<T>>
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    let (dim, ncols) = matrix.shape();
    assert_eq!(dim, ncols);

    let r = if dim % 2 != 0 && is_left_handed(matrix) {
        matrix + Matrix::identity(dim)
    } else {
        matrix - Matrix::identity(dim)
    };

    if let Some(z) = r.null_space() {
        if z.ncols == 1 {
            return Some(positive_direction(&z.get_col(0)));
        }
    }

    None
}


pub struct SpaceGroup2d {
    dimension: usize,
    crystal_system: CrystalSystem2d,
    centering: Centering2d,
    full_name: String,
    group_name: String,
    extension: String,
}


pub fn identitify_spacegroup_2d<T, CSIn, CSOut>(ops: &[AffineMap<T, CSIn>])
    -> (SpaceGroup2d, CoordinateMap<T, CSIn, CSOut>)
    where
        CSIn: Clone + PartialEq,
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    let lin_ops: Vec<_> = ops.iter().map(|op| op.linear_matrix()).collect();
    let (crystal_system, basis) = crystal_system_and_basis_2d(&lin_ops);
    let to_basis = Matrix::hstack(&basis).inverse().unwrap();
    let pcell: Vec<_> = primitive_cell(ops).iter()
        .map(|b| &to_basis * b).collect();
    let (normalized, centering) = normalized_basis2d(crystal_system, &pcell);

    todo!()
}


fn crystal_system_and_basis_2d<T>(ops: &[Matrix<T>])
    -> (CrystalSystem2d, Vec<Matrix<T>>)
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    let ops_with_details: Vec<_> = ops.iter()
        .map(|op| OperatorDetails::from(op.clone()))
        .collect();
    let mirrors: Vec<_> = ops_with_details.iter()
        .filter(|opd| !opd.is_orientation_preserving)
        .collect();
    let n = ops_with_details.iter().map(|opd| opd.order).max().unwrap();
    let m = if n == 6 { 3 } else { n };

    let crystal_system = match n {
        3 | 6 => CrystalSystem2d::Hexagonal,
        4 => CrystalSystem2d::Square,
        _ => if mirrors.len() > 0 {
            CrystalSystem2d::Rectangular
        } else {
            CrystalSystem2d::Oblique
        }
    };

    let x = if mirrors.len() > 0 {
        Matrix::col(&mirrors[0].axis.clone().unwrap())
    } else {
        Matrix::col(&[T::one(), T::zero()])
    };

    let y = if m >= 3 {
        let s = ops_with_details.iter().find(|s|
            s.is_orientation_preserving && s.is_clockwise && s.order == m
        ).unwrap();
        &s.matrix * &x
    } else if mirrors.len() > 1 {
        Matrix::col(&mirrors[1].axis.clone().unwrap())
    } else {
        let t = if x[(0, 0)] == T::zero() {
            Matrix::col(&[T::one(), T::zero()])
        } else {
            Matrix::col(&[T::zero(), T::one()])
        };
        if mirrors.len() > 0 {
            &t - &mirrors[0].matrix * &t
        } else {
            t
        }
    };

    let b = Matrix::hstack(&[x.clone(), y.clone()]);
    let basis = if is_left_handed(&b) { vec![x, -y] } else { vec![x, y] };

    (crystal_system, basis)
}


fn primitive_cell<T, CS>(ops: &[AffineMap<T, CS>]) -> Vec<Matrix<T>>
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        CS: Clone + PartialEq
{
    let dim = ops[0].dim();
    let identity: Matrix<T> = Matrix::identity(dim);

    let mut cell = identity.get_rows();

    for op in ops {
        let has_shift = op.shift().iter().any(|x| !x.is_zero());
        if has_shift && op.linear_matrix() == identity {
            let s: Vec<_> = op.shift().into_iter().collect();
            extend_basis(&s, &mut cell);
        }
    }

    as_columns(cell)
}


fn normalized_basis2d<T>(crys: CrystalSystem2d, basis_in: &Vec<Matrix<T>>)
    -> (Vec<Matrix<T>>, Centering2d)
    where
        T: Coord, for <'a> &'a T: CoordPtr<T>,
        Matrix<T>: LinearAlgebra<T>
{
    let vs: Vec<Vec<T>> = from_columns(basis_in.clone());
    let b = reduced_lattice_basis(&vs, dot, &T::epsilon()).unwrap();

    match crys {
        CrystalSystem2d::Oblique => (
            as_columns(b),
            Centering2d::Primitive
        ),
        CrystalSystem2d::Rectangular => {
            if !b[0][0].is_zero() && !b[0][1].is_zero() {
                (
                    as_columns(vec![
                        vec![T::zero(), &b[0][1] + &b[0][1]],
                        vec![&b[0][0] + &b[0][0], T::zero()]
                    ]),
                    Centering2d::Centered
                )
            } else if !b[1][0].is_zero() && !b[1][1].is_zero() {
                (
                    as_columns(vec![
                        vec![T::zero(), &b[1][1] + &b[1][1]],
                        vec![&b[1][0] + &b[1][0], T::zero()]
                    ]),
                    Centering2d::Centered
                )
            } else if b[0][1].is_zero() {
                (
                    as_columns(vec![
                        b[1].clone(),
                        vec![-&b[0][0], -&b[0][1]]
                    ]),
                    Centering2d::Primitive
                )
            } else {
                (as_columns(b), Centering2d::Primitive)
            }
        },
        CrystalSystem2d::Square => (
            as_columns(vec![
                b[0].clone(),
                vec![-&b[0][1], b[0][0].clone()]
            ]),
            Centering2d::Primitive
        ),
        CrystalSystem2d::Hexagonal => (
            as_columns(vec![
                b[0].clone(),
                vec![-&b[0][1], &b[0][0] - &b[0][1]]
            ]),
            Centering2d::Primitive
        ),
    }
}


fn from_columns<T>(cols: Vec<Matrix<T>>) -> Vec<Vec<T>>
    where T: Clone
{
    cols.iter().map(|v| v.iter().cloned().collect()).collect()
}


fn as_columns<T>(vs: Vec<Vec<T>>) -> Vec<Matrix<T>>
    where T: Clone
{
    vs.iter().map(|v| Matrix::col(v)).collect()
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::arithmetic::matrices::Matrix;

    use super::*;

    impl Coord for BigRational {
        fn epsilon() -> Self {
            BigRational::zero()
        }
    }

    impl CoordPtr<BigRational> for &BigRational {
    }

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    #[test]
    fn test_spacegroups_operator_details_1d() {
        let m = Matrix::new(1, &[r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 1);
        assert_eq!(&opd.axis, &Some(vec![r(1)]));
        assert_eq!(opd.order, 1);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(1, &[r(-1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 1);
        assert_eq!(&opd.axis, &Some(vec![r(1)]));
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, true);
    }

    #[test]
    fn test_spacegroups_operator_details_2d() {
        let m = Matrix::new(2, &[r(1), r(0), r(0), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 1);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(2, &[r(1), r(0), r(0), r(-1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(0)]));
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);

        let m = Matrix::new(2, &[r(0), r(1), r(-1), r(0)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 4);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(2, &[r(0), r(-1), r(1), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 6);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, false);

        let m = Matrix::new(2, &[r(-1), r(0), r(1), r(1)]);
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 2);
        assert_eq!(&opd.axis, &Some(vec![r(0), r(1)]));
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);
    }

    #[test]
    fn test_spacegroups_operator_details_3d() {
        let m = Matrix::new(
            3, &[r(1), r(0), r(0), r(0), r(1), r(0), r(0), r(0), r(1)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 1);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(-1), r(0), r(0), r(0), r(-1), r(0), r(0), r(0), r(-1)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(opd.axis, None);
        assert_eq!(opd.order, 2);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(0), r(0), r(1), r(1), r(0), r(0), r(0), r(1), r(0)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(1), r(1)]));
        assert_eq!(opd.order, 3);
        assert_eq!(opd.is_orientation_preserving, true);
        assert_eq!(opd.is_clockwise, true);

        let m = Matrix::new(
            3, &[r(0), r(0), r(-1), r(-1), r(0), r(0), r(0), r(-1), r(0)]
        );
        let opd = OperatorDetails::from(m.clone());
        assert_eq!(&opd.matrix, &m);
        assert_eq!(opd.dimension, 3);
        assert_eq!(&opd.axis, &Some(vec![r(1), r(1), r(1)]));
        assert_eq!(opd.order, 6);
        assert_eq!(opd.is_orientation_preserving, false);
        assert_eq!(opd.is_clockwise, false);
    }

    #[test]
    fn test_spacegroups_crystal_system_2d() {
        let ops = vec![
            Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
            Matrix::new(2, &[r(-1), r(0), r(0), r(-1)])
        ];
        let (cs, b) = crystal_system_and_basis_2d(&ops);
        assert_eq!(CrystalSystem2d::Oblique, cs);
        assert_eq!(
            &b,
            &vec![Matrix::col(&[r(1), r(0)]),Matrix::col(&[r(0), r(1)])]
        );

        let ops = vec![
            Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
            Matrix::new(2, &[r(1), r(0), r(0), r(-1)])
        ];
        let (cs, _) = crystal_system_and_basis_2d(&ops);
        assert_eq!(CrystalSystem2d::Rectangular, cs);

        let ops = vec![
            Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
            Matrix::new(2, &[r(0), r(1), r(-1), r(-1)]).transpose(),
            Matrix::new(2, &[r(-1), r(-1), r(1), r(0)]).transpose(),
        ];
        let (cs, _) = crystal_system_and_basis_2d(&ops);
        assert_eq!(CrystalSystem2d::Hexagonal, cs);

        let ops = vec![
            Matrix::new(2, &[r(1), r(0), r(0), r(1)]),
            Matrix::new(2, &[r(0), r(1), r(-1), r(0)]).transpose(),
            Matrix::new(2, &[r(-1), r(0), r(0), r(-1)]).transpose(),
            Matrix::new(2, &[r(0), r(-1), r(1), r(0)]).transpose(),
        ];
        let (cs, _) = crystal_system_and_basis_2d(&ops);
        assert_eq!(CrystalSystem2d::Square, cs);
    }
}
