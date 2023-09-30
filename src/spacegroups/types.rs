use std::{ops::{Add, Sub, Mul, Neg, SubAssign, AddAssign, MulAssign}, fmt::Display};

use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use num_rational::{Ratio, BigRational};

use crate::arithmetic::linear_algebra::Scalar;


pub trait Coord:
    Clone + PartialOrd + Signed + Scalar
    + for <'a> Add<&'a Self, Output=Self>
    + for <'a> AddAssign<&'a Self>
    + for <'a> SubAssign<&'a Self>
    + for <'a> MulAssign<&'a Self>
{
    fn round(&self) -> Self;
    fn div_rounded(&self, other: &Self) -> Self;
    fn epsilon() -> Self;
    fn promote(q: Ratio<i32>) -> Self;
}


pub trait CoordPtr<T>:
    Sized
    + Add<Output=T> + Sub<Output=T> + Sub<T, Output=T> + Neg<Output=T>
    + Mul<Output=T>
{
}


impl Coord for BigRational {
    fn round(&self) -> Self {
        BigRational::round(self)
    }

    fn div_rounded(&self, other: &Self) -> Self {
        BigRational::round(&(self / other))
    }

    fn epsilon() -> Self {
        BigRational::zero()
    }

    fn promote(q: Ratio<i32>) -> Self {
        Ratio::new(BigInt::from(*q.numer()), BigInt::from(*q.denom()))
    }
}

impl CoordPtr<BigRational> for &BigRational {
}


#[derive(Debug, PartialEq)]
pub enum CrystalSystem2d {
    Oblique,
    Rectangular,
    Square,
    Hexagonal,
}

pub enum CrystalSystem3d {
    Cubic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}

pub enum Centering2d {
    Primitive,
    Centered,
}

pub enum Centering3d {
    Primitive,
    FaceCentered,
    BodyCentered,
    Rhombohedral,
    AFaceCentered,
    BFaceCentered,
    CFaceCentered,
}

pub struct SpaceGroup2d {
    dimension: usize,
    crystal_system: CrystalSystem2d,
    centering: Centering2d,
    full_name: String,
    group_name: String,
    extension: String,
}


impl Display for CrystalSystem2d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CrystalSystem2d::Oblique => write!(f, "oblique")?,
            CrystalSystem2d::Rectangular => write!(f, "rectangular")?,
            CrystalSystem2d::Square => write!(f, "square")?,
            CrystalSystem2d::Hexagonal => write!(f, "hexagonal")?,
        }
        Ok(())
    }
}


impl Display for CrystalSystem3d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CrystalSystem3d::Cubic => write!(f, "cubic")?,
            CrystalSystem3d::Orthorhombic => write!(f, "orthorhombic")?,
            CrystalSystem3d::Hexagonal => write!(f, "hexagonal")?,
            CrystalSystem3d::Tetragonal => write!(f, "tetragonal")?,
            CrystalSystem3d::Trigonal => write!(f, "trigonal")?,
            CrystalSystem3d::Monoclinic => write!(f, "monoclinic")?,
            CrystalSystem3d::Triclinic => write!(f, "triclinic")?,
        }
        Ok(())
    }
}


impl Display for Centering2d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Centering2d::Primitive => write!(f, "p")?,
            Centering2d::Centered => write!(f, "c")?,
        }
        Ok(())
    }
}


impl Display for Centering3d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Centering3d::Primitive => write!(f, "P")?,
            Centering3d::FaceCentered => write!(f, "F")?,
            Centering3d::BodyCentered => write!(f, "I")?,
            Centering3d::Rhombohedral => write!(f, "R")?,
            Centering3d::AFaceCentered => write!(f, "A")?,
            Centering3d::BFaceCentered => write!(f, "B")?,
            Centering3d::CFaceCentered => write!(f, "C")?,
        }
        Ok(())
    }
}
