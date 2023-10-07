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
pub enum CrystalSystem {
    Oblique2d,
    Rectangular2d,
    Square2d,
    Hexagonal2d,
    Cubic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}


#[derive(Debug, PartialEq)]
pub enum Centering {
    Primitive2d,
    Centered2d,
    Primitive,
    FaceCentered,
    BodyCentered,
    Rhombohedral,
    AFaceCentered,
    BFaceCentered,
    CFaceCentered,
}


pub struct SpaceGroup {
    dimension: usize,
    crystal_system: CrystalSystem,
    centering: Centering,
    full_name: String,
    group_name: String,
    extension: String,
}


impl Display for CrystalSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CrystalSystem::Oblique2d => write!(f, "oblique")?,
            CrystalSystem::Rectangular2d => write!(f, "rectangular")?,
            CrystalSystem::Square2d => write!(f, "square")?,
            CrystalSystem::Hexagonal2d => write!(f, "hexagonal")?,
            CrystalSystem::Cubic => write!(f, "cubic")?,
            CrystalSystem::Orthorhombic => write!(f, "orthorhombic")?,
            CrystalSystem::Hexagonal => write!(f, "hexagonal")?,
            CrystalSystem::Tetragonal => write!(f, "tetragonal")?,
            CrystalSystem::Trigonal => write!(f, "trigonal")?,
            CrystalSystem::Monoclinic => write!(f, "monoclinic")?,
            CrystalSystem::Triclinic => write!(f, "triclinic")?,
        }
        Ok(())
    }
}


impl Display for Centering {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Centering::Primitive2d => write!(f, "p")?,
            Centering::Centered2d => write!(f, "c")?,
            Centering::Primitive => write!(f, "P")?,
            Centering::FaceCentered => write!(f, "F")?,
            Centering::BodyCentered => write!(f, "I")?,
            Centering::Rhombohedral => write!(f, "R")?,
            Centering::AFaceCentered => write!(f, "A")?,
            Centering::BFaceCentered => write!(f, "B")?,
            Centering::CFaceCentered => write!(f, "C")?,
        }
        Ok(())
    }
}
