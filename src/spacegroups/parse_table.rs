use std::collections::HashMap;
use std::io::{BufRead, Read, BufReader, Error};

use num_rational::Ratio;
use num_traits::{Zero, One};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::linear_algebra::Scalar;
use crate::arithmetic::matrices::Matrix;

use super::parse_operator::parse_operator;
use super::types::{CrystalSystem3d, Centering3d};


#[derive(Clone, Debug, PartialEq)]
struct Group {}


#[derive(Clone, Debug, PartialEq)]
struct Setting {}


type Fraction = Ratio<i32>;


struct Lookup {
    name: String,
    system: CrystalSystem3d,
    centering: Centering3d,
    from_std: CoordinateMap<Fraction, Group, Setting>,
}


impl CrystalSystem3d {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            "cubic" => Some(CrystalSystem3d::Cubic),
            "orthorhombic" => Some(CrystalSystem3d::Orthorhombic),
            "hexagonal" => Some(CrystalSystem3d::Hexagonal),
            "tetragonal" => Some(CrystalSystem3d::Tetragonal),
            "trigonal" => Some(CrystalSystem3d::Trigonal),
            "monoclinic" => Some(CrystalSystem3d::Monoclinic),
            "triclinic" => Some(CrystalSystem3d::Triclinic),
            _ => None,
        }
    }
}


impl Centering3d {
    pub fn from_string(s: &str) -> Option<Self> {
        match s {
            "P" => Some(Centering3d::Primitive),
            "F" => Some(Centering3d::FaceCentered),
            "I" => Some(Centering3d::BodyCentered),
            "R" => Some(Centering3d::Rhombohedral),
            "A" => Some(Centering3d::AFaceCentered),
            "B" => Some(Centering3d::BFaceCentered),
            "C" => Some(Centering3d::CFaceCentered),
            _ => None,
        }
    }
}


impl Scalar for Ratio<i32> {
    fn clear_column(col: usize, v: &mut Vec<Self>, b: &mut Vec<Self>) {
        let f = v[col] / b[col];
        v[col] = Self::zero();

        for k in (col + 1)..v.len() {
            v[k] -= &b[k] * &f;
        }
    }

    fn normalize_column(col: usize, v: &mut Vec<Self>) {
        let f = v[col];
        v[col] = Self::one();

        for k in (col + 1)..v.len() {
            v[k] /= &f;
        }
    }

    fn reduce_column(col: usize, v: &mut Vec<Self>, b: &Vec<Self>) {
        let f = v[col];
        v[col] = Self::zero();

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


impl<CS: Clone> AffineMap<Fraction, CS> {
    pub fn from_string(s: &str) -> Option<Self> {
        let (_, (linear, shift)) = parse_operator(s).ok()?;

        if linear.len() != shift.len() {
            None
        } else {
            let rows: Vec<_> = linear.iter().map(|r| Matrix::row(r)).collect();
            let a = Matrix::vstack(&rows);
            let t = Vector::new(&shift);

            Some(AffineMap::new(&a, &t))
        }
    }
}


pub fn parse_space_group_table<T: Read>(input: T) -> Result<(), Error> {
    let mut alias = HashMap::new();
    let mut lookup = vec![];

    for line in BufReader::new(input).lines() {
        let line = line?;
        let content = line.trim();

        if content.is_empty() || content.starts_with('#') {
            continue;
        } else if line.starts_with(' ') {
            let op = parse_operator(content);
            // add op to operator list in current setting
        } else {
            let fields: Vec<_> = content.split_whitespace().collect();

            if fields[0].to_lowercase() == "alias" {
                let key = fields[1].to_string();
                let val = fields[2].to_string();
                alias.insert(key, val);
            } else if fields[0].to_lowercase() == "lookup" {
                lookup.push(Lookup {
                    name: fields[1].to_string(),
                    system: CrystalSystem3d::from_string(fields[2]).unwrap(),
                    centering: Centering3d::from_string(fields[3]).unwrap(),
                    from_std: CoordinateMap::new(
                        &AffineMap::from_string(&fields[4..].join(" ")).unwrap()
                    ),
                });
            }
        }
    }

    Ok(())
}
