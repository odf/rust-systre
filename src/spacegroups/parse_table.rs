use std::collections::HashMap;
use std::io::{BufRead, Read, BufReader, Error};

use num_rational::{Ratio, BigRational};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::matrices::Matrix;

use super::parse_operator::parse_operator;
use super::types::{
    Coord, CrystalSystem2d, Centering2d, CrystalSystem3d, Centering3d
};


#[derive(Clone, Debug, PartialEq)]
struct Group {}


#[derive(Clone, Debug, PartialEq)]
struct Setting {}


struct TableEntry {
    name: String,
    canonical_name: String,
    transform: CoordinateMap<BigRational, Group, Setting>,
    operators: Vec<AffineMap<BigRational, Setting>>,
}


enum Lookup {
    Entry2d {
        name: String,
        system: CrystalSystem2d,
        centering: Centering2d,
        from_std: CoordinateMap<BigRational, Group, Setting>,
    },
    Entry3d {
        name: String,
        system: CrystalSystem3d,
        centering: Centering3d,
        from_std: CoordinateMap<BigRational, Group, Setting>,
    },
}


impl CrystalSystem2d {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            _ => None,
        }
    }
}


impl Centering2d {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            _ => None,
        }
    }
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


impl<CS: Clone> AffineMap<BigRational, CS> {
    pub fn from_string(s: &str) -> Option<Self> {
        let (_, (linear, shift)) = parse_operator(s).ok()?;

        if linear.len() != shift.len() {
            None
        } else {
            let linear: Vec<_> = linear.iter().map(promote_vec).collect();
            let rows: Vec<_> = linear.iter().map(|r| Matrix::row(r)).collect();
            let a = Matrix::vstack(&rows);
            let t = Vector::new(&promote_vec(&shift));

            Some(AffineMap::new(&a, &t))
        }
    }
}


fn promote_vec<T: Coord>(v: &Vec<Ratio<i32>>) -> Vec<T> {
    v.iter().map(|x| Coord::promote(*x)).collect()
}


impl<CSIn: Clone, CSOut: Clone>  CoordinateMap<BigRational, CSIn, CSOut> {
    pub fn from_string(s: &str) -> Option<Self> {
        Some(CoordinateMap::new(&AffineMap::from_string(s)?))
    }
}


pub fn parse_space_group_table<T: Read>(input: T) -> Result<(), Error> {
    let mut table: HashMap<String, TableEntry> = HashMap::new();
    let mut alias = HashMap::new();
    let mut lookup = vec![];
    let mut canonical_name = None;
    let mut current_name = None;

    for line in BufReader::new(input).lines() {
        let line = line?;
        let content = line.trim();

        if content.is_empty() || content.starts_with('#') {
            continue;
        } else if line.starts_with(' ') {
            let entry = table.get_mut(&current_name.clone().unwrap()).unwrap();
            entry.operators.push(AffineMap::from_string(content).unwrap());
        } else {
            let fields: Vec<_> = content.split_whitespace().collect();

            if fields[0].to_lowercase() == "alias" {
                let key = fields[1].to_string();
                let val = fields[2].to_string();
                alias.insert(key, val);
            } else if fields[0].to_lowercase() == "lookup" {
                let name = fields[1].to_string();
                let op = fields[4..].join(" ");
                let from_std = CoordinateMap::from_string(&op).unwrap();

                if from_std.dim() == 2 {
                    lookup.push(Lookup::Entry2d {
                        name,
                        system: CrystalSystem2d::from_string(fields[2]).unwrap(),
                        centering: Centering2d::from_string(fields[3]).unwrap(),
                        from_std,
                    });
                } else if from_std.dim() == 3 {
                    lookup.push(Lookup::Entry3d {
                        name,
                        system: CrystalSystem3d::from_string(fields[2]).unwrap(),
                        centering: Centering3d::from_string(fields[3]).unwrap(),
                        from_std,
                    });
                }
            } else {
                let op = fields[1..].join(" ");
                let transform: CoordinateMap<BigRational, Group, Setting> =
                    CoordinateMap::from_string(&op).unwrap();

                current_name = Some(fields[0].to_string());
                if transform == CoordinateMap::new(
                    &AffineMap::identity(transform.dim())
                ) {
                    canonical_name = current_name.clone();
                }

                table.insert(
                    current_name.clone().unwrap(),
                    TableEntry{
                        name: current_name.clone().unwrap(),
                        canonical_name: canonical_name.clone().unwrap(),
                        transform,
                        operators: vec![],
                    }
                );
            }
        }
    }

    Ok(())
}
