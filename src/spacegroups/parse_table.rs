use std::collections::HashMap;
use std::fmt::Display;
use std::io::{BufRead, Read, BufReader};
use std::iter;

use num_rational::{Ratio, BigRational};
use num_traits::{Zero, Signed, One};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::matrices::Matrix;

use super::parse_operator::parse_operator;
use super::types::{
    Coord, CrystalSystem2d, Centering2d, CrystalSystem3d, Centering3d
};


#[derive(Clone, Debug, PartialEq)]
pub struct Group {}


#[derive(Clone, Debug, PartialEq)]
pub struct Setting {}


pub struct TableEntry {
    name: String,
    canonical_name: String,
    transform: CoordinateMap<BigRational, Group, Setting>,
    operators: Vec<AffineMap<BigRational, Setting>>,
}


impl Display for TableEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Setting Name: {}", self.name)?;
        writeln!(f, "Group name: {}", self.canonical_name)?;
        write!(f, "Transform: ")?;
        write_op(f, &self.transform.forward())?;
        writeln!(f, "")?;
        writeln!(f, "Operators:")?;
        for op in &self.operators {
            write_op(f, &op)?;
            writeln!(f, "")?;
        }
        Ok(())
    }
}


fn write_op<CS>(
    f: &mut std::fmt::Formatter<'_>,
    op: &AffineMap<BigRational, CS>
)
    -> std::fmt::Result
    where CS: Clone
{
    if op.dim() > 3 {
        writeln!(f, "{}", op)?;
    } else {
        for i in 0..op.dim() {
            if i > 0 {
                write!(f, ",")?;
            }
            write_op_row(f, op.linear_matrix().get_row(i), &op.shift()[i])?;
        }
    }

    Ok(())
}


fn write_op_row(
    f: &mut std::fmt::Formatter<'_>,
    row: Vec<BigRational>,
    scalar: &BigRational
) -> Result<(), std::fmt::Error>
{
    let parts: Vec<_> = (0..row.len())
        .map(|j| (&row[j], ["x", "y", "z"][j]))
        .chain(iter::once((scalar, "")))
        .filter(|(val, _)| !val.is_zero())
        .collect();

    if parts.is_empty() {
        write!(f, "0")?;
    } else {
        for (k, (val, axis)) in parts.iter().enumerate() {
            if k > 0 && val.is_positive() {
                write!(f, "+")?;
            } else if val.is_negative() {
                write!(f, "-")?;
            }
            if !val.abs().is_one() {
                write!(f, "{}", val.abs())?;
            }
            write!(f, "{}", axis)?;
        }
    }

    Ok(())
}


pub enum Lookup {
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


impl Display for Lookup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Lookup::Entry2d { name, system, centering, from_std } => {
                writeln!(f, "Lookup name: {}", name)?;
                writeln!(f, "Crystal system: {}", system)?;
                writeln!(f, "Centering: {}", centering)?;
                write!(f, "Transform: ")?;
                write_op(f, &from_std.forward())?;
                writeln!(f)?;
            },
            Lookup::Entry3d { name, system, centering, from_std } => {
                writeln!(f, "Lookup name: {}", name)?;
                writeln!(f, "Crystal system: {}", system)?;
                writeln!(f, "Centering: {}", centering)?;
                write!(f, "Transform: ")?;
                write_op(f, &from_std.forward())?;
                writeln!(f)?;
            },
        }
        Ok(())
    }
}


pub struct Tables {
    pub settings: HashMap<String, TableEntry>,
    pub alias: HashMap<String, String>,
    pub lookup: Vec<Lookup>,
}


impl CrystalSystem2d {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            "square" => Some(CrystalSystem2d::Square),
            "hexagonal" => Some(CrystalSystem2d::Hexagonal),
            "rectangular" => Some(CrystalSystem2d::Rectangular),
            "oblique" => Some(CrystalSystem2d::Oblique),
            _ => None,
        }
    }
}


impl Centering2d {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            "p" => Some(Centering2d::Primitive),
            "c" => Some(Centering2d::Centered),
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


pub fn parse_space_group_table<T: Read>(input: T) -> Option<Tables> {
    let mut settings: HashMap<String, TableEntry> = HashMap::new();
    let mut alias = HashMap::new();
    let mut lookup = vec![];
    let mut canonical_name = "".to_string();
    let mut current_name = "".to_string();

    for line in BufReader::new(input).lines() {
        let line = line.ok().unwrap();
        let content = line.trim();

        if content.is_empty() || content.starts_with('#') {
            continue;
        } else if line.starts_with(' ') {
            let entry = settings.get_mut(&current_name).unwrap();
            entry.operators.push(AffineMap::from_string(content).unwrap());
        } else {
            let fields: Vec<_> = content.split_whitespace().collect();

            if fields[0].to_lowercase() == "alias" {
                let key = fields[1].to_string();
                let val = fields[2].to_string();
                alias.insert(key, val);
            } else if fields[0].to_lowercase() == "lookup" {
                lookup.push(make_lookup_entry(&fields));
            } else {
                let name = fields[0].to_string();
                current_name = name.clone();

                let op = fields[1..].join(" ");
                let transform = CoordinateMap::from_string(&op).unwrap();
                if transform == CoordinateMap::new(
                    &AffineMap::identity(transform.dim())
                ) {
                    canonical_name = name.clone();
                }

                let canonical_name = canonical_name.clone(); // local copy
                let operators = vec![];
                let entry = TableEntry{ name, canonical_name, transform, operators };
                settings.insert(current_name.clone(), entry);
            }
        }
    }

    Some(Tables { settings, alias, lookup })
}


fn make_lookup_entry(fields: &Vec<&str>) -> Lookup {
    let name = fields[1].to_string();
    let op = fields[4..].join(" ");
    let from_std = CoordinateMap::from_string(&op).unwrap();

    if from_std.dim() == 2 {
        Lookup::Entry2d {
            name,
            system: CrystalSystem2d::from_string(fields[2]).unwrap(),
            centering: Centering2d::from_string(fields[3]).unwrap(),
            from_std,
        }
    } else {
        Lookup::Entry3d {
            name,
            system: CrystalSystem3d::from_string(fields[2]).unwrap(),
            centering: Centering3d::from_string(fields[3]).unwrap(),
            from_std,
        }
    }
}
