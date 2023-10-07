use std::collections::HashMap;
use std::fmt::Display;
use std::io::{BufRead, Read, BufReader};
use std::iter;

use num_rational::{Ratio, BigRational};
use num_traits::{Zero, Signed, One};

use crate::arithmetic::geometry::{AffineMap, Vector, CoordinateMap};
use crate::arithmetic::matrices::Matrix;

use super::parse_operator::parse_operator;
use super::types::{ Coord, CrystalSystem, Centering };


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


pub struct Lookup {
    name: String,
    system: CrystalSystem,
    centering: Centering,
    from_std: CoordinateMap<BigRational, Group, Setting>,
}


impl Display for Lookup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Lookup name: {}", self.name)?;
        writeln!(f, "Crystal system: {}", self.system)?;
        writeln!(f, "Centering: {}", self.centering)?;
        write!(f, "Transform: ")?;
        write_op(f, &self.from_std.forward())?;
        writeln!(f)?;
        Ok(())
    }
}


pub struct Tables {
    pub settings: HashMap<String, TableEntry>,
    pub alias: HashMap<String, String>,
    pub lookup: Vec<Lookup>,
}


impl Tables {
    pub fn lookup_settings(&self, s: CrystalSystem, c: Centering)
        -> Vec<(String, CoordinateMap<BigRational, Group, Setting>)>
    {
        self.lookup.iter().filter_map(|lkp|
            if s == lkp.system && c == lkp.centering {
                Some((lkp.name.clone(), lkp.from_std.clone()))
            } else {
                None
            }
        ).collect()
    }
}


impl CrystalSystem {
    pub fn from_string(s: &str) -> Option<Self> {
        match &s.to_lowercase()[..] {
            "square" => Some(CrystalSystem::Square),
            "rectangular" => Some(CrystalSystem::Rectangular),
            "oblique" => Some(CrystalSystem::Oblique),
            "cubic" => Some(CrystalSystem::Cubic),
            "orthorhombic" => Some(CrystalSystem::Orthorhombic),
            "hexagonal" => Some(CrystalSystem::Hexagonal),
            "tetragonal" => Some(CrystalSystem::Tetragonal),
            "trigonal" => Some(CrystalSystem::Trigonal),
            "monoclinic" => Some(CrystalSystem::Monoclinic),
            "triclinic" => Some(CrystalSystem::Triclinic),
            _ => None,
        }
    }
}


impl Centering {
    pub fn from_string(s: &str) -> Option<Self> {
        match s {
            "p" => Some(Centering::Primitive2d),
            "c" => Some(Centering::Centered2d),
            "P" => Some(Centering::Primitive),
            "F" => Some(Centering::FaceCentered),
            "I" => Some(Centering::BodyCentered),
            "R" => Some(Centering::Rhombohedral),
            "A" => Some(Centering::AFaceCentered),
            "B" => Some(Centering::BFaceCentered),
            "C" => Some(Centering::CFaceCentered),
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
        let line = line.ok()?;
        let content = line.trim();

        if content.is_empty() || content.starts_with('#') {
            continue;
        } else if line.starts_with(' ') {
            let op = AffineMap::from_string(content)?;
            settings.get_mut(&current_name)?.operators.push(op);
        } else {
            let fields: Vec<_> = content.split_whitespace().collect();

            if fields[0].to_lowercase() == "alias" {
                alias.insert(fields[1].into(), fields[2].into());
            } else if fields[0].to_lowercase() == "lookup" {
                lookup.push(make_lookup_entry_2d(&fields)?);
            } else {
                let name = fields[0].to_string();
                current_name = name.clone();

                let spec = fields[1..].join(" ");
                let transform = CoordinateMap::from_string(&spec)?;
                if transform.is_identity() {
                    canonical_name = name.clone();
                }

                let canonical_name = canonical_name.clone(); // local copy
                let operators = vec![];
                settings.insert(
                    current_name.clone(),
                    TableEntry{ name, canonical_name, transform, operators }
                );
            }
        }
    }

    Some(Tables { settings, alias, lookup })
}


fn make_lookup_entry_2d(fields: &Vec<&str>) -> Option<Lookup> {
    let name = fields[1].to_string();
    let spec = fields[4..].join(" ");
    let from_std = CoordinateMap::from_string(&spec)?;

    let system = CrystalSystem::from_string(fields[2])?;
    let centering = Centering::from_string(fields[3])?;
    Some(Lookup { name, system, centering, from_std })
}
