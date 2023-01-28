use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use std::ops::Neg;


pub trait LabelVector: Neg<Output = Self> + Copy + Eq + Hash {
    fn dim() -> usize;
    fn zero() -> Self;
    fn is_negative(&self) -> bool;
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct LabelVector2d {
    x: i16,
    y: i16,
}

impl LabelVector2d {
    pub fn new(x: i16, y: i16) -> Self {
        Self { x, y }
    }
}

impl LabelVector for LabelVector2d {
    fn dim() -> usize {
        2
    }

    fn zero() -> Self {
        Self { x: 0, y: 0 }
    }

    fn is_negative(&self) -> bool {
        self.x < 0 ||
        self.x == 0 && self.y < 0
    }
}

impl Neg for LabelVector2d {
    type Output = Self;

    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y }
    }
}

impl Display for LabelVector2d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {})", self.x, self.y))
    }
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct LabelVector3d {
    x: i16,
    y: i16,
    z: i16,
}

impl LabelVector3d {
    pub fn new(x: i16, y: i16, z: i16) -> Self {
        Self { x, y, z }
    }
}

impl LabelVector for LabelVector3d {
    fn dim() -> usize {
        3
    }

    fn zero() -> Self {
        Self { x: 0, y: 0,z: 0 }
    }

    fn is_negative(&self) -> bool {
        self.x < 0 ||
        self.x == 0 && self.y < 0 ||
        self.x == 0 && self.y == 0 && self.z < 0
    }
}

impl Neg for LabelVector3d {
    type Output = Self;

    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, z: -self.z }
    }
}

impl Display for LabelVector3d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {}, {})", self.x, self.y, self.z))
    }
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct VectorLabelledEdge<T>
{
    pub head: usize,
    pub tail: usize,
    pub shift: T,
}

impl<T> VectorLabelledEdge<T>
{
    pub fn new(head: usize, tail: usize, shift: T) -> Self {
        Self { head, tail, shift }
    }
}

impl<T> Neg for VectorLabelledEdge<T>
    where T: LabelVector
{
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            head: self.tail,
            tail: self.head,
            shift: -self.shift,
        }
    }

}

impl<T> VectorLabelledEdge<T>
    where T: LabelVector
{
    pub fn canonical(self) -> Self {
        if self.tail < self.head ||
            (self.tail == self.head && self.shift.is_negative())
        {
            -self
        } else {
            self
        }
    }
}

impl<T> Display for VectorLabelledEdge<T>
    where T: Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "[{} --{}-> {}]", self.head, self.shift, self.tail
        ))
    }
}


#[derive(Debug)]
pub struct Graph<T> {
    edges: Vec<VectorLabelledEdge<T>>,
}

impl<T> Graph<T>
    where T: LabelVector
{
    pub fn dim() -> usize { T::dim() }

    pub fn new(raw_edges: &[VectorLabelledEdge<T>]) -> Self {
        let mut edges = vec![];
        let mut seen = HashSet::new();

        for e in raw_edges {
            let e = e.canonical();

            if !seen.contains(&e) {
                seen.insert(e);
                edges.push(e);
            }
        }

        Graph { edges }
    }
}

impl<T> Display for Graph<T> 
    where T: Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut text = String::new();

        for e in &self.edges {
            text.push_str(&format!("{}", e)[..]);
            text.push('\n');
        }

        f.write_str(&text[..])
    }
}
