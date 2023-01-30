use std::cell::UnsafeCell;
use std::collections::{HashSet, BTreeSet, BTreeMap};
use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Neg, Add, Sub};


pub trait LabelVector:
    Copy + Eq + Hash +
    Neg<Output = Self> +
    Add<Self, Output = Self> +
    Sub<Self, Output = Self>
{
    fn dim() -> u8;
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn is_negative(&self) -> bool;
    fn is_positive(&self) -> bool;
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
    fn dim() -> u8 {
        2
    }

    fn zero() -> Self {
        Self { x: 0, y: 0 }
    }

    fn is_zero(&self) -> bool {
        self.x == 0 && self.y == 0
    }

    fn is_negative(&self) -> bool {
        self.x < 0 ||
        self.x == 0 && self.y < 0
    }

    fn is_positive(&self) -> bool {
        self.x > 0 ||
        self.x == 0 && self.y > 0
    }
}

impl Neg for LabelVector2d {
    type Output = Self;

    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y }
    }
}

impl Add<Self> for LabelVector2d {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self { x: self.x + rhs.x, y: self.y + rhs.y }
    }
}

impl Sub<Self> for LabelVector2d {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self { x: self.x - rhs.x, y: self.y - rhs.y }
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
    fn dim() -> u8 {
        3
    }

    fn zero() -> Self {
        Self { x: 0, y: 0,z: 0 }
    }

    fn is_zero(&self) -> bool {
        self.x == 0 && self.y == 0 && self.z == 0
    }

    fn is_negative(&self) -> bool {
        self.x < 0 ||
        self.x == 0 && self.y < 0 ||
        self.x == 0 && self.y == 0 && self.z < 0
    }

    fn is_positive(&self) -> bool {
        self.x > 0 ||
        self.x == 0 && self.y > 0 ||
        self.x == 0 && self.y == 0 && self.z > 0
    }
}

impl Neg for LabelVector3d {
    type Output = Self;

    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, z: -self.z }
    }
}

impl Add<Self> for LabelVector3d {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
    }
}

impl Sub<Self> for LabelVector3d {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

impl Display for LabelVector3d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {}, {})", self.x, self.y, self.z))
    }
}


type Vertex = u32;


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct ShiftedVertex<T>
{
    pub vertex: Vertex,
    pub shift: T,
}

impl<T> ShiftedVertex<T>
{
    pub fn new(vertex: Vertex, shift: T) -> Self {
        Self { vertex, shift }
    }
}

impl<T> Display for ShiftedVertex<T>
    where T: Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{} {}", self.vertex, self.shift))
    }
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct VectorLabelledEdge<T>
{
    pub head: Vertex,
    pub tail: Vertex,
    pub shift: T,
}

impl<T> VectorLabelledEdge<T>
{
    pub fn new(head: Vertex, tail: Vertex, shift: T) -> Self {
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
            "{} --{}-> {}", self.head, self.shift, self.tail
        ))
    }
}


#[derive(Debug)]
pub struct Graph<T> {
    edges: Vec<VectorLabelledEdge<T>>,
    vertices: UnsafeCell<Option<Vec<Vertex>>>,
    incidences: UnsafeCell<BTreeMap<Vertex, Vec<ShiftedVertex<T>>>>,
}

impl<T> Graph<T>
    where T: LabelVector
{
    pub fn dim() -> u8 { T::dim() }

    pub fn new(raw_edges: &[VectorLabelledEdge<T>]) -> Self {
        let edges = raw_edges.iter().map(|e| e.canonical())
            .collect::<HashSet<_>>()
            .iter().cloned()
            .collect();

        Graph {
            edges,
            vertices: UnsafeCell::new(None),
            incidences: UnsafeCell::new(BTreeMap::new()),
        }
    }

    pub fn vertices(&self) -> Vec<Vertex> {
        if let Some(vertices) = unsafe { (*self.vertices.get()).clone() } {
            vertices
        } else {
            let vertices: Vec<_> = self.edges.iter()
                .flat_map(|e| [e.head, e.tail])
                .collect::<BTreeSet<_>>()
                .iter().cloned()
                .collect();

            unsafe { *self.vertices.get() = Some(vertices.clone()) };

            vertices
        }
    }

    pub fn incidences(&self, v: &Vertex) -> Vec<ShiftedVertex<T>> {
        if let Some(output) = unsafe { (*self.incidences.get()).get(v) } {
            output.clone()
        } else {
            let mut incidences = BTreeMap::new();

            for &e in &self.edges {
                for e in [e, -e] {
                    if !incidences.contains_key(&e.head) {
                        incidences.insert(e.head, vec![]);
                    }
                    incidences.get_mut(&e.head).unwrap()
                        .push(ShiftedVertex::new(e.tail, e.shift));
                }
            }

            let output = incidences[v].clone();

            unsafe { *self.incidences.get() = incidences };

            output
        }
    }

    pub fn coordination_sequence(&self, v: &Vertex)
        -> CoordinationSequence<T>
    {
        CoordinationSequence::new(self, v)
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


pub struct CoordinationSequence<'a, T> {
    graph: &'a Graph<T>,
    last_shell: HashSet<ShiftedVertex<T>>,
    this_shell: HashSet<ShiftedVertex<T>>,
}

impl <'a, T> CoordinationSequence<'a, T>
    where T: LabelVector
{
    fn new(graph: &'a Graph<T>, v0: &Vertex) -> Self {
        let last_shell = HashSet::new();
        let this_shell = HashSet::from([ShiftedVertex::new(*v0, T::zero())]);

        CoordinationSequence { graph, last_shell, this_shell }
    }
}

impl<'a, T> Iterator for CoordinationSequence<'a, T>
    where T: LabelVector
{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_shell = HashSet::new();

        for ShiftedVertex { vertex: v0, shift: s0 } in &self.this_shell {
            for ShiftedVertex { vertex, shift } in self.graph.incidences(&v0) {
                let w = ShiftedVertex::new(vertex, shift + *s0);

                if !self.last_shell.contains(&w) &&
                    !self.this_shell.contains(&w)
                {
                    next_shell.insert(w);
                }
            }
        }

        let n = next_shell.len();

        self.last_shell = self.this_shell.clone();
        self.this_shell = next_shell;

        Some(n)
    }
}
