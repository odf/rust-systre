use std::cell::UnsafeCell;
use std::collections::{BTreeSet, BTreeMap, HashSet, HashMap, VecDeque};
use std::fmt::Display;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;
use std::mem::replace;
use std::ops::{Neg, Add, Sub, Mul};
use itertools::Itertools;
use num_bigint::BigInt;
use num_rational::BigRational;

use crate::arithmetic::geometry;
use crate::arithmetic::linear_algebra::extend_basis;
use crate::partitions::Partition;
use crate::symmetries::symmetries;


pub trait LabelVector<CS>:
    Copy + Eq + Hash + Ord +
    Neg<Output = Self> +
    Add<Self, Output = Self> +
    Sub<Self, Output = Self>
{
    fn dim() -> usize;
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn is_negative(&self) -> bool;
    fn is_positive(&self) -> bool;
    fn to_vec(&self) -> Vec<i32>;
    fn from_vec(v: &[i32]) -> Self;
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct LabelVector2d<CS> {
    x: i32,
    y: i32,
    phantom: PhantomData<CS>,
}

impl<CS> LabelVector2d<CS> {
    pub fn new(x: i32, y: i32) -> Self {
        Self { x, y, phantom: PhantomData::default() }
    }
}

impl<CS> LabelVector<CS> for LabelVector2d<CS>
    where CS: Copy + Ord + Hash
{
    fn dim() -> usize {
        2
    }

    fn zero() -> Self {
        Self { x: 0, y: 0, phantom: PhantomData::default() }
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

    fn to_vec(&self) -> Vec<i32> {
        vec![self.x, self.y]
    }

    fn from_vec(v: &[i32]) -> Self {
        assert_eq!(v.len(), 2);
        Self::new(v[0], v[1])
    }
}

impl<CS> Neg for LabelVector2d<CS> {
    type Output = Self;

    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, phantom: self.phantom }
    }
}

impl<CS> Add<Self> for LabelVector2d<CS> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self { x: self.x + rhs.x, y: self.y + rhs.y, phantom: self.phantom }
    }
}

impl<CS> Sub<Self> for LabelVector2d<CS> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self { x: self.x - rhs.x, y: self.y - rhs.y, phantom: self.phantom }
    }
}

impl<CS> Display for LabelVector2d<CS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {})", self.x, self.y))
    }
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct LabelVector3d<CS> {
    x: i32,
    y: i32,
    z: i32,
    phantom: PhantomData<CS>,
}

impl<CS> LabelVector3d<CS> {
    pub fn new(x: i32, y: i32, z: i32) -> Self {
        Self { x, y, z, phantom: PhantomData::default() }
    }
}

impl<CS> LabelVector<CS> for LabelVector3d<CS>
    where CS: Copy + Ord + Hash
{
    fn dim() -> usize {
        3
    }

    fn zero() -> Self {
        Self { x: 0, y: 0,z: 0, phantom: PhantomData::default() }
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

    fn to_vec(&self) -> Vec<i32> {
        vec![self.x, self.y, self.z]
    }

    fn from_vec(v: &[i32]) -> Self {
        assert_eq!(v.len(), 3);
        Self::new(v[0], v[1], v[2])
    }
}

impl<CS> Neg for LabelVector3d<CS> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            phantom: self.phantom,
        }
    }
}

impl<CS> Add<Self> for LabelVector3d<CS> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            phantom: self.phantom,
        }
    }
}

impl<CS> Sub<Self> for LabelVector3d<CS> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            phantom: self.phantom,
        }
    }
}

impl<CS> Display for LabelVector3d<CS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}, {}, {})", self.x, self.y, self.z))
    }
}


pub type Vertex = u32;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub struct InputCS {}

pub type Point = geometry::Point<BigRational, InputCS>;
pub type Vector = geometry::Vector<BigRational, InputCS>;
pub type AffineMap = geometry::AffineMap<BigRational, InputCS>;


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Edge<T, CS>
{
    pub head: Vertex,
    pub tail: Vertex,
    pub shift: T,
    phantom: PhantomData<CS>,
}

impl<T, CS> Edge<T, CS>
{
    pub fn new(head: Vertex, tail: Vertex, shift: T) -> Self {
        Self { head, tail, shift, phantom: PhantomData::default() }
    }
}

impl<T, CS> Neg for &Edge<T, CS>
    where T: LabelVector<CS>
{
    type Output = Edge<T, CS>;

    fn neg(self) -> Self::Output {
        Self::Output {
            head: self.tail,
            tail: self.head,
            shift: -self.shift,
            phantom: PhantomData::default()
        }
    }
}

impl<T, CS> Edge<T, CS>
    where
        T: LabelVector<CS>,
        CS: Clone
{
    pub fn canonical(&self) -> Self {
        if self.tail < self.head ||
            (self.tail == self.head && self.shift.is_negative())
        {
            -self
        } else {
            self.clone()
        }
    }
}

impl<T, CS> Display for Edge<T, CS>
    where T: Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "{} --{}-> {}", self.head, self.shift, self.tail
        ))
    }
}


#[derive(Debug)]
pub struct Graph<T, CS>
    where
        T: Eq + Hash,
        CS: Eq + Hash
{
    edges: Vec<Edge<T, CS>>,
    vertices: UnsafeCell<Option<Vec<Vertex>>>,
    incidences: UnsafeCell<BTreeMap<Vertex, Vec<Edge<T, CS>>>>,
    positions: UnsafeCell<BTreeMap<Vertex, Point>>,
    edge_lookup: UnsafeCell<
        BTreeMap<Vertex, HashMap<Vector, Edge<T, CS>>>
    >,
    symmetries: UnsafeCell<Vec<Automorphism<T, CS>>>,
}

impl<T, CS> Graph<T, CS>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    pub fn dim() -> usize { T::dim() }

    pub fn new(raw_edges: &[Edge<T, CS>]) -> Self {
        let edges: Vec<Edge<T, CS>> = raw_edges.iter()
            .map(|e| e.canonical())
            .collect::<HashSet<_>>()
            .into_iter()
            .sorted()
            .collect();

        Graph {
            edges,
            vertices: UnsafeCell::new(None),
            incidences: UnsafeCell::new(BTreeMap::new()),
            positions: UnsafeCell::new(BTreeMap::new()),
            edge_lookup: UnsafeCell::new(BTreeMap::new()),
            symmetries: UnsafeCell::new(vec![]),
        }
    }

    pub fn vertices(&self) -> Vec<Vertex> {
        if let Some(vertices) = unsafe { (*self.vertices.get()).clone() } {
            vertices
        } else {
            let vertices: Vec<_> = self.edges.iter()
                .flat_map(|e| [e.head, e.tail])
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();

            unsafe { *self.vertices.get() = Some(vertices.clone()) };

            vertices
        }
    }

    pub fn directed_edges(&self) -> Vec<Edge<T, CS>>
    {
        self.vertices().iter().flat_map(|v| self.incidences(v)).collect()
    }

    pub fn incidences(&self, v: &Vertex) -> Vec<Edge<T, CS>> {
        if let Some(output) = unsafe { (*self.incidences.get()).get(v) } {
            output.clone()
        } else {
            let mut incidences = BTreeMap::new();

            for e in &self.edges {
                for e in [e.clone(), -e] {
                    if !incidences.contains_key(&e.head) {
                        incidences.insert(e.head, vec![]);
                    }
                    incidences.get_mut(&e.head).unwrap().push(
                        Edge::new(e.head, e.tail, e.shift)
                    );
                }
            }

            let output = incidences[v].clone();

            unsafe { *self.incidences.get() = incidences };

            output
        }
    }

    pub fn position(&self, v: &Vertex) -> Point {
        if let Some(output) = unsafe { (*self.positions.get()).get(v) } {
            output.clone()
        } else {
            let positions = barycentric_placement(self);
            let output = positions[v].clone();
            unsafe { *self.positions.get() = positions };
            output
        }
    }

    pub fn edge_by_unique_delta(&self, v: &Vertex, delta: &Vector)
        -> Option<Edge<T, CS>>
    {
        if (unsafe { (*self.edge_lookup.get()).get(v) }).is_none() {
            let data = edges_by_unique_deltas(self);
            unsafe { *self.edge_lookup.get() = data };
        }

        let lookup = unsafe { self.edge_lookup.get().as_ref() };
        Some(lookup?.get(v)?.get(delta)?.clone())
    }

    pub fn symmetries(&self) -> Vec<Automorphism<T, CS>> {
        let syms = unsafe { self.symmetries.get().as_ref().unwrap() };
        if syms.is_empty() {
            let syms = symmetries(&self);
            unsafe { *self.symmetries.get() = syms.clone() }
            syms
        } else {
            syms.clone()
        }
    }

    pub fn position_normalized(&self, v: &Vertex) -> Point {
        self.position(&v).iter().map(|q| q - q.floor()).collect()
    }

    pub fn edge_vector(&self, e: &Edge<T, CS>)
        -> Vector
    {
        let d = T::dim();
        let p = self.position(&e.head);
        let q = self.position(&e.tail);
        let s: Vec<_> = e.shift.to_vec();

        (0..d).map(|i| &q[i] + BigInt::from(s[i]) - &p[i]).collect()
    }

    pub fn is_stable(&self) -> bool {
        let mut seen = HashSet::new();

        for v in self.vertices() {
            let p = self.position_normalized(&v);
            if seen.contains(&p) {
                return false;
            } else {
                seen.insert(p);
            }
        }
        true
    }

    pub fn is_locally_stable(&self) -> bool {
        for v in self.vertices() {
            let mut seen = HashSet::new();

            for ngb in self.incidences(&v) {
                let q = self.edge_vector(&ngb);

                if seen.contains(&q) {
                    return false;
                } else {
                    seen.insert(q);
                }
            }
        }
        true
    }

    pub fn has_second_order_collisions(&self) -> bool {
        let mut seen = HashSet::new();

        for v in self.vertices() {
            let mut deltas: Vec<Vec<_>> = vec![];
            for e in self.incidences(&v) {
                let v = self.edge_vector(&e);
                deltas.push(v.into_iter().collect());
            }
            deltas.sort();
            let p = self.position_normalized(&v);
            deltas.push(p.into_iter().collect());

            if seen.contains(&deltas) {
                return true;
            } else {
                seen.insert(deltas);
            }
        }
        false
    }

    pub fn coordination_sequence(&self, v: &Vertex)
        -> CoordinationSequence<T, CS>
    {
        CoordinationSequence::new(self, v)
    }

    pub fn shift_normalized(&self) -> Self {
        let mut edges = vec![];
        let mut seen = HashSet::new();

        for v in self.vertices() {
            if !seen.contains(&v) {
                for e in traverse_with_shift_adjustments(&self, &v) {
                    seen.insert(e.head);
                    seen.insert(e.tail);
                    edges.push(e);
                }
            }
        }

        Self::new(&edges)
    }

    pub fn is_connected(&self) -> bool {
        if let Some(v0) = self.vertices().first() {
            let (size, _, multiplicity) = graph_component_measures(self, v0);
            size == self.vertices().len() && multiplicity == Some(1)
        } else {
            true
        }
    }

    pub fn vertex_orbits(&self) -> Vec<Vec<Vertex>> {
        let mut p = Partition::new();

        for phi in self.symmetries() {
            for v in self.vertices() {
                p.unite(&v, &phi.vertex_map.get(&v).unwrap())
            }
        }

        p.classes(&self.vertices())
    }
}

impl<T, CS> Display for Graph<T, CS> 
    where
        T: Display + Eq + Hash,
        CS: Eq + Hash
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


pub struct CoordinationSequence<'a, T, CS>
    where
        T: Eq + Hash,
        CS: Eq + Hash
{
    graph: &'a Graph<T, CS>,
    last_shell: HashSet<(Vertex, T)>,
    this_shell: HashSet<(Vertex, T)>,
}

impl <'a, T, CS> CoordinationSequence<'a, T, CS>
    where
        T: LabelVector<CS> + Eq + Hash,
        CS: Eq + Hash
{
    fn new(graph: &'a Graph<T, CS>, v0: &Vertex) -> Self {
        let last_shell = HashSet::new();
        let this_shell = HashSet::from([(*v0, T::zero())]);

        CoordinationSequence { graph, last_shell, this_shell }
    }
}

impl<'a, T, CS> Iterator for CoordinationSequence<'a, T, CS>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_shell = HashSet::new();

        for (vertex, shift) in &self.this_shell {
            for e in self.graph.incidences(&vertex) {
                let w = (e.tail, e.shift + *shift);

                if !self.last_shell.contains(&w) &&
                    !self.this_shell.contains(&w)
                {
                    next_shell.insert(w);
                }
            }
        }

        self.last_shell = replace(&mut self.this_shell, next_shell);
        Some(self.this_shell.len())
    }
}


#[derive(Clone, Eq, PartialEq)]
pub struct Automorphism<T, CS>
    where
        T: Eq + Hash,
        CS: Eq + Hash
{
    pub(crate) vertex_map: HashMap<Vertex, Vertex>,
    pub(crate) edge_map: HashMap<Edge<T, CS>, Edge<T, CS>>,
    pub(crate) transform: AffineMap,
}


impl<T, CS> Automorphism<T, CS>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    pub fn new(
        vertex_map: HashMap<Vertex, Vertex>,
        edge_map: HashMap<Edge<T, CS>, Edge<T, CS>>,
        transform: AffineMap,
    ) -> Automorphism<T, CS>
    {
        Automorphism { vertex_map, edge_map, transform }
    }

    pub fn identity(graph: &Graph<T, CS>) -> Automorphism<T, CS> {
        let transform = AffineMap::identity(T::dim());

        let vertex_map = graph.vertices().iter()
            .map(|v| (*v, *v))
            .collect();

        let edge_map = graph.vertices().iter()
            .flat_map(|v| graph.incidences(v))
            .map(|e| (e.clone(), e))
            .collect();

        Automorphism { vertex_map, edge_map, transform }
    }
}


impl<T, CS> Mul<&Automorphism<T, CS>> for &Automorphism<T, CS>
    where
        T: Clone + Eq + Hash,
        CS: Clone + Eq + Hash,
        for <'a> &'a AffineMap: Mul<&'a AffineMap, Output=AffineMap>
{
    type Output = Automorphism<T, CS>;

    fn mul(self, rhs: &Automorphism<T, CS>) -> Self::Output {
        let transform = &self.transform * &rhs.transform;

        let vertex_map = self.vertex_map.keys().map(|v| (
            *v,
            *self.vertex_map.get(rhs.vertex_map.get(v).unwrap()).unwrap()
        )).collect();

        let edge_map = self.edge_map.keys().map(|e| (
            e.clone(),
            self.edge_map.get(rhs.edge_map.get(e).unwrap()).unwrap().clone()
        )).collect();

        Automorphism { vertex_map, edge_map, transform }
    }
}


impl<T, CS> Hash for Automorphism<T, CS>
    where
        T: Eq + Hash,
        CS: Eq + Hash
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.transform.hash(state);
        self.vertex_map.iter().sorted().collect::<Vec<_>>().hash(state)
    }
}


fn traverse_with_shift_adjustments<T, CS>(g: &Graph<T, CS>, v0: &Vertex)
    -> Vec<Edge<T, CS>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut shifts = BTreeMap::from([(*v0, T::zero())]);
    let mut seen = HashSet::new();
    let mut queue = VecDeque::from([*v0]);
    let mut result = Vec::new();

    while let Some(v) = queue.pop_front() {
        for e in g.incidences(&v) {
            if !shifts.contains_key(&e.tail) {
                shifts.insert(e.tail, shifts[&v] + e.shift);
                queue.push_back(e.tail);
            }

            let e = e.canonical();
            if !seen.contains(&e) {
                let shift = shifts[&e.head] + e.shift - shifts[&e.tail];
                result.push(Edge::new(e.head, e.tail, shift));
                seen.insert(e);
            }
        }
    }

    result
}


pub fn graph_component_measures<T, CS>(g: &Graph<T, CS>, v0: &Vertex)
    -> (usize, usize, Option<i32>) // TODO return a struct?
    where
        T: LabelVector<CS>,
        CS: Clone + Hash + Ord
{
    let edges = traverse_with_shift_adjustments(g, v0);

    let mut vertices = BTreeSet::new();
    for e in &edges {
        vertices.insert(e.head);
        vertices.insert(e.tail);
    }

    let mut basis = vec![];
    for e in &edges {
        extend_basis(&e.shift.to_vec(), &mut basis);
    }

    let size = vertices.len();
    let rank = basis.len();

    let multiplicity = if rank == T::dim() {
        let mut d: i32 = 1;
        for i in 0..basis.len() {
            d *= basis[i][i];
        }
        Some(d.abs())
    } else {
        None
    };

    (size, rank, multiplicity)
}


fn barycentric_placement<T, CS>(g: &Graph<T, CS>)
    -> BTreeMap<Vertex, Point>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let verts = g.vertices();
    let vidcs: BTreeMap<_, _> =
        verts.iter().enumerate().map(|(i, &e)| (e, i)).collect();

    let n = verts.len();
    let d = T::dim();

    let mut a = vec![vec![0 as i64; n]; n];
    let mut t = vec![vec![0 as i64; d]; n];
    a[0][0] = 1;

    for i in 1..n {
        for ngb in g.incidences(&verts[i]) {
            let j = vidcs[&ngb.tail];
            a[i][j] -= 1;
            a[i][i] += 1;

            let s = ngb.shift.to_vec();
            for k in 0..d {
                t[i][k] += s[k] as i64;
            }
        }
    }

    let p = crate::arithmetic::modular_solver::solve(&a, &t).unwrap();

    let mut result = BTreeMap::new();
    for i in 0..n {
        result.insert(verts[i], Point::new(&p[i]));
    }

    result
}


fn edges_by_unique_deltas<T, CS>(g: &Graph<T, CS>)
    -> BTreeMap<Vertex, HashMap<Vector, Edge<T, CS>>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut result = BTreeMap::new();

    for head in g.vertices() {
        let mut seen = HashSet::new();
        let mut to_edge = HashMap::new();

        for e in g.incidences(&head) {
            let delta = g.edge_vector(&e);
            if seen.contains(&delta) {
                to_edge.remove(&delta);
            } else {
                to_edge.insert(delta.clone(), e);
                seen.insert(delta);
            }
        }

        result.insert(head, to_edge);
    }

    result
}
