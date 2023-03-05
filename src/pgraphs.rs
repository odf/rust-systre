use std::cell::UnsafeCell;
use std::collections::{BTreeSet, BTreeMap, HashSet, HashMap, VecDeque};
use std::fmt::Display;
use std::hash::Hash;
use std::mem::replace;
use std::ops::{Neg, Mul};
use itertools::Itertools;
use num_bigint::BigInt;
use num_rational::BigRational;

use crate::arithmetic::geometry;
use crate::arithmetic::linear_algebra::extend_basis;


pub type Vertex = u32;

pub type Point<CS> = geometry::Point<BigRational, CS>;
pub type Vector<CS> = geometry::Vector<BigRational, CS>;
pub type Shift<CS> = geometry::Vector<i32, CS>;
pub type AffineMap<CS> = geometry::AffineMap<BigRational, CS>;


#[derive(Clone, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Edge<CS>
{
    pub head: Vertex,
    pub tail: Vertex,
    pub shift: Shift<CS>,
}

impl<CS> Edge<CS>
{
    pub fn new(head: Vertex, tail: Vertex, shift: Shift<CS>) -> Self {
        Self { head, tail, shift }
    }
}

impl<CS> Neg for Edge<CS> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            head: self.tail,
            tail: self.head,
            shift: -self.shift,
        }
    }
}

impl<CS> Edge<CS> {
    pub fn canonical(self) -> Self
        where CS: Clone
    {
        let negative_shift = match self.shift.iter().find(|&&x| x != 0) {
            Some(&x) if x < 0 => true,
            _ => false,
        };
        if self.tail < self.head || (self.tail == self.head && negative_shift) {
            -self
        } else {
            self
        }
    }
}

impl<CS> Display for Edge<CS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.shift.iter().join(", ");
        f.write_fmt(format_args!(
            "{} --({})-> {}", self.head, s, self.tail
        ))
    }
}


#[derive(Debug)]
pub struct Graph<CS> {
    dim: usize,
    edges: Vec<Edge<CS>>,
    vertices: UnsafeCell<Option<Vec<Vertex>>>,
    incidences: UnsafeCell<BTreeMap<Vertex, Vec<Edge<CS>>>>,
    positions: UnsafeCell<BTreeMap<Vertex, Point<CS>>>,
    edge_lookup: UnsafeCell<BTreeMap<Vertex, HashMap<Vector<CS>, Edge<CS>>>>,
}

impl<CS> Graph<CS>
    where CS: Eq + Hash + Ord + Clone
{
    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn new(raw_edges: &[Edge<CS>]) -> Self {
        assert!(raw_edges.len() > 0);

        let dim = raw_edges[0].shift.dim();

        let edges: Vec<Edge<CS>> = raw_edges.iter()
            .map(|e| e.clone().canonical())
            .collect::<HashSet<_>>()
            .into_iter()
            .sorted()
            .collect();

        Graph {
            dim,
            edges,
            vertices: UnsafeCell::new(None),
            incidences: UnsafeCell::new(BTreeMap::new()),
            positions: UnsafeCell::new(BTreeMap::new()),
            edge_lookup: UnsafeCell::new(BTreeMap::new()),
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

    pub fn directed_edges(&self) -> Vec<Edge<CS>>
    {
        self.vertices().iter().flat_map(|v| self.incidences(v)).collect()
    }

    pub fn incidences(&self, v: &Vertex) -> Vec<Edge<CS>> {
        if let Some(output) = unsafe { (*self.incidences.get()).get(v) } {
            output.clone()
        } else {
            let mut incidences = BTreeMap::new();

            for e in &self.edges {
                for e in [e.clone(), -e.clone()] {
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

    pub fn position(&self, v: &Vertex) -> Point<CS> {
        if let Some(output) = unsafe { (*self.positions.get()).get(v) } {
            output.clone()
        } else {
            let positions = barycentric_placement(self);
            let output = positions[v].clone();
            unsafe { *self.positions.get() = positions };
            output
        }
    }

    pub fn edge_by_unique_delta(&self, v: &Vertex, delta: &Vector<CS>)
        -> Option<Edge<CS>>
    {
        if (unsafe { (*self.edge_lookup.get()).get(v) }).is_none() {
            let data = edges_by_unique_deltas(self);
            unsafe { *self.edge_lookup.get() = data };
        }

        let lookup = unsafe { self.edge_lookup.get().as_ref() };
        Some(lookup?.get(v)?.get(delta)?.clone())
    }

    pub fn position_normalized(&self, v: &Vertex) -> Point<CS> {
        self.position(&v).iter().map(|q| q - q.floor()).collect()
    }

    pub fn edge_vector(&self, e: &Edge<CS>) -> Vector<CS> {
        let d = self.dim();
        let p = self.position(&e.head);
        let q = self.position(&e.tail);
        let s: Vec<_> = e.shift.iter().collect();

        (0..d).map(|i| &q[i] + BigInt::from(*s[i]) - &p[i]).collect()
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
        -> CoordinationSequence<CS>
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
}

impl<CS> Display for Graph<CS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut text = String::new();

        for e in &self.edges {
            text.push_str(&format!("{}", e)[..]);
            text.push('\n');
        }

        f.write_str(&text[..])
    }
}


pub struct CoordinationSequence<'a, CS> {
    graph: &'a Graph<CS>,
    last_shell: HashSet<(Vertex, Shift<CS>)>,
    this_shell: HashSet<(Vertex, Shift<CS>)>,
}

impl <'a, CS> CoordinationSequence<'a, CS>
    where CS: Eq + Hash + Ord + Clone
{
    fn new(graph: &'a Graph<CS>, v0: &Vertex) -> Self {
        let last_shell = HashSet::new();
        let this_shell = HashSet::from([(*v0, Shift::zero(graph.dim()))]);

        CoordinationSequence { graph, last_shell, this_shell }
    }
}

impl<'a, CS> Iterator for CoordinationSequence<'a, CS>
    where CS: Eq + Hash + Ord + Clone
{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_shell = HashSet::new();

        for (vertex, shift) in &self.this_shell {
            for e in self.graph.incidences(&vertex) {
                let w = (e.tail, e.shift + shift);

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


pub struct Automorphism<CS> {
    vertex_map: HashMap<Vertex, Vertex>,
    edge_map: HashMap<Edge<CS>, Edge<CS>>,
    transform: AffineMap<CS>,
}


impl<CS> Automorphism<CS>
    where CS: Clone + Eq + Hash + Ord
{
    pub fn new(
        vertex_map: HashMap<Vertex, Vertex>,
        edge_map: HashMap<Edge<CS>, Edge<CS>>,
        transform: AffineMap<CS>,
    ) -> Automorphism<CS>
    {
        Automorphism { vertex_map, edge_map, transform }
    }

    pub fn identity(graph: &Graph<CS>) -> Automorphism<CS> {
        let transform = AffineMap::identity(graph.dim());

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


impl<CS> Mul<&Automorphism<CS>> for &Automorphism<CS>
    where
        CS: Clone + Eq + Hash,
        for <'a> &'a AffineMap<CS>: Mul<&'a AffineMap<CS>, Output=AffineMap<CS>>
{
    type Output = Automorphism<CS>;

    fn mul(self, rhs: &Automorphism<CS>) -> Self::Output {
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


fn traverse_with_shift_adjustments<CS>(g: &Graph<CS>, v0: &Vertex)
    -> Vec<Edge<CS>>
    where CS: Clone + Eq + Hash + Ord
{
    let mut shifts = BTreeMap::from([(*v0, Shift::zero(g.dim()))]);
    let mut seen = HashSet::new();
    let mut queue = VecDeque::from([*v0]);
    let mut result = Vec::new();

    while let Some(v) = queue.pop_front() {
        for e in g.incidences(&v) {
            if !shifts.contains_key(&e.tail) {
                shifts.insert(*(&e.tail), &shifts[&v] + &e.shift);
                queue.push_back(e.tail);
            }

            let e = e.canonical();
            if !seen.contains(&e) {
                let shift = &shifts[&e.head] + &e.shift - &shifts[&e.tail];
                result.push(Edge::new(e.head, e.tail, shift));
                seen.insert(e);
            }
        }
    }

    result
}


pub fn graph_component_measures<CS>(g: &Graph<CS>, v0: &Vertex)
    -> (usize, usize, Option<i32>) // TODO return a struct?
    where CS: Clone + Eq + Hash + Ord
{
    let edges = traverse_with_shift_adjustments(g, v0);

    let mut vertices = BTreeSet::new();
    for e in &edges {
        vertices.insert(e.head);
        vertices.insert(e.tail);
    }

    let mut basis = vec![];
    for e in &edges {
        extend_basis(&e.shift.iter().cloned().collect::<Vec<_>>(), &mut basis);
    }

    let size = vertices.len();
    let rank = basis.len();

    let multiplicity = if rank == g.dim() {
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


fn barycentric_placement<CS>(g: &Graph<CS>) -> BTreeMap<Vertex, Point<CS>>
    where CS: Clone + Eq + Hash + Ord
{
    let verts = g.vertices();
    let vidcs: BTreeMap<_, _> =
        verts.iter().enumerate().map(|(i, &e)| (e, i)).collect();

    let n = verts.len();
    let d = g.dim();

    let mut a = vec![vec![0 as i64; n]; n];
    let mut t = vec![vec![0 as i64; d]; n];
    a[0][0] = 1;

    for i in 1..n {
        for ngb in g.incidences(&verts[i]) {
            let j = vidcs[&ngb.tail];
            a[i][j] -= 1;
            a[i][i] += 1;

            let s = ngb.shift;
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


fn edges_by_unique_deltas<CS>(g: &Graph<CS>)
    -> BTreeMap<Vertex, HashMap<Vector<CS>, Edge<CS>>>
    where CS: Clone + Eq + Hash + Ord
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
