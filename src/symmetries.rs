use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, VecDeque, HashSet};
use std::hash::{Hash, Hasher};

use itertools::Itertools;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero, ToPrimitive, Signed};
use num_integer::Integer;

use crate::arithmetic::matrices::Matrix;
use crate::partitions::Partition;
use crate::pgraphs::*;
use crate::arithmetic::linear_algebra::{extend_basis, LinearAlgebra};


pub fn ladder_symmetries<T, CS>(graph: &Graph<T, CS>) -> Vec<Automorphism<T, CS>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    assert!(graph.is_locally_stable(), "graph must be locally stable");

    let orbit = &raw_translational_orbits(&graph)[0];
    let p0 = mod1(&graph.position(&orbit[0]));
    let id = AffineMap::identity(T::dim());

    orbit.iter().skip(1)
        .filter(|v| mod1(&graph.position(v)) == p0)
        .map(|v| automorphism(&graph, &orbit[0], &v, &id).unwrap())
        .collect::<Vec<_>>()
}


pub fn is_ladder<T, CS>(graph: &Graph<T, CS>) -> bool
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    assert!(graph.is_locally_stable(), "graph must be locally stable");

    if graph.is_stable() {
        false
    } else {
        let orbit = &raw_translational_orbits(&graph)[0];
        let p0 = mod1(&graph.position(&orbit[0]));

        orbit.iter().skip(1).any(|v| mod1(&graph.position(v)) == p0)
    }
}


pub fn is_minimal<T, CS>(graph: &Graph<T, CS>) -> bool
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    assert!(graph.is_locally_stable(), "graph must be locally stable");
    translational_orbits(&graph)[0].len() == 1
}


pub fn minimal_image<TIn, TOut, CSIn, CSOut>(graph: &Graph<TIn, CSIn>)
    -> Option<Graph<TOut, CSOut>>
    where
        TIn: LabelVector<CSIn>,
        TOut: LabelVector<CSOut>,
        CSIn: Clone + Eq + Hash + Ord,
        CSOut: Clone + Eq + Hash + Ord,
{
    assert!(graph.is_locally_stable(), "graph must be locally stable");

    let orbits = translational_orbits(graph);

    if orbits[0].len() == 1 {
        None
    } else {
        let b = extended_translation_basis(graph, &orbits[0]);
        let rows = b.iter().map(|v| Matrix::row(v)).collect::<Vec<_>>();
        let m = Matrix::vstack(&rows).inverse().unwrap();
        let aff = &AffineMap::new(&m, &Vector::zero(TIn::dim()));

        let mut imgs = HashMap::new();
        let mut shifts = HashMap::new();

        for i in 0..orbits.len() {
            let p0 = graph.position(&orbits[i][0]);
            for v in &orbits[i] {
                imgs.insert(v, i + 1);
                shifts.insert(v, &graph.position(&v) - &p0);
            }
        }

        let mut edges: Vec<Edge<TOut, CSOut>> = vec![];
        for e in graph.directed_edges() {
            if e == e.canonical() {
                let head = *imgs.get(&e.head).unwrap() as u32;
                let tail = *imgs.get(&e.tail).unwrap() as u32;
                let s: Vector = e.shift.to_vec().iter()
                    .map(|x| BigRational::from(BigInt::from(*x)))
                    .collect();
                let shift = aff * (
                    s +
                    shifts.get(&e.tail).unwrap() -
                    shifts.get(&e.head).unwrap()
                );
                let shift = shift.iter()
                    .map(|x| r_to_i32(x).unwrap())
                    .collect::<Vec<_>>();

                edges.push(Edge::new(head, tail, TOut::from_vec(&shift)));
            }
        }
        Some(Graph::new(&edges))
    }
}


pub fn translational_orbits<T, CS>(graph: &Graph<T, CS>) -> Vec<Vec<Vertex>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    assert!(graph.is_locally_stable(), "graph must be locally stable");
    translational_equivalences(&graph).classes(&graph.vertices())
}


pub fn symmetries<T, CS>(graph: &Graph<T, CS>) -> Vec<Automorphism<T, CS>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    // TODO tag edge list equivalence classes that make no automorphism
    assert!(graph.is_locally_stable(), "graph must be locally stable");

    let edge_lists = characteristic_edge_lists(&graph);
    let keys: Vec<_> = edge_lists.iter().map(hash_value).collect();
    let seeds: Vec<_> = map_vec(|es| es[0].head, &edge_lists);
    let bases: Vec<_> = map_vec(
        |es| matrix_from_edge_list(&es, &graph),
        &edge_lists
    );

    let inv_b = bases[0].inverse().unwrap();

    let mut gens = vec![];
    let mut p = Partition::<u64>::new();

    for i in 0..edge_lists.len() {
        if p.find(&keys[i]) != p.find(&keys[0]) {
            let m = &bases[i] * &inv_b;

            if is_unimodular(&m) {
                let aff = &AffineMap::new(&m, &Vector::zero(T::dim()));
                let iso = automorphism(&graph, &seeds[0], &seeds[i], &aff);

                if let Some(iso) = iso {
                    for k in 0..edge_lists.len() {
                        let mapped = map_vec(
                            |e| iso.edge_map.get(&e).unwrap(),
                            &edge_lists[k]
                        );
                        p.unite(&keys[k], &hash_value(&mapped));
                    }
                    gens.push(iso);
                }
            }
        }
    }

    let mut result = HashSet::from([Automorphism::identity(&graph)]);
    let mut queue = VecDeque::from([Automorphism::identity(&graph)]);

    while let Some(phi) = queue.pop_front() {
        for psi in &gens {
            let product = &phi * psi;
            if !result.contains(&product) {
                queue.push_back(product.clone());
                result.insert(product);
            }
        }
    }

    result.into_iter().collect()
}


fn matrix_from_edge_list<T, CS>(
    es: &Vec<Edge<T, CS>>,
    graph: &Graph<T, CS>
)
    -> Matrix<num_rational::Ratio<BigInt>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    Matrix::vstack(&map_vec(
        |e| Matrix::row(&graph.edge_vector(e).into_iter().collect::<Vec<_>>()),
        &es
    ))
}


fn hash_value<T: Hash>(x: T) -> u64 {
    let mut s = DefaultHasher::new();
    x.hash(&mut s);
    s.finish()
}


fn is_unimodular(m: &Matrix<BigRational>) -> bool {
    m.iter().all(|x| x % BigRational::one() == BigRational::zero())
    &&
    m.determinant().abs() == BigRational::one()
}


fn translational_equivalences<T, CS>(graph: &Graph<T, CS>) -> Partition<u32>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let equivs = raw_translational_equivalences(&graph);
    let orbit = &equivs.classes(&graph.vertices())[0];
    let p0 = mod1(&graph.position(&orbit[0]));
    let id = AffineMap::identity(T::dim());

    if orbit.iter().skip(1).all(|v| mod1(&graph.position(v)) != p0) {
        equivs
    } else {
        let mut p = Partition::new();

        for b in extended_translation_basis(graph, &orbit) {
            let b: Vector = b.iter().cloned().collect();
            let v = orbit.iter()
                .find(|v| mod1(&(graph.position(v) + &b)) == p0)
                .unwrap();
            let iso = automorphism(&graph, &orbit[0], &v, &id).unwrap();
            for w in &graph.vertices() {
                p.unite(w, &iso.vertex_map[w]);
            }
        }
        p
    }
}


fn raw_translational_orbits<T, CS>(graph: &Graph<T, CS>) -> Vec<Vec<Vertex>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    raw_translational_equivalences(&graph).classes(&graph.vertices())
}


fn raw_translational_equivalences<T, CS>(graph: &Graph<T, CS>) -> Partition<u32>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let id = &AffineMap::identity(T::dim());
    let verts = graph.vertices();
    let mut p = Partition::new();

    for v in &verts {
        if p.find(&verts[0]) != p.find(&v) {
            if let Some(iso) = automorphism(&graph, &verts[0], &v, id) {
                for w in &verts {
                    p.unite(&w, &iso.vertex_map[w]);
                }
            }
        }
    }

    p
}


fn mod1(p: &Point) -> Point {
    p.iter().map(|x| x % BigInt::one()).collect()
}


fn extended_translation_basis<T, CS>(graph: &Graph<T, CS>, orbit: &[Vertex]) 
    -> Vec<Vec<BigRational>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let one = BigInt::one();
    let p0 = &graph.position(&orbit[0]);

    let deltas: Vec<_> = orbit.iter().skip(1)
        .map(|v| (graph.position(&v) - p0).iter()
            .map(|x| x % &one)
            .collect::<Vec<_>>()
        )
        .collect();

    let mut common_denom = one;
    for d in &deltas {
        for x in d {
            common_denom *= x.denom() / common_denom.gcd(x.denom());
        }
    }
    let common_denom = common_denom.to_i32().unwrap();
    let f = BigRational::from(BigInt::from(common_denom));

    let id: Matrix<i32> = Matrix::identity(T::dim()) * common_denom;
    let mut basis: Vec<_> = (0..T::dim()).map(|i| id.get_row(i)).collect();

    for d in deltas {
        let d: Vec<_> = d.iter().map(|x| r_to_i32(&(x * &f)).unwrap()).collect();
        extend_basis(&d, &mut basis);
    }

    basis.iter()
        .map(|b| b.iter()
            .map(|x| BigRational::from(BigInt::from(*x)) / &f)
            .collect::<Vec<_>>()
        )
        .collect()
}


fn r_to_i32(x: &BigRational) -> Option<i32> {
    if x.is_integer() {
        x.to_integer().to_i32()
    } else {
        None
    }
}


fn automorphism<T, CS>(
    graph: &Graph<T, CS>,
    seed_src: &Vertex,
    seed_img: &Vertex,
    transform: &AffineMap
)
    -> Option<Automorphism<T, CS>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut vertex_map = HashMap::new();
    let mut edge_map = HashMap::new();
    let mut queue = VecDeque::from([(*seed_src, *seed_img)]);

    while let Some((vsrc, vimg)) = queue.pop_front() {
        match vertex_map.get(&vsrc) {
            Some(&w) if w != vimg => return None,
            None => {
                vertex_map.insert(vsrc, vimg);
                for esrc in graph.incidences(&vsrc) {
                    let dsrc = graph.edge_vector(&esrc);
                    let dimg = transform * &dsrc;

                    match graph.edge_by_unique_delta(&vimg, &dimg) {
                        None => return None,
                        Some(eimg) => {
                            queue.push_back((esrc.tail, eimg.tail));
                            edge_map.insert(esrc, eimg);
                        }
                    };
                }
            },
            _ => {},
        };
    }

    Some(Automorphism::new(vertex_map, edge_map, transform.clone()))
}


fn characteristic_edge_lists<T, CS>(graph: &Graph<T, CS>)
    -> Vec<Vec<Edge<T, CS>>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut result = vec![];

    for v in graph.vertices() {
        for es in good_combinations(&graph.incidences(&v), graph) {
            result.push(es);
        }
    }

    if result.is_empty() {
        result = good_edge_chains(graph);
    }

    if result.is_empty() {
        result = good_combinations(&graph.directed_edges(), graph);
    }
    result
}


fn good_edge_chains<T, CS>(graph: &Graph<T, CS>)
    -> Vec<Vec<Edge<T, CS>>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut result = vec![];
    for e in graph.directed_edges() {
        generate_edge_chain_extensions(vec![e], graph, &mut result);
    }
    result
}


fn generate_edge_chain_extensions<T, CS>(
    edges: Vec<Edge<T, CS>>,
    graph: &Graph<T, CS>,
    result: &mut Vec<Vec<Edge<T, CS>>>
)
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    if edges.len() == T::dim() {
        result.push(edges);
    } else {
        for e in graph.incidences(&edges[edges.len() - 1].tail) {
            let mut next = edges.clone();
            next.push(e);

            if are_linearly_independent(&next, graph) {
                generate_edge_chain_extensions(next, graph, result);
            }
        }
    }
}


fn good_combinations<T, CS>(edges: &Vec<Edge<T, CS>>, graph: &Graph<T, CS>)
    -> Vec<Vec<Edge<T, CS>>>
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut result = vec![];
    for es in edges.iter().combinations(T::dim()) {
        let es = es.iter().map(|&e| e.clone()).collect::<Vec<_>>();
        if are_linearly_independent(&es, graph) {
            for es in es.iter().permutations(T::dim()) {
                result.push(es.iter().map(|&e| e.clone()).collect());
            }
        }
    }
    result
}


fn are_linearly_independent<T, CS>(
    edges: &Vec<Edge<T, CS>>, graph: &Graph<T, CS>
)
    -> bool
    where
        T: LabelVector<CS>,
        CS: Clone + Eq + Hash + Ord
{
    let mut basis = vec![];
    for e in edges {
        let v: Vec<_> = graph.edge_vector(&e)
            .into_iter().collect();
        extend_basis(&v, &mut basis);
    }
    basis.len() == edges.len()
}


fn map_vec<'a, A, B, F>(fun: F, input: &'a Vec<A>) -> Vec<B>
    where F: FnMut(&'a A) -> B
{
    input.iter().map(fun).collect()
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::arithmetic::matrices::Matrix;

    use super::*;

    #[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
    struct World {}

    type Label2d = LabelVector2d<World>;

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    fn graph2d(spec: &[[i32; 4]]) -> Graph<LabelVector2d<World>, World> {
        let mut edges = vec![];

        for [u, v, x, y] in spec {
            edges.push(Edge::new(
                *u as u32,
                *v as u32,
                LabelVector2d::new(*x, *y)
            ));
        }

        Graph::new(&edges)
    }

    fn graph3d(spec: &[[i32; 5]]) -> Graph<LabelVector3d<World>, World> {
        let mut edges = vec![];

        for [u, v, x, y, z] in spec {
            edges.push(Edge::new(
                *u as u32,
                *v as u32,
                LabelVector3d::new(*x, *y, *z)
            ));
        }

        Graph::new(&edges)
    }

    fn sql() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 1, 1, 0],
            [1, 1, 0, 1],
        ])
    }

    fn pcu() -> Graph<LabelVector3d<World>, World> {
        graph3d(&[
            [1, 1, 1, 0, 0],
            [1, 1, 0, 1, 0],
            [1, 1, 0, 0, 1],
        ])
    }

    fn hcb() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 2, 0, 0],
            [1, 2, 1, 0],
            [1, 2, 0, 1],
        ])
    }

    fn dia() -> Graph<LabelVector3d<World>, World> {
        graph3d(&[
            [1, 2, 0, 0, 0],
            [1, 2, 1, 0, 0],
            [1, 2, 0, 1, 0],
            [1, 2, 0, 0, 1],
        ])
    }

    fn sql2() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 2, 0, 0],
            [2, 1, 1, 0],
            [1, 1, 0, 1],
            [2, 2, 0, 1],
        ])
    }

    fn sql4() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 2, 0, 0],
            [2, 1, 1, 0],
            [3, 4, 0, 0],
            [4, 3, 1, 0],
            [1, 3, 0, 0],
            [3, 1, 0, 1],
            [2, 4, 0, 0],
            [4, 2, 0, 1],
        ])
    }

    fn sql_c2() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 1, 1, 0],
            [1, 1, 0, 1],
            [2, 2, 1, 0],
            [2, 2, 0, 1],
            [1, 2, 0, 0],
        ])
    }

    fn dia_c2() -> Graph<LabelVector3d<World>, World> {
        graph3d(&[
            [1, 2, 0, 0, 0],
            [1, 2, 1, 0, 0],
            [1, 2, 0, 1, 0],
            [1, 2, 0, 0, 1],
            [3, 4, 0, 0, 0],
            [3, 4, 1, 0, 0],
            [3, 4, 0, 1, 0],
            [3, 4, 0, 0, 1],
            [1, 3, 0, 0, 0],
            [2, 4, 0, 0, 0],
        ])
    }

    fn sql2_c2() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [1, 2, 0, 0],
            [2, 1, 1, 0],
            [1, 1, 0, 1],
            [2, 2, 0, 1],
            [3, 4, 0, 0],
            [4, 3, 1, 0],
            [3, 3, 0, 1],
            [4, 4, 0, 1],
            [1, 3, 0, 0],
            [2, 4, 0, 0],
        ])
    }

    fn sql_c4() -> Graph<LabelVector2d<World>, World> {
        graph2d(&[
            [ 1, 1, 1, 0 ],
            [ 1, 1, 0, 1 ],
            [ 2, 2, 1, 0 ],
            [ 2, 2, 0, 1 ],
            [ 3, 3, 1, 0 ],
            [ 3, 3, 0, 1 ],
            [ 4, 4, 1, 0 ],
            [ 4, 4, 0, 1 ],
            [ 1, 2, 0, 0 ],
            [ 2, 3, 0, 0 ],
            [ 3, 4, 0, 0 ],
            [ 4, 1, 0, 0 ],
        ])
    }

    #[test]
    fn test_syms_translational_automorphism() {
        let g = sql();
        let a = automorphism(&g, &1, &1, &AffineMap::identity(2)).unwrap();
        assert_eq!(a.transform, AffineMap::identity(2));
        assert_eq!(a.vertex_map, HashMap::from([(1, 1)]));
        assert_eq!(a.edge_map.len(), 4);

        let g = sql2();
        let a = automorphism(&g, &1, &2, &AffineMap::identity(2)).unwrap();
        assert_eq!(a.transform, AffineMap::identity(2));
        assert_eq!(a.vertex_map, HashMap::from([(1, 2), (2, 1)]));
        assert_eq!(a.edge_map.len(), 8);

        let g = sql4();

        let a = automorphism(&g, &1, &2, &AffineMap::identity(2)).unwrap();
        assert_eq!(a.transform, AffineMap::identity(2));
        assert_eq!(
            a.vertex_map,
            HashMap::from([(1, 2), (2, 1), (3, 4), (4, 3)])
        );
        assert_eq!(a.edge_map.len(), 16);

        let a = automorphism(&g, &1, &4, &AffineMap::identity(2)).unwrap();
        assert_eq!(a.transform, AffineMap::identity(2));
        assert_eq!(
            a.vertex_map,
            HashMap::from([(1, 4), (4, 1), (2, 3), (3, 2)])
        );
        assert_eq!(a.edge_map.len(), 16);
    }

    #[test]
    fn test_syms_nontranslational_automorphism() {
        let aff = AffineMap::new(
            &Matrix::new(2, &[r(0), r(-1), r(1), r(0)]),
            &Vector::new(&[r(1), r(0)])
        );
        let a = automorphism(&sql(), &1, &1, &aff).unwrap();

        assert_eq!(a.transform, aff);
        assert_eq!(a.vertex_map, HashMap::from([(1, 1)]));
        assert_eq!(a.edge_map.len(), 4);
    }

    #[test]
    fn test_syms_characteristic_edge_lists() {
        let left = Edge::new(1, 1, LabelVector2d::<World>::new(-1, 0));
        let right = Edge::new(1, 1, LabelVector2d::<World>::new(1, 0));
        let up = Edge::new(1, 1, LabelVector2d::<World>::new(0, -1));
        let down = Edge::new(1, 1, LabelVector2d::<World>::new(0, 1));

        let mut expected: Vec<Vec<_>> = vec![
            vec![left, up], vec![left, down],
            vec![up, left], vec![up, right],
            vec![down, left], vec![down, right],
            vec![right, up], vec![right, down],
        ];
        expected.sort();

        let g = sql();

        let mut els1 = characteristic_edge_lists(&g);
        els1.sort();

        let mut els2 = good_edge_chains(&g);
        els2.sort();

        let mut els3 = good_combinations(&g.directed_edges(), &g);
        els3.sort();

        assert_eq!(els1, expected);
        assert_eq!(els2, expected);
        assert_eq!(els3, expected);
    }

    #[test]
    fn test_syms_raw_equivalence_classes() {
        let g = sql();
        let p = raw_translational_equivalences(&g);
        let cl = p.classes(&g.vertices());
        assert_eq!(cl, [[1]]);

        let g = hcb();
        let p = raw_translational_equivalences(&g);
        let cl = p.classes(&g.vertices());
        assert_eq!(cl, [[1], [2]]);

        let g = sql2();
        let p = raw_translational_equivalences(&g);
        let cl = p.classes(&g.vertices());
        assert_eq!(cl, [[1, 2]]);

        let g = sql4();
        let p = raw_translational_equivalences(&g);
        let cl = p.classes(&g.vertices());
        assert_eq!(cl, [[1, 2, 3, 4]]);

        let g = sql2_c2();
        let p = raw_translational_equivalences(&g);
        let cl = p.classes(&g.vertices());
        assert_eq!(cl, [[1, 2, 3, 4]]);
    }

    #[test]
    fn test_syms_equivalence_classes() {
        let g = sql();
        let vs = g.vertices();
        assert_eq!(
            raw_translational_equivalences(&g).classes(&vs),
            translational_equivalences(&g).classes(&vs),
        );

        let g = sql2();
        let vs = g.vertices();
        assert_eq!(
            raw_translational_equivalences(&g).classes(&vs),
            translational_equivalences(&g).classes(&vs),
        );

        let g = sql4();
        let vs = g.vertices();
        assert_eq!(
            raw_translational_equivalences(&g).classes(&vs),
            translational_equivalences(&g).classes(&vs),
        );
    }

    #[test]
    fn test_syms_equivalence_classes_ladder() {
        let g = sql2_c2();
        let vs = g.vertices();
        assert_eq!(
            translational_equivalences(&g).classes(&vs),
            [[1, 2], [3, 4]],
        );
    }

    #[test]
    fn test_syms_ladder_symmetries() {
        assert_eq!(ladder_symmetries(&sql()).len(), 0);
        assert_eq!(ladder_symmetries(&sql2()).len(), 0);
        assert_eq!(ladder_symmetries(&sql4()).len(), 0);
        assert_eq!(ladder_symmetries(&sql2_c2()).len(), 1);
    }

    #[test]
    #[should_panic]
    fn test_syms_ladder_symmetries_not_locally_stable() {
        assert_eq!(ladder_symmetries(&sql_c4()).len(), 4);
    }

    #[test]
    fn test_syms_is_ladder() {
        assert_eq!(is_ladder(&sql()), false);
        assert_eq!(is_ladder(&hcb()), false);
        assert_eq!(is_ladder(&sql2()), false);
        assert_eq!(is_ladder(&sql4()), false);
        assert_eq!(is_ladder(&sql_c2()), true);
        assert_eq!(is_ladder(&sql2_c2()), true);
    }

    #[test]
    fn test_syms_is_minimal() {
        assert_eq!(is_minimal(&sql()), true);
        assert_eq!(is_minimal(&hcb()), true);
        assert_eq!(is_minimal(&sql2()), false);
        assert_eq!(is_minimal(&sql4()), false);
        assert_eq!(is_minimal(&sql_c2()), true);
        assert_eq!(is_minimal(&sql2_c2()), false);
    }

    #[test]
    fn test_syms_minimal_image() {
        assert!(minimal_image::<_, Label2d, World, World>(&sql()).is_none());
        assert!(minimal_image::<_, Label2d, World, World>(&hcb()).is_none());
        assert!(minimal_image::<_, Label2d, World, World>(&sql_c2()).is_none());

        assert_eq!(
            minimal_image::<_, Label2d, World, World>(&sql2()).unwrap().to_string(),
            sql().to_string()
        );

        assert_eq!(
            minimal_image::<_, Label2d, World, World>(&sql2_c2()).unwrap().to_string(),
            sql_c2().to_string()
        );

        let g = minimal_image::<_, Label2d, World, World>(&sql4()).unwrap();
        assert_eq!(g.vertices().len(), 1);
        assert_eq!(g.directed_edges().len(), 4);
    }

    #[test]
    fn test_syms_symmetries() {
        assert_eq!(symmetries(&sql()).len(), 8);
        assert_eq!(symmetries(&pcu()).len(), 48);
        assert_eq!(symmetries(&hcb()).len(), 12);
        assert_eq!(symmetries(&dia()).len(), 48);
        assert_eq!(symmetries(&sql_c2()).len(), 16);
        assert_eq!(symmetries(&dia_c2()).len(), 96);
    }
}
