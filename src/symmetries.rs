use std::collections::{HashMap, VecDeque};

use itertools::Itertools;

use crate::pgraphs::*;
use crate::arithmetic::linear_algebra::extend_basis;


fn automorphism<T: LabelVector + std::fmt::Debug>(
    graph: &Graph<T>,
    seed_src: &Vertex,
    seed_img: &Vertex,
    transform: &AffineMap
)
    -> Option<Automorphism<T>>
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

                    match graph.edge_by_unique_delta(&vsrc, &dimg) {
                        None => return None,
                        Some(eimg) => {
                            edge_map.insert(esrc, eimg);
                            queue.push_back((esrc.tail, eimg.tail));
                        }
                    };
                }
            },
            _ => {},
        };
    }

    Some(Automorphism::new(vertex_map, edge_map, transform.clone()))
}


fn characteristic_edge_lists<T: LabelVector>(graph: &Graph<T>)
    -> Vec<Vec<Edge<T>>>
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


fn good_edge_chains<T>(graph: &Graph<T>) -> Vec<Vec<Edge<T>>>
    where T: LabelVector
{
    let mut result = vec![];
    for e in graph.directed_edges() {
        generate_edge_chain_extensions(vec![e], graph, &mut result);
    }
    result
}


fn generate_edge_chain_extensions<T: LabelVector>(
    edges: Vec<Edge<T>>,
    graph: &Graph<T>,
    result: &mut Vec<Vec<Edge<T>>>
) {
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


fn good_combinations<T: LabelVector>(
    edges: &Vec<Edge<T>>, graph: &Graph<T>
) -> Vec<Vec<Edge<T>>>
{
    let mut result = vec![];
    for es in edges.iter().combinations(T::dim()) {
        let es = es.iter().map(|&&e| e).collect::<Vec<_>>();
        if are_linearly_independent(&es, graph) {
            for es in es.iter().permutations(T::dim()) {
                result.push(es.into_iter().cloned().collect());
            }
        }
    }
    result
}


fn are_linearly_independent<T: LabelVector>(
    edges: &Vec<Edge<T>>, graph: &Graph<T>
) -> bool
{
    let mut basis = vec![];
    for e in edges {
        let v: Vec<_> = graph.edge_vector(&e)
            .into_iter().collect();
        extend_basis(&v, &mut basis);
    }
    basis.len() == edges.len()
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::arithmetic::matrices::Matrix;

    use super::*;

    #[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
    struct World {}

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    fn graph2d(spec: &[[i32; 4]]) -> Graph<LabelVector2d<World>> {
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

    fn sql() -> Graph<LabelVector2d<World>> {
        graph2d(&[
            [1, 1, 1, 0],
            [1, 1, 0, 1],
        ])
    }

    #[test]
    fn test_identity_automorphism() {
        let a = automorphism(&sql(), &1, &1, &AffineMap::identity(2)).unwrap();

        assert_eq!(a.transform, AffineMap::identity(2));
        assert_eq!(a.vertex_map, HashMap::from([(1, 1)]));
        assert_eq!(a.edge_map.len(), 4);
    }

    #[test]
    fn test_nontrivial_automorphism() {
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
    fn test_characteristic_edge_lists() {
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
}
