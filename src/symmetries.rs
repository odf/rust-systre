use std::collections::{HashMap, VecDeque};
use std::hash::Hash;

use itertools::Itertools;

use crate::pgraphs::*;
use crate::arithmetic::linear_algebra::extend_basis;


fn automorphism<CS>(
    graph: &Graph<CS>,
    seed_src: &Vertex,
    seed_img: &Vertex,
    transform: AffineMap<CS>
)
    -> Option<Automorphism<CS>>
    where CS: Clone + Eq + Hash + Ord
{
    let mut vertex_map = HashMap::from([(*seed_src, *seed_img)]);
    let mut edge_map = HashMap::new();
    let mut queue = VecDeque::from([(*seed_src, *seed_img)]);

    while let Some((vsrc, vimg)) = queue.pop_front() {
        match vertex_map.get(&vsrc) {
            Some(&w) if w != vimg => return None,
            None => {
                vertex_map.insert(vsrc, vimg);
                for esrc in graph.incidences(&vsrc) {
                    let dsrc = graph.edge_vector(&esrc);
                    let dimg = &transform * &dsrc;

                    match graph.edge_by_unique_delta(&vsrc, &dimg) {
                        None => return None,
                        Some(eimg) => {
                            edge_map.insert(esrc.clone(), eimg.clone());
                            queue.push_back((esrc.tail, eimg.tail));
                        }
                    };
                }
            },
            _ => {},
        };
    }

    Some(Automorphism::new(vertex_map, edge_map, transform))
}


fn characteristic_edge_lists<CS>(graph: &Graph<CS>)
    -> Vec<Vec<Edge<CS>>>
    where CS: Clone + Eq + Hash + Ord
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


fn good_edge_chains<CS>(graph: &Graph<CS>)
    -> Vec<Vec<Edge<CS>>>
    where CS: Clone + Eq + Hash + Ord
{
    let mut result = vec![];
    for e in graph.directed_edges() {
        generate_edge_chain_extensions(vec![e], graph, &mut result);
    }
    result
}


fn generate_edge_chain_extensions<CS>(
    edges: Vec<Edge<CS>>,
    graph: &Graph<CS>,
    result: &mut Vec<Vec<Edge<CS>>>
)
    where CS: Clone + Eq + Hash + Ord
{
    if edges.len() == graph.dim() {
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


fn good_combinations<CS>(edges: &Vec<Edge<CS>>, graph: &Graph<CS>)
    -> Vec<Vec<Edge<CS>>>
    where CS: Clone + Eq + Hash + Ord
{
    let mut result = vec![];
    for es in edges.iter().combinations(graph.dim()) {
        let es = es.into_iter().cloned().collect::<Vec<_>>();
        if are_linearly_independent(&es, graph) {
            for es in es.iter().permutations(graph.dim()) {
                result.push(es.into_iter().cloned().collect());
            }
        }
    }
    result
}


fn are_linearly_independent<CS>(edges: &Vec<Edge<CS>>, graph: &Graph<CS>)
    -> bool
    where CS: Clone + Eq + Hash + Ord
{
    let mut basis = vec![];
    for e in edges {
        let v: Vec<_> = graph.edge_vector(&e)
            .into_iter().collect();
        extend_basis(&v, &mut basis);
    }
    basis.len() == edges.len()
}
