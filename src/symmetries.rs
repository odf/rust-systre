use std::collections::{HashMap, VecDeque};

use itertools::Itertools;

use crate::pgraphs::*;
use crate::arithmetic::linear_algebra::extend_basis;


fn automorphism<T: LabelVector>(
    graph: &Graph<T>,
    seed_src: &Vertex,
    seed_img: &Vertex,
    transform: AffineMap
)
    -> Option<Automorphism<T>>
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
                            edge_map.insert(esrc, eimg);
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
        result = good_combinations(&directed_edge(graph), graph);
    }
    result
}


fn good_edge_chains<T>(graph: &Graph<T>) -> Vec<Vec<Edge<T>>>
    where T: LabelVector
{
    let mut result = vec![];
    for e in directed_edge(graph) {
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


fn directed_edge<T: LabelVector>(graph: &Graph<T>) -> Vec<Edge<T>>
{
    graph.vertices().iter().flat_map(|v| graph.incidences(v)).collect()
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
