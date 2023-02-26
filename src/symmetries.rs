use itertools::Itertools;

use crate::pgraphs::{Graph, LabelVector, VectorLabelledEdge, Vertex};
use crate::arithmetic::linear_algebra::extend_basis;


fn characteristic_edge_lists<T: LabelVector>(graph: &Graph<T>)
    -> Vec<Vec<VectorLabelledEdge<T>>>
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


fn good_edge_chains<T>(graph: &Graph<T>) -> Vec<Vec<VectorLabelledEdge<T>>>
    where T: LabelVector
{
    let mut result = vec![];
    for e in directed_edge(graph) {
        generate_edge_chain_extensions(vec![e], graph, &mut result);
    }
    result
}


fn generate_edge_chain_extensions<T: LabelVector>(
    edges: Vec<VectorLabelledEdge<T>>,
    graph: &Graph<T>,
    result: &mut Vec<Vec<VectorLabelledEdge<T>>>
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
    edges: &Vec<VectorLabelledEdge<T>>, graph: &Graph<T>
) -> Vec<Vec<VectorLabelledEdge<T>>>
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


fn directed_edge<T: LabelVector>(graph: &Graph<T>) -> Vec<VectorLabelledEdge<T>>
{
    graph.vertices().iter().flat_map(|v| graph.incidences(v)).collect()
}


fn are_linearly_independent<T: LabelVector>(
    edges: &Vec<VectorLabelledEdge<T>>, graph: &Graph<T>
) -> bool
{
    let mut basis = vec![];
    for VectorLabelledEdge { head, tail, shift } in edges {
        extend_basis(&graph.edge_vector(&head, &tail, &shift), &mut basis);
    }
    basis.len() == edges.len()
}
