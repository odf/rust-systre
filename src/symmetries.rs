use crate::pgraphs::*;
use crate::arithmetic::linear_algebra::extend_basis;


impl<T: LabelVector> Graph<T>
{
    fn good_edge_chains(&self) -> Vec<Vec<VectorLabelledEdge<T>>>
        where T: LabelVector
    {
        let mut result = vec![];
        for &e in &self.edges {
            for e in [e, -e] {
                generate_edge_chain_extensions(vec![e], &self, &mut result);
            }
        }
        result
    }
}


fn generate_edge_chain_extensions<T: LabelVector>(
    edges: Vec<VectorLabelledEdge<T>>,
    graph: &Graph<T>,
    result: &mut Vec<Vec<VectorLabelledEdge<T>>>
) {
    if edges.len() == T::dim() as usize {
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
