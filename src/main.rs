use std::fmt::Display;
use std::hash::Hash;

use rust_systre::pgraphs::*;


fn main() {
    test_graph_examples();
}


fn test_graph_examples() {
    test_graph(graph2d(&[
        [1, 2, 0, 0],
        [1, 2, 1, 0],
        [1, 2, 0, 1],
    ]));

    test_graph(graph3d(&[
        [1, 2, 0, 0, 0],
        [1, 2, 0, 1, 0],
        [1, 2, 0, 0, 1],
    ]));

    test_graph(graph3d(&[
        [1, 2, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [1, 2, 0, 1, 0],
        [1, 2, 0, 0, 1],
    ]));

    test_graph(graph2d(&[
        [1, 1, 1, 0],
        [1, 1, 0, 1],
        [1, 2, 0, 0],
        [1, 2, 1, 1],
        [1, 3, 0, 0],
        [1, 3, 1, -1],
    ]));

    test_graph(graph2d(&[
        [1, 1, 1, 0],
        [1, 1, 0, 1],
        [1, 2, 0, 0],
        [1, 2, 1, 1],
        [1, 3, 0, 0],
        [1, 3, 1, -1],
        [1, 4, 0, 0],
        [1, 4, 1, -1],
    ]));

    test_graph(graph3d(&[
        [1, 1, -1, 1, 1],
        [1, 1, 0, -1, 1],
        [1, 1, 0, 0, -1],
    ]));

    test_graph(graph3d(&[
        [1, 2, 0, 0, 0],
        [1, 2, 2, 0, 0],
        [1, 2, 0, 2, 0],
        [1, 2, 0, 0, 2],
    ]));

    test_graph(graph3d(&[
        [1, 3, 0, 0, 0],
        [1, 3, 2, 0, 0],
        [1, 3, 0, 2, 0],
        [1, 3, 0, 0, 2],
        [2, 4, 0, 0, 0],
        [2, 4, 2, 0, 0],
        [2, 4, 0, 2, 0],
        [2, 4, 0, 0, 2],
        ]));

    test_graph(graph3d(&[
        [1, 2, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [4, 5, 0, 0, 0],
        [5, 6, 0, 0, 0],
        [6, 1, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [2, 3, 0, 1, 0],
        [3, 4, 0, 0, 1],
        [4, 5, -1, 0, 0],
        [5, 6, 0, -1, 0],
        [6, 1, 0, 0, -1],
    ]));
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
struct World {}


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


fn test_graph<T, CS>(g: Graph<T, CS>)
    where
        T: LabelVector<CS> + Display,
        CS: Clone + Eq + Hash + Ord
{
    println!("Graph:");
    println!("{}", g);

    println!("With normalized shifts:");
    println!("{}", g.shift_normalized());

    println!("Coordination sequences:");
    for v in g.vertices() {
        print!("{}: ", v);
        for n in g.coordination_sequence(&v).take(20) {
            print!(" {}", n);
        }
        println!();
    }
    println!();

    println!("Connected: {}", g.is_connected());
    println!();

    if g.is_connected() {
        for v in g.vertices() {
            println!("pos({}) = ({})", v, point_to_string(&g.position(&v)));
        }
        println!();
        println!("Stable: {}", g.is_stable());
        println!("Locally stable: {}", g.is_locally_stable());
        println!(
            "Second order collisions: {}", g.has_second_order_collisions()
        );
        println!();

        for v in g.vertices() {
            for ngb in g.incidences(&v) {
                let delta = g.edge_vector(&ngb);
                if let Some(e) = g.edge_by_unique_delta(&v, &delta) {
                    println!("{} => ({}) => {}", v, vector_to_string(&delta), e);
                } else {
                    println!("{} => ({}) => None", v, vector_to_string(&delta));
                }
            }
        }
    }
    println!();

    println!("==========");
    println!();
}

fn point_to_string(p: &Point) -> String {
    p.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
}

fn vector_to_string(p: &Vector) -> String {
    p.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
}
