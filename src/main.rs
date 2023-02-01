use std::fmt::Display;

use rust_systre::pgraphs::*;


fn main() {
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


fn graph2d(spec: &[[i32; 4]]) -> Graph<LabelVector2d> {
    let mut edges = vec![];

    for [u, v, x, y] in spec {
        edges.push(VectorLabelledEdge::new(
            *u as u32,
            *v as u32,
            LabelVector2d::new(*x, *y)
        ));
    }

    Graph::new(&edges)
}


fn graph3d(spec: &[[i32; 5]]) -> Graph<LabelVector3d> {
    let mut edges = vec![];

    for [u, v, x, y, z] in spec {
        edges.push(VectorLabelledEdge::new(
            *u as u32,
            *v as u32,
            LabelVector3d::new(*x, *y, *z)
        ));
    }

    Graph::new(&edges)
}


fn test_graph<T>(g: Graph<T>)
    where T: LabelVector + Display
{
    println!("Edges:");
    println!("{}", g);

    println!("Is connected: {}", g.is_connected());
    println!();

    println!("Vertices: {:?}", g.vertices());
    println!();

    println!("Incidences:");
    for v in g.vertices() {
        print!("{}: ", v);
        for a in g.incidences(&v) {
            print!(" {},", a);
        }
        println!();
    }
    println!();

    println!("Coordination sequences:");
    for v in g.vertices() {
        print!("{}: ", v);
        for n in g.coordination_sequence(&v).take(20) {
            print!(" {}", n);
        }
        println!();
    }
    println!();

    let g = g.shift_normalized();

    println!("With normalized shifts:");
    println!("{}", g);

    let (size, rank, m) = graph_component_measures(&g, &1);
    println!("Component: size = {}, rank = {}, multiplicity = {:?}", size, rank, m);
    println!();

    println!("==========");
    println!();
}
