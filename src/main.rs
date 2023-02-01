use std::fmt::Display;

use rust_systre::pgraphs::*;


fn main() {
    test_graph(graph2d(&[
        [1, 1, 1, 2],
        [1, 1, 2, 1],
    ]));

    test_graph(graph3d(&[
        [1, 2, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [1, 2, 0, 1, 0],
        [1, 2, 0, 0, 1],
    ]));

    let a = 39;
    let b = 51;
    println!("gcdx({}, {}) = {:?}", a, b, gcdx(a, b));


    let a = 39;
    let b = 78;
    println!("gcdx({}, {}) = {:?}", a, b, gcdx(a, b));
    println!();

    let mut basis = vec![];
    for v in [[8, 13], [3, 5]] {
        extend_basis(&v, &mut basis);
        println!("add vector {:?} => {:?}", &v, &basis);
    }
}


fn graph2d(spec: &[[i16; 4]]) -> Graph<LabelVector2d> {
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


fn graph3d(spec: &[[i16; 5]]) -> Graph<LabelVector3d> {
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
    where T: LabelVector<Item = i16> + Display
{
    println!("Edges:");
    println!("{}", g);

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
