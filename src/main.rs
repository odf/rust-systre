use rust_systre::pgraphs::*;


fn main() {
    let edge = |u, v, x, y, z| {
        VectorLabelledEdge::new(u, v, LabelVector3d::new(x, y, z))
    };

    let g = Graph::new(&[
        edge(1, 2, 0, 0, 0),
        edge(1, 2, 1, 0, 0),
        edge(1, 2, 0, 1, 0),
        edge(1, 2, 0, 0, 1),
    ]);

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
    println!();

    let a = 51;
    let b = 69;
    println!("gcdx({}, {}) = {:?}", a, b, gcdx(a, b));
}
