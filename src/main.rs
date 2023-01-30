use rust_systre::pgraphs::*;


fn main() {
    let edge = |u, v, x, y| {
        VectorLabelledEdge::new(u, v, LabelVector2d::new(x, y))
    };

    let g = Graph::new(&[
        edge(2, 1, 0, 0),
        edge(1, 2, 1, 0),
        edge(1, 2, 0, 1),
        edge(2, 1, -1, 0),
    ]);

    println!("Edges:");
    println!("{}", g);

    println!("Vertices: {:?}", g.vertices());
    println!();

    println!("Incidences:");
    let adj = g.incidences();
    for v in g.vertices() {
        print!("{}: ", v);
        for a in &adj[&v] {
            print!(" {},", a);
        }
        println!();
    }
}
