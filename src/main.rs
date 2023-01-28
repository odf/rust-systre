use rust_systre::pgraphs::*;


fn main() {
    let e1 = VectorLabelledEdge::new(1, 1, LabelVector2d::zero());
    let e2 = VectorLabelledEdge::new(1, 1, LabelVector2d::new(1, 0));

    let g = Graph::new(&[e1, e2]);
    println!("{}", g);
}
