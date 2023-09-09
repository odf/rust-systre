use num_rational::Ratio;
use rust_systre::arithmetic::matrices::Matrix;
use rust_systre::arithmetic::geometry::{
    Vector, Point, ScalarProduct, AffineMap, CoordinateMap
};

#[derive(Clone)]
struct World {}

fn main() {
    let r = |n, d| Ratio::new(n, d);

    println!(
        "Matrix:\n{}",
        Matrix::new(3, &[1, 2, 3, 4, 5, 6])
    );
    println!(
        "Column:\n{}",
        Matrix::col(&[1.2, 1.4, 1.6])
    );
    println!(
        "Vector:\n{}",
        Vector::<_, World>::new(&[0.8, 2.4, 7.2])
    );
    println!(
        "Point:\n{}",
        Point::<_, World>::new(&[r(3, 4), r(4, 3), r(1, 2)])
    );
    println!(
        "ScalarProduct:\n{}",
        ScalarProduct::<Ratio<i32>, World>::default(4)
    );
    println!(
        "AffineMap:\n{}",
        AffineMap::<Ratio<i32>, World>::identity(4)
    );

    println!(
        "CoordinateMap:\n{}",
        CoordinateMap::<i32, World, World>::new(
            &AffineMap::<i32, World>::identity(4)
        )
    );
}
