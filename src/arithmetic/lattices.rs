use std::ops::{Add, Sub, Mul, Neg};

use num_traits::{Zero, One, Signed};


pub trait Coord:
    Clone + PartialOrd + Zero + One + Signed
    + for <'a> Add<&'a Self, Output=Self>
{
    fn div_rounded(&self, other: &Self) -> Self;
}


pub trait CoordPtr<T>:
    Sized
    + Add<Output=T> + Sub<Output=T> + Sub<T, Output=T> + Neg<Output=T>
    + Mul<Output=T>
{
}


pub fn dirichlet_vectors<T>(basis: &[Vec<T>], epsilon: T) -> Vec<Vec<T>>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    let n = basis.len();
    assert!(basis.iter().all(|v| v.len() == n));

    match n {
        0 | 1 => basis.to_vec(),
        2 => {
            let (u, v) = gauss_reduced(&basis[0], &basis[1], epsilon);
            let s = (0..2).map(|i| &u[i] + &v[i]).collect();
            vec![u, v, s]
        }
        3 => {
            let (u, v, w) = selling_reduced(
                &basis[0], &basis[1], &basis[2], epsilon
            );
            let uv = (0..2).map(|i| &u[i] + &v[i]).collect();
            let uw = (0..2).map(|i| &u[i] + &w[i]).collect();
            let vw = (0..2).map(|i| &v[i] + &w[i]).collect();
            let uvw = (0..2).map(|i| &u[i] + &v[i] + &w[i]).collect();
            vec![u, v, w, uv, uw, vw, uvw]
        }
        _ => panic!("dimension {} not supported", n)
    }
}


fn gauss_reduced<T>(u: &[T], v: &[T], epsilon: T) -> (Vec<T>, Vec<T>)
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    let norm = |v: &[T]| dot(v, v);
    let fudge = &T::one() - &epsilon;
    let mut u = u.to_vec();
    let mut v = v.to_vec();

    if norm(&v) < norm(&u) {
        (u, v) = (v, u);
    }

    while norm(&u) < &norm(&v) * &fudge {
        let t = dot(&u, &v).div_rounded(&norm(&u));
        let w: Vec<_> = (0..u.len()).map(|i| &v[i] - &t * &u[i]).collect();
        (u, v) = (w, u);
    }

    if dot(&u, &v).is_positive() {
        v = v.into_iter().map(|x| -x).collect();
    }

    (u, v)
}


fn selling_reduced<T>(u: &[T], v: &[T], w: &[T], epsilon: T)
    -> (Vec<T>, Vec<T>, Vec<T>)
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    let s: Vec<_> = (0..u.len()).map(|i| -(&u[i] + &v[i] + &w[i])).collect();
    let mut vs = [u.to_vec(), v.to_vec(), w.to_vec(), s];

    loop {
        let mut changed = false;

        for i in 0..3 {
            for j in (i + 1)..4 {
                if dot(&vs[i], &vs[j]) > epsilon {
                    for k in 0..4 {
                        if k != i && k != j {
                            vs[k] = (0..u.len()).map(|mu|
                                &vs[k][mu] + &vs[i][mu]
                            ).collect();
                        }
                    }
                    vs[i] = vs[i].iter().map(|x| -x).collect();
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    (vs[0].clone(), vs[1].clone(), vs[2].clone())
}


fn dot<T>(v: &[T], w: &[T]) -> T
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    assert_eq!(v.len(), w.len());

    let mut s = T::zero();
    for i in 0..v.len() {
        s = s + &v[i] * &w[i];
    }

    s
}


fn positive_direction<T>(v: &[T]) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    if let Some(x) = v.iter().find(|x| !x.is_zero()) {
        if x.is_negative() {
            return v.iter().map(|x| -x).collect();
        }
    }
    v.to_vec()
}
