use std::cmp::Ordering;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Signed, abs};
use num_integer::Integer;

use crate::arithmetic::linear_algebra::extend_basis;
use super::types::{Coord, CoordPtr};


pub fn shift_for_dirichlet_domain<T, F>(
    p: &[T], vs: &[Vec<T>], dot: F, epsilon: &T
)
    -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>, F: Fn(&[T], &[T]) -> T
{
    let pos = p.to_vec();
    let mut shift = vec![T::zero(); p.len()];
    let mut proceed = true;

    while proceed {
        proceed = false;

        for v in vs {
            let t = dot(&add(&pos, &shift), v) / dot(v, v);
            let t2 = &t + &t;
            if t2 < -T::one() || t2 > T::one() + epsilon {
                shift = sub(&shift, &mul(v, &t.round()));
                proceed = true;
            }
        }
    }

    shift
}


pub fn reduced_lattice_basis<T, F>(vs: &[Vec<T>], dot: F, epsilon: &T)
    -> Option<Vec<Vec<T>>>
    where T: Coord, for <'a> &'a T: CoordPtr<T>, F: Fn(&[T], &[T]) -> T
{
    let dim = vs[0].len();
    let mut a: Vec<Vec<T>> = vec![];

    let mut ws = dirichlet_vectors(vs, &dot, epsilon);
    ws.sort_by(|v, w| {
        let d = dot(v, v).partial_cmp(&dot(w, w)).unwrap();
        if d <= Ordering::Equal {
            d
        } else {
            positive_direction(w).partial_cmp(&positive_direction(v)).unwrap()
        }
    });

    for w in &ws {
        let w: Vec<_> = w.iter().map(|x|
                if &abs(x.clone()) < epsilon { T::zero() } else { x.clone() }
            ).collect();
        a.push(if a.len() > 0 && dot(&a[0], &w).is_positive() {
            neg(&w)
        } else {
            w.clone()
        });

        if rank(&a) < a.len() {
            a.pop();
        } else if a.len() >= dim {
            return Some(a);
        }
    }

    None
}


pub fn dirichlet_vectors<T, F>(basis: &[Vec<T>], dot: F, epsilon: &T)
    -> Vec<Vec<T>>
    where T: Coord, for <'a> &'a T: CoordPtr<T>, F: Fn(&[T], &[T]) -> T
{
    let n = basis.len();
    assert!(basis.iter().all(|v| v.len() == n));

    match n {
        0 | 1 => basis.to_vec(),
        2 => {
            let (u, v) = lagrange_reduced(&basis[0], &basis[1], dot, epsilon);
            let w = add(&u, &v);
            vec![u, v, w]
        }
        3 => {
            let (u, v, w) = selling_reduced(
                &basis[0], &basis[1], &basis[2], dot, epsilon
            );
            let uv = add(&u, &v);
            let uw = add(&u, &w);
            let vw = add(&v, &w);
            let uvw = add(&uv, &w);
            vec![u, v, w, uv, uw, vw, uvw]
        }
        _ => panic!("dimension {} not supported", n)
    }
}


fn lagrange_reduced<T, F>(u: &[T], v: &[T], dot: F, epsilon: &T)
    -> (Vec<T>, Vec<T>)
    where T: Coord, for <'a> &'a T: CoordPtr<T>, F: Fn(&[T], &[T]) -> T
{
    let norm = |v: &[T]| dot(v, v);
    let fudge = &T::one() - epsilon;
    let mut u = u.to_vec();
    let mut v = v.to_vec();

    if norm(&v) < norm(&u) {
        (u, v) = (v, u);
    }

    while norm(&u) < &norm(&v) * &fudge {
        let t = dot(&u, &v).div_rounded(&norm(&u));
        let w = sub(&v, &mul(&u, &t));
        (u, v) = (w, u);
    }

    if dot(&u, &v).is_positive() {
        v = neg(&v);
    }

    (u, v)
}


fn selling_reduced<T, F>(u: &[T], v: &[T], w: &[T], dot: F, epsilon: &T)
    -> (Vec<T>, Vec<T>, Vec<T>)
    where T: Coord, for <'a> &'a T: CoordPtr<T>, F: Fn(&[T], &[T]) -> T
{
    let s = neg(&add(&add(&u, &v), &w));
    let mut vs = [u.to_vec(), v.to_vec(), w.to_vec(), s];

    loop {
        let mut changed = false;

        for i in 0..3 {
            for j in (i + 1)..4 {
                if &dot(&vs[i], &vs[j]) > epsilon {
                    for k in 0..4 {
                        if k != i && k != j {
                            vs[k] = add(&vs[k], &vs[i]);
                        }
                    }
                    vs[i] = neg(&vs[i]);
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


fn rank<T>(vs: &[Vec<T>]) -> usize
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    let mut basis = vec![];
    for v in vs {
        extend_basis(&v, &mut basis);
    }
    basis.len()
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
            return neg(v);
        }
    }
    v.to_vec()
}


fn neg<T>(v: &[T]) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    v.iter().map(|x| -x).collect()
}


fn add<T>(v: &[T], w: &[T]) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    (0..v.len()).map(|i| &v[i] + &w[i]).collect()
}


fn sub<T>(v: &[T], w: &[T]) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    (0..v.len()).map(|i| &v[i] - &w[i]).collect()
}


fn mul<T>(v: &[T], t: &T) -> Vec<T>
    where T: Coord, for <'a> &'a T: CoordPtr<T>
{
    (0..v.len()).map(|i| &v[i] * t).collect()
}


pub fn rational_lattice_basis(vectors: &Vec<Vec<BigRational>>)
    -> Vec<Vec<BigRational>>
{
    let dim = vectors[0].len();
    let mut common_denom = BigInt::one();

    for d in vectors {
        for x in d {
            common_denom *= x.denom() / common_denom.gcd(x.denom());
        }
    }

    let common_denom = common_denom.to_i32().unwrap();
    let f = BigRational::from(BigInt::from(common_denom));

    let mut basis = vec![];

    for k in 0..dim {
        let mut b = vec![0; dim];
        b[k] = common_denom;
        extend_basis(&b, &mut basis);
    }

    for d in vectors {
        let d: Vec<_> = d.iter().map(|x| r_to_i32(&(x * &f)).unwrap()).collect();
        extend_basis(&d, &mut basis);
    }

    basis.iter()
        .map(|b| b.iter()
            .map(|x| BigRational::from(BigInt::from(*x)) / &f)
            .collect::<Vec<_>>()
        )
        .collect()
}


fn r_to_i32(x: &BigRational) -> Option<i32> {
    if x.is_integer() {
        x.to_integer().to_i32()
    } else {
        None
    }
}


#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::*;

    fn r(x: i32) -> BigRational {
        BigRational::from(BigInt::from(x))
    }

    #[test]
    fn test_lattice_reduced2d() {
        let eps = r(0);
        let vs = [vec![r(3), r(2)], vec![r(4), r(3)]];

        let (u, v) = lagrange_reduced(&vs[0], &vs[1], dot, &eps);
        let w = neg(&add(&u, &v));
        assert!(!dot(&u, &v).is_positive());
        assert!(!dot(&u, &w).is_positive());
        assert!(!dot(&v, &w).is_positive());

        let ws = reduced_lattice_basis(&vs, dot, &eps).unwrap();
        let (u, v) = (ws[0].clone(), ws[1].clone());
        assert_eq!(rank(&vec![u.clone(), v.clone()]), 2);
        assert!(dot(&u, &u) <= dot(&v, &v));
        assert!(!dot(&u, &v).is_positive());
    }

    #[test]
    fn test_lattice_shift2d() {
        let eps = r(0);
        let p = vec![r(32) / r(10), -r(14) / r(10)];
        let vs = vec![vec![r(1), r(0)], vec![r(0), r(1)], vec![r(1), r(1)]];
        let s = shift_for_dirichlet_domain(&p, &vs, dot, &eps);
        let q = add(&p, &s);
        assert_eq!(&q, &vec![r(2) / r(10), -r(4) / r(10)]);
    }

    #[test]
    fn test_lattice_reduced3d() {
        let eps = r(0);
        let vs = [
            vec![r(3), r(2), r(5)],
            vec![r(4), r(3), r(2)],
            vec![r(7), r(5), r(8)],
        ];

        let (u, v, w) = selling_reduced(&vs[0], &vs[1], &vs[2], dot, &eps);
        let t = neg(&add(&add(&u, &v), &w));
        assert!(!dot(&u, &v).is_positive());
        assert!(!dot(&u, &w).is_positive());
        assert!(!dot(&u, &t).is_positive());
        assert!(!dot(&v, &w).is_positive());
        assert!(!dot(&v, &t).is_positive());
        assert!(!dot(&w, &t).is_positive());

        let ws = reduced_lattice_basis(&vs, dot, &eps).unwrap();
        let (u, v, w) = (ws[0].clone(), ws[1].clone(), ws[2].clone());
        assert_eq!(rank(&vec![u.clone(), v.clone(), w.clone()]), 3);
        assert!(dot(&u, &u) <= dot(&v, &v));
        assert!(dot(&v, &v) <= dot(&w, &w));
        assert!(!dot(&u, &v).is_positive());
        assert!(!dot(&u, &w).is_positive());
    }

    #[test]
    fn test_lattice_shift3d() {
        let eps = r(0);
        let p = vec![r(36) / r(10), -r(16) / r(10), r(24) / r(10)];
        let vs = vec![
            vec![r(1), r(0), r(0)],
            vec![r(0), r(1), r(0)],
            vec![r(0), r(0), r(1)],
            vec![r(1), r(1), r(0)],
            vec![r(1), r(0), r(1)],
            vec![r(0), r(1), r(1)],
            vec![r(1), r(1), r(1)],
        ];
        let s = shift_for_dirichlet_domain(&p, &vs, dot, &eps);
        let q = add(&p, &s);
        assert_eq!(&q, &vec![r(-4) / r(10), r(4) / r(10), r(4) / r(10)]);
    }
}
