use num_bigint::BigInt;
use num_rational::BigRational;


pub fn modular_inverse(argument: i64, modulus: i64) -> Option<i64> {
    let (mut t, mut t1) = (0, 1);
    let (mut r, mut r1) = (modulus, argument);

    while r1 != 0 {
        let q = r / r1;
        (t, t1) = (t1, t - q * t1);
        (r, r1) = (r1, r - q * r1);
    }

    if r == 1 {
        Some(if t < 0 { t + modulus } else { t })
    } else {
        None
    }
}


type Matrix = Vec<Vec<i64>>;

pub fn modular_row_echelon_form(matrix: &Matrix, prime: i64) -> Matrix {
    if matrix.len() == 0 {
        return matrix.clone();
    }

    let mut a: Matrix = vec![];

    for row in matrix {
        assert!(a.len() == 0 || row.len() == a[0].len());
        let row = row.iter()
            .map(|&n| n % prime + if n < 0 { prime } else { 0 })
            .collect();
        a.push(row);
    }

    let (nrows, ncols) = (a.len(), a[0].len());

    let mut irow = 0;
    for icol in 0..ncols {
        if let Some(r) = (irow..nrows).find(|&r| a[r][icol] != 0) {
            if r != irow {
                swap_rows(&mut a, irow, r);
            }

            let f = modular_inverse(a[irow][icol], prime).unwrap();
            for j in icol..ncols {
                a[irow][j] = (a[irow][j] * f) % prime;
            }

            for i in 0..nrows {
                if i != irow && a[i][icol] != 0 {
                    let f = a[i][icol];
                    for j in icol..ncols {
                        if a[irow][j] != 0 {
                            a[i][j] = (
                                prime - (a[irow][j] * f) % prime + a[i][j]
                            ) % prime;
                        }
                    }
                }
            }
        }

        irow += 1;
    }

    a
}


fn swap_rows<T>(a: &mut Vec<Vec<T>>, i: usize, j: usize) {
    let row_i = std::mem::take(&mut a[i]);
    let row_j = std::mem::take(&mut a[j]);
    a[i] = row_j;
    a[j] = row_i;
}


pub fn modular_matrix_inverse(matrix: &Matrix, prime: i64) -> Option<Matrix> {
    let n = matrix.len();

    if n == 0 {
        return None;
    }

    let mut a = vec![];
    for (i, row) in matrix.iter().enumerate() {
        assert!(row.len() == n);
        let mut row = row.clone();
        row.extend(std::iter::repeat(0).take(n));
        row[n + i] = 1;
        a.push(row);
    }

    let a = modular_row_echelon_form(&a, prime);

    for i in 0..n {
        for j in 0..n {
            if a[i][j] != (i == j) as i64 {
                return None;
            }
        }
    }

    Some(a.iter().map(|row| row[n..].to_vec()).collect())
}


pub fn modular_matrix_product(a: &Matrix, b: &Matrix, prime: i64) -> Matrix {
    let ncols_a = a[0].len();
    let ncols_b = b[0].len();

    assert!(ncols_a == b.len());

    let mut product = vec![];
    for row in a {
        let mut row_out = vec![];
        for j in 0..ncols_b {
            let p = (0..ncols_a)
                .map(|i| (row[i] * b[i][j]) % prime)
                .reduce(|p, q| (p + q) % prime)
                .unwrap();
            row_out.push(p);
        }
        product.push(row_out);
    }

    product
}


pub fn integer_matrix_product(a: &Matrix, b: &Matrix) -> Matrix {
    let ncols_a = a[0].len();
    let ncols_b = b[0].len();

    assert!(ncols_a == b.len());

    let mut product = vec![];
    for row in a {
        let mut row_out = vec![];
        for j in 0..ncols_b {
            let p = (0..ncols_a).map(|i| row[i] * b[i][j]).sum();
            row_out.push(p);
        }
        product.push(row_out);
    }

    product
}


fn number_of_p_adic_steps_needed(a: &Matrix, b: &Matrix, prime: i64) -> u64
{
    let mut log_norms: Vec<_> = (0..a[0].len())
        .map(|j| column_norm(a, j).ln())
        .collect();
    log_norms.push((0..b[0].len())
        .map(|j| column_norm(b, j).ln())
        .max_by(|a, b| a.total_cmp(b))
        .unwrap());

    log_norms.sort_by(|a, b| a.total_cmp(b));

    let log_delta: f64 = log_norms.iter().skip(1).sum();
    let golden_ratio = (1.0 + (5 as f64).sqrt()) / 2.0;

    (2.0 * (log_delta + golden_ratio.ln()) / (prime as f64).ln()).ceil() as u64
}


fn column_norm(a: &Matrix, j: usize) -> f64 {
    (0..a.len()).map(|i| (a[i][j] as f64).powf(2.0)).sum::<f64>().sqrt()
}


fn rational_reconstruction(s: &BigInt, h: &BigInt) -> BigRational {
    let (mut u, mut u1) = (h.clone(), s.clone());
    let (mut v, mut v1) = (BigInt::from(0), BigInt::from(1));
    let mut sign = BigInt::from(1);

    while &u1.pow(2) > &h {
        let (q, r) = (&u / &u1, &u % &u1);

        (u, u1) = (u1, r);

        let v1_next = v + q * &v1;
        (v, v1) = (v1, v1_next);

        sign = -sign;
    }

    BigRational::new(sign * u1, v1)
}


pub fn solve(a: &Matrix, b: &Matrix) -> Option<Vec<Vec<BigRational>>> {
    let prime = 9999991;

    if let Some(c) = modular_matrix_inverse(a, prime) {
        let nr_steps = number_of_p_adic_steps_needed(a, b, prime);
        let nrows = b.len();
        let ncols = b[0].len();

        let mut p = BigInt::from(1);
        let mut b = b.clone();

        let mut s: Vec<Vec<_>> = (0..nrows).map(|_| {
            (0..ncols).map(|_| { BigInt::from(0) }).collect()
        }).collect();

        for step in 0..nr_steps {
            let x = modular_matrix_product(&c, &b, prime);
            for i in 0..nrows {
                for j in 0..ncols {
                    s[i][j] += &p * x[i][j];
                }
            }

            p *= prime;

            if step + 1 < nr_steps {
                let ax = integer_matrix_product(a, &x);
                for i in 0..nrows {
                    for j in 0..ncols {
                        b[i][j] = (b[i][j] - ax[i][j]) / prime;
                    }
                }
            }
        }

        Some(s.iter().map(|row| {
            row.iter().map(|x| { rational_reconstruction(x, &p) }).collect()
        }).collect())
    } else {
        None
    }
}
