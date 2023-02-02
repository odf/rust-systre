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
