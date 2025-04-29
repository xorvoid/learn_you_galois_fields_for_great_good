#![allow(non_snake_case)]
use crate::gf_256::GF;
use std::ops::{Index, IndexMut};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    dat: Vec<GF>, // row-major-order
}

impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        let mut dat = vec![];
        dat.resize(rows * cols, GF::new(0));
        Self { rows, cols, dat }
    }

    pub fn from_data(rows: usize, cols: usize, dat: Vec<GF>) -> Self {
        assert_eq!(rows * cols, dat.len());
        Self { rows, cols, dat }
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = GF;
    fn index(&self, idx: (usize, usize)) -> &GF {
        assert!(idx.0 < self.rows);
        assert!(idx.1 < self.cols);
        &self.dat[idx.0 * self.cols + idx.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut GF {
        assert!(idx.0 < self.rows);
        assert!(idx.1 < self.cols);
        &mut self.dat[idx.0 * self.cols + idx.1]
    }
}

impl Matrix {
    fn matmul(&self, other: &Matrix) -> Matrix {
        // Multiply two matrices: (m x n) * (n x p) => (m x p)
        let m = self.rows;
        let n = self.cols;
        let p = other.cols;
        assert_eq!(self.cols, other.rows);

        let mut out = Matrix::zeros(m, p);
        for i in 0..m {
            for j in 0..p {
                let mut elem = GF::new(0);
                for k in 0..n {
                    elem = elem + self[(i,k)] * other[(k,j)];
                }
                out[(i,j)] = elem;
            }
        }

        out
    }
}

impl Matrix {
    fn vandermonde(n: usize, g: GF) -> Matrix {
        let mut out = Matrix::zeros(n, n);
        let mut a = g;

        // for each row
        for i in 0..n {
            let mut elt = GF::new(1);
            // for each col
            for j in 0..n {
                out[(i,j)] = elt;
                elt = elt * a;
            }
            a = a * g;
        }

        out
    }
}

// LU Factorization
struct LU {
    lower: Matrix,
    upper: Matrix,
    permute: Vec<usize>, // row permutation map (new_row => old_row)
}

impl Matrix {
    fn factor_lu(&self) -> Option<LU> {
        // Sanity check the shape (square matrix required)
        assert_eq!(self.rows, self.cols);

        let mut A = self.clone();
        let rows = A.rows;
        let cols = A.cols;
        let mut L = Matrix::zeros(rows, cols);
        let mut U = Matrix::zeros(rows, cols);
        let mut P: Vec<_> = (0..rows).collect();

        // Loop over columns
        for k in 0..cols {
            // Search for a non-zero pivot. Unlike floating-point operations, we don't
            // have to worry about numerical issues. We only need to pivot because non-singular matrices
            // can still have zeros in the pivots. If the matrix is non-singular, then the pivot column
            // must have at least one non-zero entry. So, partial pivoting is always sufficent.
            let mut found = false;
            for i in k..rows {
                if A[(P[i], k)] != GF::new(0) {
                    // Perform a swap of the permutation map
                    let save = P[k];
                    P[k] = P[i];
                    P[i] = save;

                    // All good!
                    found = true;
                    break;
                }
            }
            if !found {
                return None; // Matrix is singular!
            }

            // Retrieve the pivot element
            let pivot = A[(P[k], k)];

            // Copy the pivot row to the U matrix
            for j in k..cols {
                U[(k,j)] = A[(P[k],j)];
            }


            // Compute the multipliers and store in column of L
            L[(k,k)] = GF::new(1);
            for i in (k+1)..rows {
                L[(i,k)] = A[(P[i], k)] / pivot;
            }

            // Apply the transform (row subtraction to the submatrix
            for i in (k+1)..rows {
                let m = L[(i,k)];
                for j in (k+1)..cols {
                    A[(P[i],j)] = A[(P[i],j)] - m * A[(P[k],j)];
                }
            }
        }

        Some(LU {
            lower: L,
            upper: U,
            permute: P,
        })
    }
}

// Ly = b: Forward solve with lower triangular matrix (using the permutation)
fn solve_lower_tringular(L: &Matrix, b: &Matrix, permute: &[usize]) -> Matrix {
    // Sanity check the shape
    assert!(L.rows == L.cols);
    assert!(b.rows == L.rows);
    assert!(b.cols == 1);

    let n = L.rows;
    let mut b = b.clone();
    let mut y = Matrix::zeros(n, 1);
    for j in 0..n { // columns in L
        let elt = b[(permute[j],0)] / L[(j,j)];
        y[(j,0)] = elt;

        for i in (j+1)..n { // walk down the column, subtracting off
            b[(permute[i],0)] = b[(permute[i], 0)] - L[(i,j)] * elt;
        }
    }
    y
}

// Ux = y: Backwards solve with upper triangular matrix
fn solve_upper_tringular(U: &Matrix, y: &Matrix) -> Matrix {
    // Sanity check the shape
    assert!(U.rows == U.cols);
    assert!(y.rows == U.rows);
    assert!(y.cols == 1);

    let n = U.rows;
    let mut y = y.clone();
    let mut x = Matrix::zeros(n, 1);
    for j in (0..n).rev() { // columns in L
        let elt = y[(j,0)] / U[(j,j)];
        x[(j,0)] = elt;

        for i in (0..j).rev() { // walk up the columns, subtracting off
            y[(i,0)] = y[(i, 0)] - U[(i,j)] * elt;
        }
    }
    x
}

impl LU {
    // Solve LUx = b
    fn solve(&self, b: &Matrix) -> Matrix {
        // Sanity check the shape
        assert!(self.lower.rows == b.rows);
        assert!(self.upper.rows == b.rows);
        assert!(self.upper.rows == self.upper.cols);
        assert!(self.lower.rows == self.lower.cols);
        assert!(b.cols == 1);
        let n = self.lower.rows;

        // Ly = b: Forward solve with lower triangular matrix (using the permutation)
        let y = solve_lower_tringular(&self.lower, b, &self.permute);

        // Ux = y: Backward solve with upper triangular matrix
        let x = solve_upper_tringular(&self.upper, &y);

        x
    }
}

impl Matrix {
    // Solve Ax = b with this matrix
    fn solve(&self, b: &Matrix) -> Option<Matrix> {
        let lu = self.factor_lu()?;
        Some(lu.solve(b))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_matrix_zeros() {
        assert_eq!(Matrix::zeros(1, 2).dat.len(), 2);
        assert_eq!(Matrix::zeros(2, 1).dat.len(), 2);
        assert_eq!(Matrix::zeros(4, 5).dat.len(), 20);
        assert_eq!(Matrix::zeros(9, 3).dat.len(), 27);
    }

    #[test]
    fn test_matrix_indexing() {
        let mut mat = Matrix::zeros(2, 3);
        mat[(0, 1)] = GF::new(5);
        mat[(1, 2)] = GF::new(1);
        assert_eq!(mat.dat, vec![
            GF::new(0), GF::new(5), GF::new(0),
            GF::new(0), GF::new(0), GF::new(1),
        ]);

        assert_eq!(mat[(0,0)], GF::new(0));
        assert_eq!(mat[(0,1)], GF::new(5));
        assert_eq!(mat[(0,2)], GF::new(0));
        assert_eq!(mat[(1,0)], GF::new(0));
        assert_eq!(mat[(1,1)], GF::new(0));
        assert_eq!(mat[(1,2)], GF::new(1));
    }

    #[test]
    fn test_matrix_multiply() {
        let mut a = Matrix::from_data(3, 2, vec![
            GF::new(45), GF::new(89),
            GF::new(2), GF::new(123),
            GF::new(11), GF::new(200),
        ]);

        let mut b = Matrix::from_data(2, 3, vec![
            GF::new(8), GF::new(77), GF::new(211),
            GF::new(139), GF::new(3), GF::new(197),
        ]);

        let mut c_expected = Matrix::from_data(3, 3, vec![
            GF::new(31), GF::new(180), GF::new(187),
            GF::new(161), GF::new(23), GF::new(17),
            GF::new(186), GF::new(202), GF::new(2),
        ]);

        let c = a.matmul(&b);
        assert_eq!(c.rows, c_expected.rows);
        assert_eq!(c.cols, c_expected.cols);
        assert_eq!(c.dat, c_expected.dat);
    }

    #[test]
    fn test_vandermonde() {
        let m = Matrix::vandermonde(5, GF::new(2));
        assert_eq!(m.rows, 5);
        assert_eq!(m.cols, 5);
        assert_eq!(m.dat, vec![
            GF::new(1), GF::new(2),  GF::new(4),   GF::new(8),   GF::new(16),
            GF::new(1), GF::new(4),  GF::new(16),  GF::new(64),  GF::new(27),
            GF::new(1), GF::new(8),  GF::new(64),  GF::new(54),  GF::new(171),
            GF::new(1), GF::new(16), GF::new(27),  GF::new(171), GF::new(94),
            GF::new(1), GF::new(32), GF::new(108), GF::new(47),  GF::new(151)
        ]);
    }

    #[test]
    fn test_factor_lu() {
        let m = Matrix::from_data(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);
        let lu = m.factor_lu().unwrap();

        let L = &lu.lower;
        assert_eq!(L, &Matrix::from_data(5, 5, vec![
            GF::new(1), GF::new(0), GF::new(0), GF::new(0), GF::new(0),
            GF::new(179), GF::new(1), GF::new(0), GF::new(0), GF::new(0),
            GF::new(127), GF::new(137), GF::new(1), GF::new(0), GF::new(0),
            GF::new(131), GF::new(227), GF::new(135), GF::new(1), GF::new(0),
            GF::new(199), GF::new(90), GF::new(104), GF::new(106), GF::new(1)
        ]));

        let U = &lu.upper;
        assert_eq!(U, &Matrix::from_data(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(0),  GF::new(87), GF::new(106), GF::new(80),  GF::new(191),
            GF::new(0),  GF::new(0),  GF::new(18),  GF::new(68),  GF::new(111),
            GF::new(0),  GF::new(0),  GF::new(0),   GF::new(44),  GF::new(111),
            GF::new(0),  GF::new(0),  GF::new(0),   GF::new(0),   GF::new(130)
        ]));

        let P = &lu.permute;
        assert_eq!(P, &[0, 1, 2, 3, 4]);

        assert_eq!(L.matmul(U), m);
    }

    #[test]
    fn test_factor_lu_permute() {
        // LU factor a matrix that requires a pivot
        let m = Matrix::from_data(2, 2, vec![
            GF::new(0),  GF::new(5),
            GF::new(10), GF::new(20),
        ]);

        let lu = m.factor_lu().unwrap();

        let L = &lu.lower;
        assert_eq!(L, &Matrix::from_data(2, 2, vec![
            GF::new(1), GF::new(0),
            GF::new(0), GF::new(1),
        ]));

        let U = &lu.upper;
        assert_eq!(U, &Matrix::from_data(2, 2, vec![
            GF::new(10), GF::new(20),
            GF::new(0), GF::new(5),
        ]));

        let P = &lu.permute;
        assert_eq!(P, &[1, 0]);
    }

    #[test]
    fn test_solve_lower_triangular() {
        let L = Matrix::from_data(4, 4, vec![
            GF::new(100), GF::new(0),  GF::new(0),   GF::new(0),
            GF::new(23),  GF::new(8),  GF::new(0),   GF::new(0),
            GF::new(41),  GF::new(10), GF::new(33),  GF::new(0),
            GF::new(201), GF::new(55), GF::new(192), GF::new(111),
        ]);

        let x = Matrix::from_data(4, 1, vec![
            GF::new(73), GF::new(117), GF::new(202), GF::new(244),
        ]);

        let b = L.matmul(&x);
        assert_eq!(b, Matrix::from_data(4, 1, vec![
            GF::new(1), GF::new(157), GF::new(89), GF::new(35),
        ]));

        let xp = solve_lower_tringular(&L, &b, &[0, 1, 2, 3]);
        assert_eq!(x, xp);
    }

    #[test]
    fn test_solve_upper_triangular() {
        let U = Matrix::from_data(4, 4, vec![
            GF::new(201), GF::new(55), GF::new(192), GF::new(111),
            GF::new(0),   GF::new(41), GF::new(10),  GF::new(33),
            GF::new(0),   GF::new(0),  GF::new(23),  GF::new(8),
            GF::new(0),   GF::new(0),  GF::new(0),   GF::new(100),
        ]);

        let x = Matrix::from_data(4, 1, vec![
            GF::new(73), GF::new(117), GF::new(202), GF::new(244),
        ]);

        let b = U.matmul(&x);
        assert_eq!(b, Matrix::from_data(4, 1, vec![
            GF::new(35), GF::new(10), GF::new(181), GF::new(29),
        ]));

        let xp = solve_upper_tringular(&U, &b);
        assert_eq!(x, xp);
    }


    #[test]
    fn test_factor_lu_and_solve() {
        let m = Matrix::from_data(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);

        let x1 = Matrix::from_data(5, 1, vec![
            GF::new(67), GF::new(79), GF::new(119), GF::new(200), GF::new(151),
        ]);
        let b1 = m.matmul(&x1);

        let x2 = Matrix::from_data(5, 1, vec![
            GF::new(146), GF::new(37), GF::new(155), GF::new(190), GF::new(131),
        ]);
        let b2 = m.matmul(&x2);

        let lu = m.factor_lu().unwrap();

        let x = lu.solve(&b1);
        assert_eq!(x, x1);

        let x = lu.solve(&b2);
        assert_eq!(x, x2);
    }

    #[test]
    fn test_factor_lu_and_solve_permute() {
        // LU factor a matrix that requires a pivot
        let m = Matrix::from_data(2, 2, vec![
            GF::new(0),  GF::new(5),
            GF::new(10), GF::new(20),
        ]);

        let lu = m.factor_lu().unwrap();

        let x1 = Matrix::from_data(2, 1, vec![
            GF::new(67), GF::new(79),
        ]);
        let b1 = m.matmul(&x1);

        let x2 = Matrix::from_data(2, 1, vec![
            GF::new(146), GF::new(37),
        ]);
        let b2 = m.matmul(&x2);

        let lu = m.factor_lu().unwrap();

        let x = lu.solve(&b1);
        assert_eq!(x, x1);

        let x = lu.solve(&b2);
        assert_eq!(x, x2);
    }
}
