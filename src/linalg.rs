#![allow(non_snake_case)]
use crate::gf_256::GF;
use std::ops::{Index, IndexMut, Range};

//// #### Matrix structure

//// First, we will define a matrix type. Each matrix may have a different number of rows and columns. We'll store the matrix elements
//// in row-major ordering. We will use the shape parameters to calculate offsets into the flattened array.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    dat: Vec<GF>, // elements stored in row-major order
}

//// We'll have two primitive ways to construct a matrix: (1) from explicit data elements or
//// (2) initialized to all zeros
impl Matrix {
    pub fn new(rows: usize, cols: usize, dat: Vec<GF>) -> Self {
        assert_eq!(rows * cols, dat.len());
        Self { rows, cols, dat }
    }

    pub fn zeros(rows: usize, cols: usize) -> Self {
        let dat = vec![GF::new(0); rows * cols];
        Self { rows, cols, dat }
    }
}

//// We'll also need a few simple accessors:
impl Matrix {
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn data(&self) -> &[GF] {
        &self.dat
    }
}

//// #### Element Indexing

//// Let's now implement matrix element indexing.
////
//// An `m x n` matrix is indexed by a tuple `(i,j)` where
//// `0 <= i < m` and `0 <= j < n`.
////
//// Syntactically:
////
//// - `mat[(0,2)]` accesses the element in row 0 and column 2
//// - `mat[(5,0)]` accesses the element in row 5 and column 0
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

//// #### Identity Matrices

//// Now that we have indexing, we can write a routine to construct identity matrices. We'll start with
//// all zeros, and then write 1 to the diagonal entries:

impl Matrix {
    pub fn identity(n: usize) -> Self {
        let mut m = Matrix::zeros(n, n);
        for i in 0..n {
            m[(i,i)] = GF::new(1);
        }
        m
    }
}

//// #### Slicing
////
//// We'll occasionally need to slice rows and columns of a matrix to form a new sub-matrix.
////
//// Syntactically:
////
//// - `mat.slice_rows(0..3)` constructs a new matrix with rows 0 through 3 of the original matrix
//// - `mat.slice_cols(2..4)` constructs a new matrix with columns 2 through 4 of the original matrix

impl Matrix {
    pub fn slice_rows(&self, r: Range<usize>) -> Matrix {
        let start = r.start;
        let end = r.end;
        assert!(end <= self.rows);

        let new_rows = end - start;
        let mut out = Matrix::zeros(new_rows, self.cols);
        for i in 0..new_rows {
            for j in 0..self.cols {
                out[(i, j)] = self[(start+i, j)];
            }
        }

        out
    }

    pub fn slice_cols(&self, r: Range<usize>) -> Matrix {
        let start = r.start;
        let end = r.end;
        assert!(end <= self.cols);

        let new_cols = end - start;
        let mut out = Matrix::zeros(self.rows, new_cols);
        for i in 0..self.rows {
            for j in 0..new_cols {
                out[(i, j)] = self[(i, start+j)];
            }
        }

        out
    }
}

//// #### Selection
////
//// Similar to slicing, we'll also want to select subsets of rows and columns by index
////
//// Syntactically:
////
//// - `mat.select_rows(&[0,1,2,4])` constructs a new matrix with rows 0, 1, 2, and 4 of the original matrix
//// - `mat.select_cols(&[3,2])` constructs a new matrix with columns 3 and 2 of the original matrix

impl Matrix {
    pub fn select_rows(&self, sel: &[usize]) -> Matrix {
        let new_rows = sel.len();
        let mut out = Matrix::zeros(new_rows, self.cols);
        for i in 0..new_rows {
            for j in 0..self.cols {
                out[(i, j)] = self[(sel[i], j)];
            }
        }
        out
    }

    pub fn select_cols(&self, sel: &[usize]) -> Matrix {
        let new_cols = sel.len();
        let mut out = Matrix::zeros(self.rows, new_cols);
        for i in 0..self.rows {
            for j in 0..new_cols {
                out[(i, j)] = self[(i, sel[j])];
            }
        }
        out
    }
}

//// #### Matrix multiplication
//// A matrix multiply can be decomposed into a sequence of multiplications and additions, and any field has
//// these two operations. This means that we can perform matrix multiplication using field arithmetic!
////
//// The implementation isn't anything special either. It's just the textbook algorithm.
impl Matrix {
    pub fn matmul(&self, other: &Matrix) -> Matrix {
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

//// #### LU Factorization

//// Now, we wish to solve linear systems of equations in the form of `Ax = b`.  We can also perform this operation over an arbitrary field!

//// The classic approach is to factor the
//// `A` matrix into upper-triangular (`U`) and lower-triangular (`L`) parts. That is, we want to find
//// `A = LU`. This can be done with Gaussian Elimination.

//// First, we'll define a struct to store the results of a factorization. In optimized implementations,
//// the `A` matrix is typically mutated in-place to store the resulting `L` and `U` factors. But we'll expose
//// them explicitly here for ease of understanding.
////
//// We'll also need a permutation vector for row exchanges. This vector maps the new row (used for
//// pivoting) to the original row in the `A` matrix.

#[derive(Debug)]
pub struct LU {
    lower: Matrix,       // L: lower triangular matrix
    upper: Matrix,       // U: upper triangular matrix
    permute: Vec<usize>, // P: row permutation map (new_row => old_row)
}

//// Let's now implement Gaussian elimination. We will need to use partial pivoting because
//// non-singular matrices can have zeros on the diagonal. But, unlike numerical linear algebra, we
//// don't have to worry about numerical stability. Finite field arithmetic is exact. This is one of
//// the many cool things about linear algebra over finite fields.
////
//// Note, the algorithm below is completely "textbook". I'm not going to explain it in detail. If you have trouble
//// following, consider reviewing LU factorization over Real Numbers. The algorithm doesn't fundamentally change when
//// using it over a finite field.

impl Matrix {
    pub fn factor_lu(&self) -> Option<LU> {
        // Sanity check the shape (a square matrix required in this implementation)
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
            // have to worry about numerical issues. We only need to pivot because
            // non-singular matrices can still have zeros in the pivots. If the matrix
            // is non-singular, then the pivot column must have at least one non-zero
            // entry. This, partial pivoting is always sufficient!
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
            L[(P[k],k)] = GF::new(1);
            for i in (k+1)..rows {
                L[(P[i],k)] = A[(P[i], k)] / pivot;
            }

            // Apply the transform (row subtraction to the sub-matrix)
            for i in (k+1)..rows {
                let m = L[(P[i],k)];
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

//// #### Solving `Ax = b` using an LU Factorization

//// Once we have an LU Factorization, we can solve `Ax = b` as `LUx = b`. This is a much easier problem.
////
//// There are two steps:
////
//// 1. Solve `Ly = b` for `y` (using forward-substitution)
//// 2. Solve `Ux = y` for `x` (using back-substitution)
////
//// But, we'll need to be careful to use our permutation map to reorder the `b` vector in the forward-substitution step.
////
//// Let's start with forward-substitution. This is also completely "textbook".

fn solve_lower_triangular(L: &Matrix, b: &Matrix, permute: &[usize]) -> Matrix {
    // Sanity check the shape
    assert!(L.rows == L.cols);
    assert!(b.rows == L.rows);
    assert!(b.cols == 1);

    let n = L.rows;
    let mut b = b.clone();
    let mut y = Matrix::zeros(n, 1);
    for j in 0..n { // columns in L
        let elt = b[(permute[j],0)] / L[(permute[j],j)];
        y[(j,0)] = elt;

        for i in (j+1)..n { // walk down the column, subtracting off
            b[(permute[i],0)] = b[(permute[i], 0)] - L[(permute[i],j)] * elt;
        }
    }
    y
}

//// Next, we'll implement back-substitution. This is nearly the same as forward-substitution, but going
//// in the opposite direction.
fn solve_upper_triangular(U: &Matrix, y: &Matrix) -> Matrix {
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

//// With these two routines, we can now solve the full `LUx = b` problem by simply combining them.

impl LU {
    // Solve LUx = b
    pub fn solve(&self, b: &Matrix) -> Matrix {
        // Sanity check the shape
        assert!(self.lower.rows == b.rows);
        assert!(self.upper.rows == b.rows);
        assert!(self.upper.rows == self.upper.cols);
        assert!(self.lower.rows == self.lower.cols);
        assert!(b.cols == 1);

        // Ly = b: Forward solve with lower triangular matrix (using the permutation)
        let y = solve_lower_triangular(&self.lower, b, &self.permute);

        // Ux = y: Backward solve with upper triangular matrix
        let x = solve_upper_triangular(&self.upper, &y);

        // Done! The result is 'x'
        x
    }
}

//// And we can add a solver routine for `Ax = b` directly
impl Matrix {
    // Solve Ax = b with this matrix
    pub fn solve(&self, b: &Matrix) -> Option<Matrix> {
        Some(self.factor_lu()?.solve(b))
    }
}

//// #### Matrix Inverses
//// Solving linear systems is a core operation in linear algebra. We now have great power!
////
//// Let's use it to construct matrix inverses. This is as simple as solving the linear system `Ax = b` where
//// `b` is each column of the identity matrix (i.e. each basis vector).
impl Matrix {
    // Compute the inversion of the matrix. Returns None if the matrix is singular.
    pub fn inverse(&self) -> Option<Matrix> {
        // Factor A into LU
        let lu = self.factor_lu()?;

        // Allocate the output inverse matrix
        let n = self.rows;
        let mut inv = Matrix::zeros(n, n);

        // Loop over each column of the output inverse matrix
        for j in 0..n {
            // Construct the ith basis vector
            let mut b = Matrix::zeros(n, 1);
            b[(j,0)] = GF::new(1);

            // Solve that basis vector
            let x = lu.solve(&b);

            // Copy the result into the ith column of the inverse matrix
            for i in 0..n {
                inv[(i, j)] = x[(i, 0)];
            }
        }

        Some(inv)
    }
}

//// ### Testing time
////
//// We omit the test listing in this article for cleanliness. We encourage you to read them in the [source code](https://github.com/xorvoid/learn_you_galois_fields_for_great_good/blob/main/src/linalg.rs).
//// <!--
#[cfg(test)]
mod tests {
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
    fn test_slicing() {
        let m = Matrix::new(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);

        assert_eq!(m.slice_rows(0..5), m);
        assert_eq!(m.slice_cols(0..5), m);

        assert_eq!(m.slice_rows(1..4), Matrix::new(3, 5, vec![
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
        ]));

        assert_eq!(m.slice_cols(1..4), Matrix::new(5, 3, vec![
            GF::new(76), GF::new(234), GF::new(100),
            GF::new(123), GF::new(34), GF::new(2),
            GF::new(33), GF::new(203), GF::new(91),
            GF::new(16), GF::new(150), GF::new(160),
            GF::new(85), GF::new(191), GF::new(230),
        ]));

        assert_eq!(m.slice_rows(3..3), Matrix::new(0, 5, vec![]));
        assert_eq!(m.slice_cols(0..0), Matrix::new(5, 0, vec![]));
    }

    #[test]
    fn test_selecting() {
        let m = Matrix::new(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);

        assert_eq!(m.select_rows(&[0, 1, 2, 3, 4]), m);
        assert_eq!(m.select_cols(&[0, 1, 2, 3, 4]), m);

        assert_eq!(m.select_rows(&[1, 2, 4]), Matrix::new(3, 5, vec![
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]));

        assert_eq!(m.select_cols(&[1, 3, 4]), Matrix::new(5, 3, vec![
            GF::new(76),  GF::new(100), GF::new(1),
            GF::new(123), GF::new(2),   GF::new(12),
            GF::new(33),  GF::new(91),  GF::new(8),
            GF::new(16),  GF::new(160), GF::new(49),
            GF::new(85),  GF::new(230), GF::new(96),
        ]));

        assert_eq!(m.select_rows(&[]), Matrix::new(0, 5, vec![]));
        assert_eq!(m.select_cols(&[]), Matrix::new(5, 0, vec![]));
    }

    #[test]
    fn test_matrix_identity() {
        assert_eq!(Matrix::identity(2), Matrix::new(2, 2, vec![
            GF::new(1), GF::new(0),
            GF::new(0), GF::new(1),
        ]));

        assert_eq!(Matrix::identity(5), Matrix::new(5, 5, vec![
            GF::new(1), GF::new(0), GF::new(0), GF::new(0), GF::new(0),
            GF::new(0), GF::new(1), GF::new(0), GF::new(0), GF::new(0),
            GF::new(0), GF::new(0), GF::new(1), GF::new(0), GF::new(0),
            GF::new(0), GF::new(0), GF::new(0), GF::new(1), GF::new(0),
            GF::new(0), GF::new(0), GF::new(0), GF::new(0), GF::new(1),
        ]));
    }

    #[test]
    fn test_matrix_multiply() {
        let a = Matrix::new(3, 2, vec![
            GF::new(45), GF::new(89),
            GF::new(2), GF::new(123),
            GF::new(11), GF::new(200),
        ]);

        let b = Matrix::new(2, 3, vec![
            GF::new(8), GF::new(77), GF::new(211),
            GF::new(139), GF::new(3), GF::new(197),
        ]);

        let c_expected = Matrix::new(3, 3, vec![
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
    fn test_factor_lu() {
        let m = Matrix::new(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);
        let lu = m.factor_lu().unwrap();

        let L = &lu.lower;
        assert_eq!(L, &Matrix::new(5, 5, vec![
            GF::new(1), GF::new(0), GF::new(0), GF::new(0), GF::new(0),
            GF::new(179), GF::new(1), GF::new(0), GF::new(0), GF::new(0),
            GF::new(127), GF::new(137), GF::new(1), GF::new(0), GF::new(0),
            GF::new(131), GF::new(227), GF::new(135), GF::new(1), GF::new(0),
            GF::new(199), GF::new(90), GF::new(104), GF::new(106), GF::new(1)
        ]));

        let U = &lu.upper;
        assert_eq!(U, &Matrix::new(5, 5, vec![
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
        // LU factor a matrix that requires a row swap (non-trivial permutation)
        let m = Matrix::new(2, 2, vec![
            GF::new(0),  GF::new(5),
            GF::new(10), GF::new(20),
        ]);

        let lu = m.factor_lu().unwrap();

        let L = &lu.lower;
        assert_eq!(L, &Matrix::new(2, 2, vec![
            GF::new(0), GF::new(1),
            GF::new(1), GF::new(0),
        ]));

        let U = &lu.upper;
        assert_eq!(U, &Matrix::new(2, 2, vec![
            GF::new(10), GF::new(20),
            GF::new(0), GF::new(5),
        ]));

        let P = &lu.permute;
        assert_eq!(P, &[1, 0]);
    }

    #[test]
    fn test_solve_lower_triangular() {
        let L = Matrix::new(4, 4, vec![
            GF::new(100), GF::new(0),  GF::new(0),   GF::new(0),
            GF::new(23),  GF::new(8),  GF::new(0),   GF::new(0),
            GF::new(41),  GF::new(10), GF::new(33),  GF::new(0),
            GF::new(201), GF::new(55), GF::new(192), GF::new(111),
        ]);

        let x = Matrix::new(4, 1, vec![
            GF::new(73), GF::new(117), GF::new(202), GF::new(244),
        ]);

        let b = L.matmul(&x);
        assert_eq!(b, Matrix::new(4, 1, vec![
            GF::new(1), GF::new(157), GF::new(89), GF::new(35),
        ]));

        let xp = solve_lower_triangular(&L, &b, &[0, 1, 2, 3]);
        assert_eq!(x, xp);
    }

    #[test]
    fn test_solve_upper_triangular() {
        let U = Matrix::new(4, 4, vec![
            GF::new(201), GF::new(55), GF::new(192), GF::new(111),
            GF::new(0),   GF::new(41), GF::new(10),  GF::new(33),
            GF::new(0),   GF::new(0),  GF::new(23),  GF::new(8),
            GF::new(0),   GF::new(0),  GF::new(0),   GF::new(100),
        ]);

        let x = Matrix::new(4, 1, vec![
            GF::new(73), GF::new(117), GF::new(202), GF::new(244),
        ]);

        let b = U.matmul(&x);
        assert_eq!(b, Matrix::new(4, 1, vec![
            GF::new(35), GF::new(10), GF::new(181), GF::new(29),
        ]));

        let xp = solve_upper_triangular(&U, &b);
        assert_eq!(x, xp);
    }

    #[test]
    fn test_factor_lu_and_solve() {
        let m = Matrix::new(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);

        let x1 = Matrix::new(5, 1, vec![
            GF::new(67), GF::new(79), GF::new(119), GF::new(200), GF::new(151),
        ]);
        let b1 = m.matmul(&x1);

        let x2 = Matrix::new(5, 1, vec![
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
        // LU factor a matrix that requires a row swap (non-trivial permutation)
        let m = Matrix::new(2, 2, vec![
            GF::new(0),  GF::new(5),
            GF::new(10), GF::new(20),
        ]);

        let x1 = Matrix::new(2, 1, vec![
            GF::new(67), GF::new(79),
        ]);
        let b1 = m.matmul(&x1);

        let x2 = Matrix::new(2, 1, vec![
            GF::new(146), GF::new(37),
        ]);
        let b2 = m.matmul(&x2);

        let lu = m.factor_lu().unwrap();

        let x = lu.solve(&b1);
        assert_eq!(x, x1);

        let x = lu.solve(&b2);
        assert_eq!(x, x2);
    }

    #[test]
    fn test_inverse() {
        let m = Matrix::new(5, 5, vec![
            GF::new(45), GF::new(76), GF::new(234), GF::new(100), GF::new(1),
            GF::new(55), GF::new(123), GF::new(34), GF::new(2), GF::new(12),
            GF::new(52), GF::new(33), GF::new(203), GF::new(91), GF::new(8),
            GF::new(6), GF::new(16), GF::new(150), GF::new(160), GF::new(49),
            GF::new(7), GF::new(85), GF::new(191), GF::new(230), GF::new(96),
        ]);

        let m_inv = m.inverse().unwrap();
        assert_eq!(m_inv, Matrix::new(5, 5, vec![
            GF::new(99),  GF::new(183), GF::new(91),  GF::new(237), GF::new(66),
            GF::new(165), GF::new(209), GF::new(224), GF::new(40),  GF::new(12),
            GF::new(246), GF::new(244), GF::new(219), GF::new(115), GF::new(69),
            GF::new(205), GF::new(198), GF::new(141), GF::new(69),  GF::new(236),
            GF::new(210), GF::new(105), GF::new(34),  GF::new(160), GF::new(127),
        ]));

        assert_eq!(m.matmul(&m_inv), Matrix::identity(5));
    }

    #[test]
    fn test_inverse_permute() {
        let m = Matrix::new(2, 2, vec![
            GF::new(0),  GF::new(5),
            GF::new(10), GF::new(20),
        ]);

        let m_inv = m.inverse().unwrap();
        assert_eq!(m_inv, Matrix::new(2, 2, vec![
            GF::new(164), GF::new(41),
            GF::new(82),  GF::new(0),
        ]));

        assert_eq!(m.matmul(&m_inv), Matrix::identity(2));
    }

    #[test]
    fn test_inverse_singular() {
        let m = Matrix::new(3, 3, vec![
            GF::new(45), GF::new(76),  GF::new(234),
            GF::new(55), GF::new(123), GF::new(34),
            GF::new(26), GF::new(55),  GF::new(200),  // This row is the sum of the previous two
        ]);

        assert!(m.inverse().is_none());
    }
}
//// -->
