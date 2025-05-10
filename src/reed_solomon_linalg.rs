#![allow(non_snake_case)]
#![allow(dead_code)]
use crate::gf_256::GF;
use crate::linalg::Matrix;

//// As before, we need some parameters:
const K: usize = 3;  // Number of data elements (before encoding)
const N: usize = 5;  // Number of encoded elements

//// We also need to specify the evaluation points. We need to have `N` of these.
//// For this variant of Reed-Solomon, it really doesn't matter what we pick. They just need to be
//// distinct elements of `GF(256)`.
const X: [GF; 5] = [
    GF::new(42), GF::new(222), GF::new(2), GF::new(8), GF::new(99)
];


//// #### Encoding and Decoding with a Matrix

//// For encoding, we'll use a generic matrix `A`. This matrix should have the required
//// invertibility properties.
////
//// We just need to perform a simple multiply (`e = A*d`):
pub fn reed_solomon_linalg_encode(A: &Matrix, data: &[GF]) -> Vec<GF> {
    // Sanity check
    assert_eq!(data.len(), K);

    // Construct the encoding matrix (A) and data vector (d)
    let d = Matrix::new(K, 1, data.to_vec());

    // Encode: e = A*d
    let e = A.matmul(&d);

    // Return the raw data
    e.data().to_vec()
}

//// Decoding will use the same `A` matrix, removing any rows for erased elements.
////
//// We just need to solve a linear system (solve `A_r * d = e_r` for `d`):
pub fn reed_solomon_linalg_decode(A: &Matrix, encoded: &[Option<GF>]) -> Option<Vec<GF>> {
    // Sanity check
    assert_eq!(encoded.len(), N);

    // First, we need to gather up evaluations that haven't been erased.
    let mut evals = vec![];
    let mut idx = vec![];
    for i in 0..N {
        if let Some(value) = encoded[i] {
            evals.push(value);
            idx.push(i);
        }
    }

    // Make sure we have enough evaluations to decode
    if evals.len() < K {
        return None; // Too many erasures, can't decode!
    }

    // Select the rows from the
    // Construct the matrix (A_r)
    let A_r = A.select_rows(&idx[..K]);
    let e_r = Matrix::new(K, 1, evals[..K].to_vec());

    // Decode by linear solve: A_r * d = e_r
    let d = A_r.solve(&e_r).unwrap();

    // Return the raw data
    Some(d.data().to_vec())
}

//// #### Vandermonde Matrices

//// We'll build a Vandermonde Matrix out of `N` evaluation points (`xs`)
pub fn vandermonde_encoding_matrix(xs: &[GF]) -> Matrix {
    // Sanity check
    assert_eq!(xs.len(), N);

    // The resulting matrix will be an N x K Vandermonde
    let mut mat = Matrix::zeros(N, K);

    // Build up one row at a time
    for i in 0..N {
        // Compute each value in the row: 1, x, x^2, x^3, etc
        let mut value = GF::new(1);
        for j in 0..K {
            mat[(i,j)] = value;
            value = value * xs[i];
        }
    }

    mat
}

//// #### Testing Time
////
//// Just like the polynomial version, we should be able to decode up to 2-element erasures.
//// But, 3-element erasures will fail.
////
//// Let's test each possibility.
#[cfg(test)]
fn encode_decode_all(A: &Matrix, data: &[GF], expected_enc: &[GF]) {

    // encode
    let enc = reed_solomon_linalg_encode(A, data);
    assert_eq!(enc, expected_enc);
    let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

    // decode with no erasures: success!
    assert_eq!(reed_solomon_linalg_decode(A, &recv_all), Some(data.to_vec()));

    // decode with all one element erasures: success!
    for i in 0..N {
        let mut recv = recv_all.clone();
        recv[i] = None;
        assert_eq!(reed_solomon_linalg_decode(A, &recv), Some(data.to_vec()));
    }

    // decode with all two element erasures: success!
    for i in 0..N {
        for j in (i+1)..N {
            let mut recv = recv_all.clone();
            recv[i] = None;
            recv[j] = None;
            assert_eq!(reed_solomon_linalg_decode(A, &recv), Some(data.to_vec()));
        }
    }

    // decode with all three element erasures: failure!
    for i in 0..N {
        for j in (i+1)..N {
            for k in (j+1)..N {
                let mut recv = recv_all.clone();
                recv[i] = None;
                recv[j] = None;
                recv[k] = None;
                assert_eq!(reed_solomon_linalg_decode(A, &recv), None);
            }
        }
    }
}

//// Now we can test a bunch of different data vectors. Note: these are the exact same ones
//// as in `reed_solomon_poly.rs`, giving an empirical example of equivalence.
#[cfg(test)]
#[test]
fn test_vandermonde_encode_decode() {
    // construct our vandermonde encoding matrix
    let A = vandermonde_encoding_matrix(&X);

    // test: trivial
    encode_decode_all(&A,
        &[GF::new(0), GF::new(0), GF::new(0)],
        &[GF::new(0), GF::new(0), GF::new(0), GF::new(0), GF::new(0)],
    );

    // test: ones
    encode_decode_all(&A,
        &[GF::new(1), GF::new(1), GF::new(1)],
        &[GF::new(3), GF::new(161), GF::new(7), GF::new(73), GF::new(160)],
    );

    // test: pattern
    encode_decode_all(&A,
        &[GF::new(100), GF::new(150), GF::new(200)],
        &[GF::new(160), GF::new(135), GF::new(94), GF::new(104), GF::new(194)],
    );

    // test: random
    encode_decode_all(&A,
        &[GF::new(216), GF::new(196), GF::new(171)],
        &[GF::new(81), GF::new(157), GF::new(209), GF::new(193), GF::new(105)],
    );
}

//// #### Systematic Matrices

//// We will now construct systematic matrices by transforming an existing `A` matrix
//// using the `A * inv(A_r)` approach:

pub fn systematic_encoding_matrix(A: &Matrix) -> Matrix {
    // Compute the inverse of the upper K x K square
    let inv = A.slice_rows(0..K).inverse().unwrap();

    // Multiply it
    A.matmul(&inv)
}

//// #### Testing again

//// We can completely reuse the encoding and decoding routines for the systematic matrices.
//// This is a very inefficient approach. In practice, you'd instead develop hardcoded
//// routines for the specific parity rows.
//// But for our pedagogical purposes, this approach is okay.
////
//// Let's demonstrate these with a test, just like before.

#[cfg(test)]
#[test]
fn test_systematic_encode_decode() {
    // construct our systematic matrix by transforming a vandermonde matrix
    let A = systematic_encoding_matrix(&vandermonde_encoding_matrix(&X));

    // The A matrix should start with a K x K identity
    assert_eq!(A.slice_rows(0..K), Matrix::identity(K));

    // The remaining N-K rows should be parity
    assert_eq!(A.slice_rows(K..N), Matrix::new(N-K, K, vec![
        GF::new(146), GF::new(30),  GF::new(141),
        GF::new(155), GF::new(137), GF::new(19),
    ]));

    // test: trivial
    encode_decode_all(&A,
        &[GF::new(0), GF::new(0), GF::new(0)],
        &[GF::new(0), GF::new(0), GF::new(0), GF::new(0), GF::new(0)],
    );

    // test: ones
    encode_decode_all(&A,
        &[GF::new(1), GF::new(1), GF::new(1)],
        &[GF::new(1), GF::new(1), GF::new(1), GF::new(1), GF::new(1)],
    );

    // test: pattern
    encode_decode_all(&A,
        &[GF::new(100), GF::new(150), GF::new(200)],
        &[GF::new(100), GF::new(150), GF::new(200), GF::new(64), GF::new(57)],
    );

    // test: random
    encode_decode_all(&A,
        &[GF::new(216), GF::new(196), GF::new(171)],
        &[GF::new(216), GF::new(196), GF::new(171), GF::new(31), GF::new(66)],
    );
}
