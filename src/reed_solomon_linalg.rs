#![allow(non_snake_case)]
use crate::gf_256::GF;
use crate::linalg::Matrix;


const K: usize = 3;  // Number of data elements (before encoding)
const N: usize = 5;  // Number of encoded elements

//// We also need to specify the evaluation points. We need to have `N` of these.
//// For this variant of Reed-Solomon, it really doesn't matter what we pick. They just need to be
//// distinct elements of `GF(256)`. We choose: `[42, 222, 2, 8, 99]`.
const X: [GF; 5] = [
    GF::new(42), GF::new(222), GF::new(2), GF::new(8), GF::new(99)
];

//// We'll also need some y's for cauchy
const Y: [GF; 3] = [
    GF::new(66), GF::new(84), GF::new(112),
];

//// Build N x K vandermonde encoding matrix

pub fn vandermonde_encoding_matrix() -> Matrix {
    let mut out = Matrix::zeros(N, K);

    // for each row
    for i in 0..N {
        let x = X[i];
        let mut elt = GF::new(1);
        // for each col
        for j in 0..K {
            out[(i,j)] = elt;
            elt = elt * x;
        }
    }

    out
}

//// Let's write the encode routine. BLAH BLAH.
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

//// Let's write the decode routine.
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

//// This passes the exact same tests as the polynomial version!

#[cfg(test)]
mod test_vandermonde {
    use super::*;

    #[test]
    fn test_encode_decode() {
        let A = vandermonde_encoding_matrix();
        let data = vec![GF::new(100), GF::new(150), GF::new(200)];

        // encode
        let enc = reed_solomon_linalg_encode(&A, &data);
        assert_eq!(enc, vec![GF::new(160), GF::new(135), GF::new(94), GF::new(104), GF::new(194)]);
        let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

        // decode with no erasures
        assert_eq!(reed_solomon_linalg_decode(&A, &recv_all), Some(data.clone()));

        // erase element 0 and decode
        let mut recv = recv_all.clone();
        recv[0] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2 and decode
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 2,4 and decode
        let mut recv = recv_all.clone();
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2,3 and decode will fail
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        recv[3] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);

        // erase elements 0,2,4 and decode will fail
        let mut recv = recv_all.clone();
        recv[0] = None;
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);
    }
}

//// ### Systematic Reed-Solomon

//// We want systematic matrices now...

pub fn systematic_encoding_matrix() -> Matrix {
    // First build the normal vandermonde
    let A = vandermonde_encoding_matrix();

    // Second compute the inverse of the upper K x K square
    let inv = A.slice_rows(0..K).inverse().unwrap();

    // Third, Multiply
    A.matmul(&inv)
}

//// #### Testing again

#[cfg(test)]
mod test_systematic {
    use super::*;

    #[test]
    fn test_systematic_encoding() {
        let A = systematic_encoding_matrix();
        assert_eq!(A.slice_rows(0..K), Matrix::identity(3));
    }

    #[test]
    fn test_encode_decode() {
        let A = systematic_encoding_matrix();
        let data = vec![GF::new(100), GF::new(150), GF::new(200)];

        // encode
        let enc = reed_solomon_linalg_encode(&A, &data);
        assert_eq!(enc, vec![GF::new(100), GF::new(150), GF::new(200), GF::new(64), GF::new(57)]);
        let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

        // decode with no erasures
        assert_eq!(reed_solomon_linalg_decode(&A, &recv_all), Some(data.clone()));

        // erase element 0 and decode
        let mut recv = recv_all.clone();
        recv[0] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2 and decode
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 2,4 and decode
        let mut recv = recv_all.clone();
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2,3 and decode will fail
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        recv[3] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);

        // erase elements 0,2,4 and decode will fail
        let mut recv = recv_all.clone();
        recv[0] = None;
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);
    }
}

//// ### Reed-Solomon with Cauchy matrices

//// We want cauchy matrices now...

pub fn cauchy_encoding_matrix() -> Matrix {
    // Sanity check
    assert_eq!(X.len(), N);
    assert_eq!(Y.len(), K);

    let mut out = Matrix::zeros(N, K);

    // for each row
    for i in 0..N {
        let x = X[i];
        for j in 0..K {
            let y = Y[j];
            out[(i,j)] = GF::new(1)/(x+y);
        }
    }

    out
}

//// #### Testing again

#[cfg(test)]
mod test_cauchy {
    use super::*;

    #[test]
    fn test_encode_decode() {
        let A = cauchy_encoding_matrix();
        let data = vec![GF::new(100), GF::new(150), GF::new(200)];

        // encode
        let enc = reed_solomon_linalg_encode(&A, &data);
        assert_eq!(enc, vec![GF::new(7), GF::new(31), GF::new(153), GF::new(171), GF::new(22)]);
        let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

        // decode with no erasures
        assert_eq!(reed_solomon_linalg_decode(&A, &recv_all), Some(data.clone()));

        // erase element 0 and decode
        let mut recv = recv_all.clone();
        recv[0] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2 and decode
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 2,4 and decode
        let mut recv = recv_all.clone();
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), Some(data.clone()));

        // erase elements 1,2,3 and decode will fail
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        recv[3] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);

        // erase elements 0,2,4 and decode will fail
        let mut recv = recv_all.clone();
        recv[0] = None;
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_linalg_decode(&A, &recv), None);
    }
}
