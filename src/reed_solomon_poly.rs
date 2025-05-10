use crate::gf_256::GF;

//// #### Polynomial Evaluation
////
//// First, we need a routine to evaluate polynomials. Our routine will take an array of coefficients and
//// an evaluation point. We will use [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method)
//// for the algorithm. This allows us to do an evaluation in `O(n)` time.

fn poly_eval_point(coeffs: &[GF], x: GF) -> GF {
    let n = coeffs.len();
    let mut result = GF::new(0);
    for i in (0..n).rev() {
        result = result * x + coeffs[i];
    }
    result
}

//// We'll also add a routine to evaluate many points at once. This routine takes `O(nm)` time
//// where `n` is the number of coefficients and `m` is the number of evaluation points

fn poly_eval(coeffs: &[GF], xs: &[GF]) -> Vec<GF> {
    let mut out = vec![];
    for x in xs {
        out.push(poly_eval_point(coeffs, *x));
        }
    out
}

//// #### Polynomial Interpolation
////
//// For polynomial interpolation, we'll recover the coefficients by constructing a Newton Polynomial. But, unlike our
//// practice examples, we'll store all polynomials in monomial form.
////
//// Note: this algorithm is a bit tricky. If you don't follow it right
//// away, I'd suggest deriving it yourself directly and comparing results. If you
//// still are having difficulty, it's okay to skim this section. The important parts are:
//// (1) understanding that polynomial interpolation is not magic, and (2) understanding the
//// runtime complexity of the algorithm
////
//// Recall that a Newton Polynomial has the form:
////
//// ```
//// p(x) = k0 + k1*(x-x0) + k2*(x-x0)*(x-x1) + k3*(x-x0)*(x-x1)*(x-x2) + ...
//// ```
////
//// For our implementation, we'll maintain two polynomials:
////
//// 1. The `p(x)` polynomial constructed so far
//// 2. The `z(x)` polynomial of `(x - c)` factors so far
////
////
//// We can add new points by combining these polynomials (using a suitable `k`):
//// ```
//// p'(x) = p(x) + k*z(x)
//// ```
////
//// For each new point `(xk, yk)`, we will use the steps:
////
//// 1. Solve for `k` with `k = (yk - p(xk))/z(xk)`
//// 2. Update `p` with  `p = p + k * z`
//// 3. Update `z` with `z = z * (x - xk)`
////
//// We will also maintain the polynomials `p(x)` and `z(x)` in monomial coefficient form. This makes them
//// easy to combine and evaluate. It also means that the final answer is
//// just `p(x)` whenever we finish adding new points.
////
//// Finally, notice that this algorithm has runtime complexity `O(n^2)` where `n` is the number of points
//// interpolated. This is because the inner loop is `O(n)` and we're looping over `n` points.
////
//// Okay, let's code it:
fn poly_interpolate(evals: &[GF], xs: &[GF]) -> Vec<GF> {
    // Sanity check
    assert_eq!(evals.len(), xs.len());

    // Initialize our polynomials to p(x) = 0 and z(x) = 1
    let n = xs.len();
    let mut p = vec![GF::new(0); n];
    let mut z = vec![GF::new(0); n];
    z[0] = GF::new(1);

    // Loop over each point
    for i in 0..n {
        // Unpack the point for this iteration
        let xk = xs[i];
        let yk = evals[i];

        // Step 1: Solve k = (yk - p(xk)) / z(xk)
        let k = (yk - poly_eval_point(&p, xk)) / poly_eval_point(&z, xk);

        // Step 2: Update the polynomial: p' = p + k * z
        for j in 0..n {
            p[j] = p[j] + k * z[j];
        }

        // Step 3: Update the z polynomial: z' = z * (x - xk)
        //
        // Let's simplify a bit first by distributing:
        //   z' = z * x - z * xk
        //
        // Now notice that multiplying by `x` is a shift of the coefficients.
        // This leads to the form we'll use:
        //   z' = shift(z) - xk * z
        //
        for j in (1..n).rev() { // update in reverse
            z[j] = z[j-1] - xk * z[j];
        }
        z[0] = GF::new(0) - xk * z[0];
    }

    // The final result is just `p`
    p
}

//// #### Reed-Solomon Encoding and Decoding
//// Now, we'll implement the encoding and decoding functions. First, we'll to define a few parameters:
const K: usize = 3;  // Number of data elements (before encoding)
const N: usize = 5;  // Number of encoded elements

//// We also need to specify the evaluation points. We need to have `N` of these.
//// For this variant of Reed-Solomon, it really doesn't matter what we pick. They just need to be
//// distinct elements of `GF(256)`:
const X: [GF; 5] = [
    GF::new(42), GF::new(222), GF::new(2), GF::new(8), GF::new(99)
];

//// The encode routine is just a simple wrapper around polynomial evaluation
pub fn reed_solomon_poly_encode(data: &[GF]) -> Vec<GF> {
    // Sanity check
    assert_eq!(data.len(), K);
    assert_eq!(X.len(), N);

    // Just treat the data as if it's polynomial coefficients and evaluate!
    poly_eval(data, &X)
}

//// For decoding, we'll specify presence/erasure using the `Option<GF>` type. If present, we have `Some(x)` and otherwise `None`.
//// This is not a common way to do it, but it illustrates the idea of erasures very well.
pub fn reed_solomon_poly_decode(encoded: &[Option<GF>]) -> Option<Vec<GF>> {
    // Sanity check
    assert_eq!(encoded.len(), N);

    // First, we need to gather up evaluations that haven't been erased.
    let mut evals = vec![];
    let mut xs = vec![];
    for i in 0..N {
        if let Some(value) = encoded[i] {
            evals.push(value);
            xs.push(X[i]);
        }
    }

    // Make sure we have enough evaluations to decode
    if evals.len() < K {
        return None; // Too many erasures, can't decode!
    }

    // Decode it with polynomial interpolation. Note that we only use the
    // first K evaluations. The additional evaluations aren't needed (we are
    // assuming that non-erased evaluations are "error-free")
    let data = poly_interpolate(&evals[..K], &xs[..K]);

    // Great success!
    Some(data)
}

//// #### Testing Time
////
////
//// First, we'll test that polynomial evaluation works correctly:
#[cfg(test)]
#[test]
fn test_poly_eval() {
    assert_eq!(poly_eval_point(&[], GF::new(1)), GF::new(0));
    assert_eq!(poly_eval_point(&[GF::new(2)], GF::new(2)), GF::new(2));
    assert_eq!(poly_eval_point(&[GF::new(2), GF::new(2)], GF::new(3)), GF::new(4));
    assert_eq!(poly_eval_point(&[GF::new(46), GF::new(198), GF::new(0), GF::new(89)], GF::new(201)), GF::new(126));
}

//// And we'll test polynomial interpolation
#[cfg(test)]
#[test]
fn test_poly_interpolate() {
    let poly = poly_interpolate(
        &[GF::new(21), GF::new(21), GF::new(19)],
        &[GF::new(4), GF::new(5), GF::new(6)]
    );
    assert_eq!(poly, vec![GF::new(1), GF::new(1), GF::new(1)]);

    let poly = poly_interpolate(
        &[GF::new(57), GF::new(56), GF::new(49)],
        &[GF::new(4), GF::new(5), GF::new(6)]
    );
    assert_eq!(poly, vec![GF::new(1), GF::new(2), GF::new(3)]);
}

//// We should be able to decode up to 2-element erasures. But, 3-element erasures will fail.
////
//// Let's test each possibility.
#[cfg(test)]
fn encode_decode_all(data: &[GF], expected_enc: &[GF]) {
    // encode
    let enc = reed_solomon_poly_encode(data);
    assert_eq!(enc, expected_enc);
    let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

    // decode with no erasures: success!
    assert_eq!(reed_solomon_poly_decode(&recv_all), Some(data.to_vec()));

    // decode with all one element erasures: success!
    for i in 0..N {
        let mut recv = recv_all.clone();
        recv[i] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), Some(data.to_vec()));
    }

    // decode with all two element erasures: success!
    for i in 0..N {
        for j in (i+1)..N {
            let mut recv = recv_all.clone();
            recv[i] = None;
            recv[j] = None;
            assert_eq!(reed_solomon_poly_decode(&recv), Some(data.to_vec()));
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
                assert_eq!(reed_solomon_poly_decode(&recv), None);
            }
        }
    }
}

//// Now we can test a bunch of different data vectors:
#[cfg(test)]
#[test]
fn test_encode_decode() {
    // test: trivial
    encode_decode_all(
        &[GF::new(0), GF::new(0), GF::new(0)],
        &[GF::new(0), GF::new(0), GF::new(0), GF::new(0), GF::new(0)],
    );

    // test: ones
    encode_decode_all(
        &[GF::new(1), GF::new(1), GF::new(1)],
        &[GF::new(3), GF::new(161), GF::new(7), GF::new(73), GF::new(160)],
    );

    // test: pattern
    encode_decode_all(
        &[GF::new(100), GF::new(150), GF::new(200)],
        &[GF::new(160), GF::new(135), GF::new(94), GF::new(104), GF::new(194)],
    );

    // test: random
    encode_decode_all(
        &[GF::new(216), GF::new(196), GF::new(171)],
        &[GF::new(81), GF::new(157), GF::new(209), GF::new(193), GF::new(105)],
    );
}
