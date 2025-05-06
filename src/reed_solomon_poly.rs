use crate::gf_256::GF;

//// #### Polynomial Evaluation
////
//// We will need a routine to evaluate polynomials. For this we'll accept an array of coefficents and
//// an evalation point. We will use [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method)
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
//// where `n` is the number of coefficents and `m` is the number of evaluation points

fn poly_eval(coeffs: &[GF], xs: &[GF]) -> Vec<GF> {
    let mut out = vec![];
    for x in xs {
        out.push(poly_eval_point(coeffs, *x));
        }
    out
}

//// #### Polynomial Interpolation
////
//// For polynomial interpolation, we'll recover the coefficents using an expansion from
//// the newton polynomial. Note: this algorithm is a bit tricky. If you don't follow it right
//// away, I'd suggest thinking about deriving it yourself directly and compare results. If you
//// still are having difficulty, it's okay to skim this section. The important things are that
//// (1) you know that polynomial interpolation is not magic, and (2) you understand where the
//// runtime complexity comes from in the algorithm.
////
//// A newton polynomial looks like:
////
//// ```
//// p(x) = k0 + k1*(x-x0) + k2*(x-x0)*(x-x1) + k3*(x-x0)*(x-x1)*(x-x2) + ...
//// ```
////
//// Our algorithm will maintain two polynomials (stored in monomial coefficent form):
////
//// 1. The `p(x)` polynomial constructed so far
//// 2. The `z(x)` polynomial of `(x - c)` factors so far
////
//// We can iteratively add a new interpolant point each step with the steps:
////
//// 1. Solve for `k` by evaluating `p(x)` and `z(x)`
//// 2. Update `p` as  `p = p + k * z`
//// 3. Update `z` as `z = z * (x - c)`
////
//// Since we'll maintain `p(x)` in monomial coefficent form, the final answer is
//// just `p(x)` whenever we finish adding new points.
////
//// This algorithm has runtime complexity `O(n^2)` where `n` is the number of evaluations
//// interpolated. This is because we're doing `O(n)` polynomial evaluations and updates
//// in the inner loop, and we're looping over `n` points.
fn poly_interpolate(evals: &[GF], xs: &[GF]) -> Vec<GF> {
    // Sanity check
    assert_eq!(evals.len(), xs.len());

    // Setup
    let n = xs.len();
    let mut p = vec![GF::new(0); n];
    let mut z = vec![GF::new(0); n];
    z[0] = GF::new(1);

    // Loop over each evaluation point
    for i in 0..n {
        // Unpack the point for this iteration
        let x0 = xs[i];
        let y0 = evals[i];

        // Step 1: Solve k = (y0 - p(x0)) / z(x0)
        let k = (y0 - poly_eval_point(&p, x0)) / poly_eval_point(&z, x0);

        // Step 2: Update the polynomial: p' = p + k * z
        for j in 0..n {
            p[j] = p[j] + k * z[j];
        }

        // Step 3: Update the "zero" polynomial: z' = z * (x - x0)
        //
        // Let's simplify a bit first by distributing:
        //   z' = z * x - z * x0
        //
        // Now notice that multiplying by `x` is a shift of the coefficents.
        // This leads to the form we'll use:
        //   z' = shift(z) + (-x0)*z
        //
        let x0_neg = x0.negate();
        let mut prev = GF::new(0);
        for j in (1..n).rev() { // update in reverse
            z[j] = z[j-1] + x0_neg * z[j];
        }
        z[0] = x0_neg * z[0];
    }

    // The final result is just `p`
    p
}

//// #### Reed-Solomon Encoding and Decoding
//// Okay, time for encoding and decoding functions. Similar to how we built fields, we'll be hardcoding the parameters for this
//// particular Reed-Solomon code. Feel free to change them and experiment!

const K: usize = 3;                     // Number of data elements (before encoding)
const N: usize = 5;                     // Number of encoded elements

//// We also need to specify the evaluation points. We need to have `N` of these.
//// For this variant of Reed-Solomon, it really doesn't matter what we pick. They just need to be
//// distinct elements of `GF(256)`. We choose: `[42, 222, 2, 8, 99]`.
const X: [GF; 5] = [
    GF::new(42), GF::new(222), GF::new(2), GF::new(8), GF::new(99)
];

//// Let's write the encode routine. This is just a trivial wrapper around polynomial evaluation!
fn reed_solomon_poly_encode(data: &[GF]) -> Vec<GF> {
    // Sanity check
    assert_eq!(data.len(), K);
    assert_eq!(X.len(), N);

    // Just treat the data as if it's polynomial coefficents and evaluate!
    poly_eval(data, &X)
}

//// Let's write the decode routine. We'll specify erasures using `Some(x)` and `None` in Rust. This is not
//// a common way to do it, but it illustrates the idea of erasures very well. Also, we have the possiblity
//// that decode could fail if too many erasures happen. This is indicated by returning None.
fn reed_solomon_poly_decode(encoded: &[Option<GF>]) -> Option<Vec<GF>> {
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
//// You always write tests for you code too... right? Riiiight!? I always do of course (**wink**)

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_eval() {
        assert_eq!(poly_eval_point(&[], GF::new(1)), GF::new(0));
        assert_eq!(poly_eval_point(&[GF::new(2)], GF::new(2)), GF::new(2));
        assert_eq!(poly_eval_point(&[GF::new(2), GF::new(2)], GF::new(3)), GF::new(4));
        assert_eq!(poly_eval_point(&[GF::new(46), GF::new(198), GF::new(0), GF::new(89)], GF::new(201)), GF::new(126));
    }

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

    #[test]
    fn test_encode_decode() {
        let data = vec![GF::new(100), GF::new(150), GF::new(200)];

        // encode
        let enc = reed_solomon_poly_encode(&data);
        assert_eq!(enc, vec![GF::new(160), GF::new(135), GF::new(94), GF::new(104), GF::new(194)]);
        let recv_all: Vec<_> = enc.iter().map(|x| Some(*x)).collect();

        // decode with no erasures
        let mut recv = recv_all.clone();
        assert_eq!(reed_solomon_poly_decode(&recv), Some(data.clone()));

        // erase element 0 and decode
        let mut recv = recv_all.clone();
        recv[0] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), Some(data.clone()));

        // erase elements 1,2 and decode
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), Some(data.clone()));

        // erase elements 2,4 and decode
        let mut recv = recv_all.clone();
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), Some(data.clone()));

        // erase elements 1,2,3 and decode will fail
        let mut recv = recv_all.clone();
        recv[1] = None;
        recv[2] = None;
        recv[3] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), None);

        // erase elements 0,2,4 and decode will fail
        let mut recv = recv_all.clone();
        recv[0] = None;
        recv[2] = None;
        recv[4] = None;
        assert_eq!(reed_solomon_poly_decode(&recv), None);
    }
}
