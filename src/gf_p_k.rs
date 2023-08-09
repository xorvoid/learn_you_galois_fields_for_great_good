//// #### An implementation of GF(p^k) where p is a prime number and k is any integer >= 1
////
//// We will implement addition, subtraction, multiplication, and division using
//// operator-overloading so we can use the normal operators: +, -, *, /
use std::ops::{Add, Mul, Sub, Div};

//// Similar to our implementation of `GF(p)` fields, we'll define the constant paramaters
//// here instead of using advanced Rust features.
////
//// In particular, we need:
////
////   - `P`: A prime value for the coefficients field `GF(p)`
////   - `K`: A polynomial order limit parameter
////   - `Q`: An irredicible polynomial for the modulus operation
////
//// We use some defaults below for a `GF(3^2)` Field.
////
//// Feel free change the parameters to get a different field. But, please do be careful to configure
//// correctly. Most notably: `P` must be prime, and `Q` must be irreducible. We will demostrate tooling
//// later to help find irredicuble polynomials.
pub const P: usize = 3;  // Default: GF(3)
pub const K: usize = 2;  // Default: All polynomials below x^2
pub const Q: usize = 10; // Default: x^2 + 1

//// Our representation for a number in `GF(p^k)` will use the vector representation with unsigned 8-bit integers (`u8`).
//// Each element in the vector is an element of `GF(p)`..
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GF([u8; K]); // An array of K elements of type u8


//// The number of elements we have in `GF(p^k)` is clearly `p^k`, but we need a little care to
//// avoid silent integer overflow.
impl GF {
  pub fn number_of_elements() -> usize {
    let p: u32 = P.try_into().unwrap(); // Abort if the number doesn't fit in 32-bits
    let k: u32 = K.try_into().unwrap(); // Abort if the number doesn't fit in 32-bits
    let p_k = p.checked_pow(k).unwrap(); // Abort if the number doesn't fit in 32-bits
    p_k as usize
  }
}

//// Let's construct a new number in `GF(p^k)`. Our function will accept a number in the "decimal representation" and will convert it into the internal vector representation.
////
//// It's easy to mix up ordering and "endianness", so we'll be explicit about how we're using  vector representation:
////
//// ```
//// index 0 => 1
//// index 1 => x
//// index 2 => x^2
//// index 3 => x^3
//// index 4 => x^4
//// ... etc ...
//// ```
////
//// Okay, let's put it all together

impl GF {
  pub fn new(mut val: u64) -> GF {
    // Sanity check!
    assert!((val as usize) < GF::number_of_elements());

    // Decompose the decimal represenation into a vector
    let p = P as u64;
    let mut vec = [0; K];
    for i in 0..K {
      vec[i] = (val % p) as u8;
      val /= p;
    }

    GF(vec)
  }
}

//// And we'll want to convert back to decimal representation also. This is effectively a
//// dot-product with a Radix-P Basis: `[1, P, P^2, P^3, ...]`
impl GF {
  pub fn value(&self) -> u64 {
    // Convert vector representation back into decimal representation
    let p = P as u64;
    let mut val = 0;
    for i in (0..K).rev() {
      val = val * p + self.0[i] as u64;
    }
    val
  }
}

//// #### Coefficient Field Implementation
////
//// We will quickly re-implement parts of `GF(p)` here for the coefficients.
////
//// We will need these internally to build the full Field.
//// We could inline them, but doing so is distracting when reading the polynomial arithmetic
//// algorithms.
////
//// Please refer to the `GF(p)` implementation if these aren't clear:
fn coeff_add(a: u8, b: u8) -> u8 {
  (((a as u16) + (b as u16)) % (P as u16)) as u8
}
fn coeff_neg(a: u8) -> u8 {
  ((P as u8) - a) % (P as u8)
}
fn coeff_mul(a: u8, b: u8) -> u8 {
  (((a as u16) * (b as u16)) % (P as u16)) as u8
}

//// #### Addition and Subtraction
////
//// Addition is very simple. It's just a simple pointwise-addition of coefficients (vector addition).
impl Add<GF> for GF {
  type Output = GF;
  fn add(self, rhs: GF) -> GF {
    let mut vec = [0; K];
    for i in 0..K {
      vec[i] = coeff_add(self.0[i], rhs.0[i]);
    }
    GF(vec)
  }
}

//// Negation (or additive-inverse) is fairly straight-forward also. Just vectorized.
impl GF {
  pub fn negate(self) -> GF {
    let mut vec = [0; K];
    for i in 0..K {
      vec[i] = coeff_neg(self.0[i]);
    }
    GF(vec)
  }
}

//// And subtraction: trivial.
////
//// Getting bored? Be careful what you wish for...
impl Sub<GF> for GF {
  type Output = GF;
  fn sub(self, rhs: GF) -> GF {
    self + rhs.negate()
  }
}

//// #### Multiplication and Division
////
//// Multiplication is a completely different beast than addition. You were warned ;-).
////
//// Recall that multiplication is a convolution of the coefficients. We can compute a convolution
//// by effecively flipping around one of the vectors and computing sliding dot-products.
////
//// Let's do some examples:
////
//// ```text
//// Consider:
////   A = [a, b, c, d]
////   B = [e, f, g, h]
////
//// We want to compute:
////   C = convolve(A, B)
////     = [c_0, c_1, c_2, c_3, c_4, c_5, c_6]
////
//// We can compute each of these coefficients (`c_i`) as follows:
////
////   c_0
////   ----------------------------
////   [d, c, b, a]
////             *
////            [e, f, g, h]
////
////   = a*h
////
////   c_1
////   ----------------------------
////   [d, c, b, a]
////          *  *
////         [e, f, g, h]
////
////   = b*e + a*f
////
////   c_2
////   ----------------------------
////   [d, c, b, a]
////       *  *  *
////      [e, f, g, h]
////
////   = c*e + b*f + a*g
////
////
////   ----------------------------
////   ... SNIP ...
////   ----------------------------
////
////
////   c_5
////   ----------------------------
////         [d, c, b, a]
////          *  *
////   [e, f, g, h]
////
////   = d*g + c*h
////
////   c_6
////   ----------------------------
////            [d, c, b, a]
////             *
////   [e, f, g, h]
////
////   = d*h
////
//// ```
////
//// As you can see, computing a convolution on two arrays is a fairly simple shift-and-multiply
//// operation. It can be a LOT of multiply operations, in fact `O(n^2)`, but it's otherwise fairly
//// formulaic.
////
//// A few things to note:
////
//// 1. Instead of actually reversing a vector, we'll just use backwards array indexing
//// 2. Observe that a convolution of two k-element vectors produces a (2k-1)-element result!
////
//// Here's some psuedocode for the algorithm:
////
//// ```text
//// for i in 0..(2k-1) {
////   c[i] = 0
////   for j in 0..k {
////     c[i] += a[i - j] * b[j]     (where out-of-bounds indicies implicitly load 0)
////   }
//// }
//// ```
////
//// And here's the real implementation that deals with all the real-world practical messiness. Compare it to the psuedocode if it's confusing. The same algorithmic structure is present.

fn poly_mul(a: &[u8; K], b: &[u8; K]) -> [u8; 2*K-1] {
  // A convolution implementation over the field GF(p)
  let mut c = [0; 2*K-1];
  for i in 0..(2*K-1) {
    // Each coefficient is the sum of many sub-terms
    for j in 0..K {
      // Ignore terms with out-of-bounds indicies (they are implicitly zero)
      if i < j || i - j >= K {
        continue;
      }
      // Add this term to the result (NOTE: poly_a is reversed!)
      c[i] = coeff_add(c[i], coeff_mul(a[i - j], b[j]));
    }
  }
  c
}

//// Great fun, huh!
////
//// Now we need to implement polynomial modulus.
////
//// This one is also tricky because long-division is tedious. Fortunately, it reduces to
//// an algorithm just as cleanly as convolution.
////
//// A key observation is that we don't care about the result of the division. We can discard that.
//// We only care about what remains after removing some number of `Q`s from `A`.
////
//// Consider:

//// ```text
//// A = 2x^5 + x^4 + 1
//// Q = x^3 + 1
//// ```
////
//// And suppose we want to compute:
////
//// ```text
//// A % Q
//// ```
////

//// We can do this by starting at the largest term (`2x^5`) and eliminating it.
//// After, we'll have just a 4th order polynomial. Then, we can repeat until we have a remainder.
////
//// For this, we'll "shift" and "scale" up `Q` and subtract it from `A`:
////
//// ```text
//// A' = A - 2x^2 * Q
////    = (2x^5 + x^4 + 1) - 2x^2 * (x^3 + 1)        [expand]
////    = (2x^5 + x^4 + 1) + (-2x^2 * (x^3 + 1))     [re-arrange]
////    = (2x^5 + x^4 + 1) + (x^2 * (x^3 + 1))       [negate in GF(3)]
////    = (2x^5 + x^4 + 1) + (x^5 + x^2)             [multiply: "shift" and "scale"]
////    = (2x^5 + x^4 + 1) + (x^5 + x^2)             [add in GF(3)]
////    = x^4 + x^2 + 1                              [result: polynomial of order 4!]
//// ```
////
//// Next, we can eliminate the `x^4` term by computing:
////
//// ```text
//// A'' = A' - x * Q
////     = (x^4 + x^2 + 1) - x * (x^3 + 1)           [expand]
//// ```
////
//// <i><u>Exercise:</u></i> Compute this to confirm your understanding.
////
//// In order to do this, we need to multiply `Q`. But this might seem like some kind of chicken-and-egg problem?
////
//// Well actually, we're doing two simpler operations:
////
//// - `shift`: Multiplying by a monomial in `(1, x, x^2, x^3, ...etc...)` is a simple shift of the coefficent array
//// - `scale`: Multiplying by any scalar is a simple pointwise-coefficent multiply.
////
//// So we don't need normal polynomial multiplication to do this.
////
//// Okay, let's code!

fn poly_mod(mut a: [u8; 2*K-1], q: &[u8; K+1]) -> [u8; K] {
  // We'll iterate from high-terms to low-terms, eliminating each:
  //   2k-2, 2k-3, ..., k+1, k
  // This will leave exactly k coefficents
  for i in (K..(2*K-1)).rev() {

    // Determine "our shift" and "scale" for `Q`
    let shift = i - K; // coefficent shift
    let scale = coeff_neg(a[i]); // coefficent scale

    // Apply to each coefficient, shifted
    for j in shift..=i { // inclusive: [shift, i]
      // Shift and Scale up one coefficient of Q
      let val = coeff_mul(scale, q[j-shift]);
      // Now "subtract" it (using addition)
      a[j] = coeff_add(a[j], val);
    }
    // The largest term (ith) should be zero now!
    assert_eq!(a[i], 0);
  }
  // Only the lowest k-terms are possibly non-zero, now
  a[..K].try_into().unwrap()
}

//// Lastly, we need a little routine to convert our irredicible polynomial into the vector
//// notation used by `poly_mod()`
fn poly_q() -> [u8; K+1] {
  // Decompose the decimal represenation into a vector
  let p = P as u64;
  let mut val = Q as u64;
  let mut vec = [0; K+1];
  for i in 0..K+1 {
    vec[i] = (val % p) as u8;
    val /= p;
  }
  vec
}

//// With these polynomial operations, our final multiply is simple to complete.
//// It's just the defintion that we gave in the previous theory
//// section: `poly_mod(poly_mul(A, B), Q)`
impl Mul<GF> for GF {
  type Output = GF;
  fn mul(self, rhs: GF) -> GF {
    let vec = poly_mod(poly_mul(&self.0, &rhs.0), &poly_q());
    GF(vec)
  }
}

//// And now I suspect the dear reader is a tad exhausted from our journey through the *Polynomial Math Wilderness*.
////
//// Well it's your lucky day! We can adapt the rest from `GF(p)` fields by reusing our simple
//// naive brute-force approach to inverses.
impl GF {
  pub fn invert(self) -> Result<GF, String> {
    // Important: Zero has no inverse, it's invalid
    if self == GF::new(0) {
      return Err("Zero has no inverse".to_string());
    }
    // Scan the numbers {1, 2, ..., P-1} until we find the inverse
    for x in 1..GF::number_of_elements() {
      let candidate = GF::new(x as u64);
      if self * candidate == GF::new(1) {
        return Ok(candidate); // Found!
      }
    }
    unreachable!("Every non-zero number has an inverse");
  }
}

//// And division can be exactly the same also.
impl Div<GF> for GF {
  type Output = Result<GF, String>;
  fn div(self, rhs: Self) -> Result<GF, String> {
    // Important: Cannot divide by zero
    if rhs == GF::new(0) {
      return Err("Cannot divide by zero".to_string());
    }
    Ok(self * rhs.invert().unwrap())
  }
}

//// #### Some final things
////
//// Just as before, we need to teach Rust a few extra tricks.
//// These are all quite similair to `GF(p)` if you'd rather skim them.
////
//// Printing out numbers:
impl std::fmt::Display for GF {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    write!(f, "{}", self.value())
  }
}

//// And converting strings into our Field's numbers:
impl std::str::FromStr for GF {
  type Err = String;
  fn from_str(s: &str) -> Result<GF, String> {
    let num: u64 = s.parse().map_err(|_| format!("Not an 64-bit integer"))?;
    // Return an error if the number is too big for the field
    let limit = GF::number_of_elements() as u64;
    if num >= limit {
      return Err(format!("Number too large, got {}, but limit is {}", num, limit));
    }
    Ok(GF::new(num))
  }
}

//// And telling Rust that we built a Field type:
impl crate::field::Field for GF {
  fn number_of_elements() -> usize {
    GF::number_of_elements()
  }
}

//// #### Testing Time
////
//// Note that these tests assume `GF(3^2)`. If you change the field, they are not expected to pass.
#[cfg(test)]
mod tests {
  use super::*;

  // TEST: Verify that decimal representations are converted to vector representation correctly
  #[test]
  fn test_convert() {
    assert_eq!(GF::new(0), GF([0, 0]));
    assert_eq!(GF::new(1), GF([1, 0]));
    assert_eq!(GF::new(2), GF([2, 0]));
    assert_eq!(GF::new(3), GF([0, 1]));
    assert_eq!(GF::new(4), GF([1, 1]));
    assert_eq!(GF::new(5), GF([2, 1]));
    assert_eq!(GF::new(6), GF([0, 2]));
    assert_eq!(GF::new(7), GF([1, 2]));
    assert_eq!(GF::new(8), GF([2, 2]));
  }

  // TEST: Verify conversion to vector representations and back to decimal representation is an identity
  #[test]
  fn test_convert_identity() {
    assert_eq!(GF::new(0).value(), 0);
    assert_eq!(GF::new(1).value(), 1);
    assert_eq!(GF::new(2).value(), 2);
    assert_eq!(GF::new(3).value(), 3);
    assert_eq!(GF::new(4).value(), 4);
    assert_eq!(GF::new(5).value(), 5);
    assert_eq!(GF::new(6).value(), 6);
    assert_eq!(GF::new(7).value(), 7);
    assert_eq!(GF::new(8).value(), 8);
  }

  // TEST: We shouldn't be able to construct numbers out of the range
  #[should_panic]
  #[test]
  fn test_invalid_numbers() {
    GF::new(9);
  }

  // TEST: Addition
  #[test]
  fn test_add() {
    // GF(3^2) | 0 + 1 = 1
    assert_eq!(GF::new(0) + GF::new(1), GF::new(1));
    // GF(3^2) | 1 + 1 = 2
    assert_eq!(GF::new(1) + GF::new(1), GF::new(2));
    // GF(3^2) | 1 + 2 = 0
    assert_eq!(GF::new(1) + GF::new(2), GF::new(0));
    // GF(3^2) | 2 + 2 = 1
    assert_eq!(GF::new(2) + GF::new(2), GF::new(1));
    // GF(3^2) | 1 + 3 = 4
    assert_eq!(GF::new(1) + GF::new(3), GF::new(4));
    // GF(3^2) | 2 + 3 = 5
    assert_eq!(GF::new(2) + GF::new(3), GF::new(5));
    // GF(3^2) | 3 + 3 = 6
    assert_eq!(GF::new(3) + GF::new(3), GF::new(6));
    // GF(3^2) | 2 + 4 = 3
    assert_eq!(GF::new(2) + GF::new(4), GF::new(3));
    // GF(3^2) | 4 + 4 = 8
    assert_eq!(GF::new(4) + GF::new(4), GF::new(8));
    // GF(3^2) | 2 + 7 = 6
    assert_eq!(GF::new(2) + GF::new(7), GF::new(6));
    // GF(3^2) | 6 + 7 = 4
    assert_eq!(GF::new(6) + GF::new(7), GF::new(4));
    // GF(3^2) | 8 + 8 = 4
    assert_eq!(GF::new(8) + GF::new(8), GF::new(4));
  }

  // TEST: Subtraction
  #[test]
  fn test_sub() {
    // GF(3^2) | 1 - 0 = 1
    assert_eq!(GF::new(1) - GF::new(0), GF::new(1));
    // GF(3^2) | 0 - 1 = 2
    assert_eq!(GF::new(0) - GF::new(1), GF::new(2));
    // GF(3^2) | 1 - 1 = 0
    assert_eq!(GF::new(1) - GF::new(1), GF::new(0));
    // GF(3^2) | 1 - 2 = 2
    assert_eq!(GF::new(1) - GF::new(2), GF::new(2));
    // GF(3^2) | 2 - 2 = 0
    assert_eq!(GF::new(2) - GF::new(2), GF::new(0));
    // GF(3^2) | 1 - 3 = 7
    assert_eq!(GF::new(1) - GF::new(3), GF::new(7));
    // GF(3^2) | 2 - 3 = 8
    assert_eq!(GF::new(2) - GF::new(3), GF::new(8));
    // GF(3^2) | 3 - 3 = 0
    assert_eq!(GF::new(3) - GF::new(3), GF::new(0));
    // GF(3^2) | 2 - 4 = 7
    assert_eq!(GF::new(2) - GF::new(4), GF::new(7));
    // GF(3^2) | 4 - 4 = 0
    assert_eq!(GF::new(4) - GF::new(4), GF::new(0));
    // GF(3^2) | 2 - 7 = 4
    assert_eq!(GF::new(2) - GF::new(7), GF::new(4));
    // GF(3^2) | 6 - 7 = 2
    assert_eq!(GF::new(6) - GF::new(7), GF::new(2));
    // GF(3^2) | 8 - 8 = 0
    assert_eq!(GF::new(8) - GF::new(8), GF::new(0));
  }

  // TEST: Multiplication
  #[test]
  fn test_mul() {
    // GF(3^2) | 0 * 1 = 0
    assert_eq!(GF::new(0) * GF::new(1), GF::new(0));
    // GF(3^2) | 0 * 2 = 0
    assert_eq!(GF::new(0) * GF::new(2), GF::new(0));
    // GF(3^2) | 1 * 1 = 1
    assert_eq!(GF::new(1) * GF::new(1), GF::new(1));
    // GF(3^2) | 1 * 2 = 2
    assert_eq!(GF::new(1) * GF::new(2), GF::new(2));
    // GF(3^2) | 2 * 2 = 1
    assert_eq!(GF::new(2) * GF::new(2), GF::new(1));
    // GF(3^2) | 2 * 3 = 6
    assert_eq!(GF::new(2) * GF::new(3), GF::new(6));
    // GF(3^2) | 3 * 3 = 2
    assert_eq!(GF::new(3) * GF::new(3), GF::new(2));
    // GF(3^2) | 2 * 4 = 8
    assert_eq!(GF::new(2) * GF::new(4), GF::new(8));
    // GF(3^2) | 3 * 4 = 5
    assert_eq!(GF::new(3) * GF::new(4), GF::new(5));
    // GF(3^2) | 4 * 4 = 6
    assert_eq!(GF::new(4) * GF::new(4), GF::new(6));
    // GF(3^2) | 2 * 7 = 5
    assert_eq!(GF::new(2) * GF::new(7), GF::new(5));
    // GF(3^2) | 6 * 7 = 8
    assert_eq!(GF::new(6) * GF::new(7), GF::new(8));
    // GF(3^2) | 8 * 8 = 6
    assert_eq!(GF::new(8) * GF::new(8), GF::new(6));
  }

  // TEST: Division
  #[test]
  fn test_div() {
    // GF(3^2) | 0 / 1 = 0
    assert_eq!(GF::new(0) / GF::new(1), Ok(GF::new(0)));
    // GF(3^2) | 0 / 2 = 0
    assert_eq!(GF::new(0) / GF::new(2), Ok(GF::new(0)));
    // GF(3^2) | 1 / 0 = ERROR
    assert!(matches!(GF::new(1) / GF::new(0), Err(_)));
    // GF(3^2) | 2 / 0 = ERROR
    assert!(matches!(GF::new(2) / GF::new(0), Err(_)));
    // GF(3^2) | 1 / 1 = 1
    assert_eq!(GF::new(1) / GF::new(1), Ok(GF::new(1)));
    // GF(3^2) | 1 / 2 = 2
    assert_eq!(GF::new(1) / GF::new(2), Ok(GF::new(2)));
    // GF(3^2) | 2 / 2 = 1
    assert_eq!(GF::new(2) / GF::new(2), Ok(GF::new(1)));
    // GF(3^2) | 2 / 3 = 3
    assert_eq!(GF::new(2) / GF::new(3), Ok(GF::new(3)));
    // GF(3^2) | 3 / 3 = 1
    assert_eq!(GF::new(3) / GF::new(3), Ok(GF::new(1)));
    // GF(3^2) | 1 / 4 = 5
    assert_eq!(GF::new(1) / GF::new(4), Ok(GF::new(5)));
    // GF(3^2) | 2 / 4 = 7
    assert_eq!(GF::new(2) / GF::new(4), Ok(GF::new(7)));
    // GF(3^2) | 3 / 4 = 8
    assert_eq!(GF::new(3) / GF::new(4), Ok(GF::new(8)));
    // GF(3^2) | 4 / 4 = 1
    assert_eq!(GF::new(4) / GF::new(4), Ok(GF::new(1)));
    // GF(3^2) | 1 / 7 = 5
    assert_eq!(GF::new(1) / GF::new(7), Ok(GF::new(8)));
    // GF(3^2) | 2 / 7 = 4
    assert_eq!(GF::new(2) / GF::new(7), Ok(GF::new(4)));
    // GF(3^2) | 6 / 7 = 5
    assert_eq!(GF::new(6) / GF::new(7), Ok(GF::new(5)));
    // GF(3^2) | 1 / 8 = 7
    assert_eq!(GF::new(1) / GF::new(8), Ok(GF::new(7)));
  }
}
