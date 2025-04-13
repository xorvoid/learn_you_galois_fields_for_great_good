use std::ops::{Add, Div, Mul, Sub};

//// #### An implementation of Computer Science Fields GF(2^k)
////
//// After implementing the general `GF(p^k)` Fields, why would we implement a special case?
//// Can't we just use that implemention, configured for `p = 2`?
////
//// Fair question. You could indeed use it. But, often it's a bit overkill and significantly
//// better implementations exist for `GF(2^k)`
////
//// Let us explain.
////
//// Recall that we used a vector representation for `GF(p^k)`. In computer memory, this was an
//// array of bytes: one byte per coefficient. What happens when we use `p = 2`?  Well, it becomes
//// a vector of bits (0 or 1). And, we have a very efficient way to store that vector in a computer:
////
//// *An Ordinary Binary Number*
////
//// Consider an 8-bit unsigned integer (`u8`), we can view this as a *vector of 8 bits*.
//// Similarly, a 16-bit unsigned integer (`u8`) can be viewed as a *vector of 16 bits*. And so on!
////
//// Nifty, huh?
////
//// This is why I call these fields the *Computer Science Fields*. These map extremely well onto
//// ordinary binary-state transistor computers. And, it turns out that common bitwise operations
//// are useful as well! We'll discuss these as we implement `GF(2^k)`
////

pub const K: usize = 3; // Must be less-or-equal 64, as we are using u64 for the bitvector
pub const Q: u64 = 11; // Default: x^3 + x + 1

//// We will store polynomials as a single u64 because all of our coefficients
//// are in {0,1}
////
//// For example x^3 + x^2 (+ 0x) + 1 can be represented as 0b1101
////
//// This is slightly different than the previous integer representation
//// conversion scheme we used for GF(P^K).
////
//// Exercize: or is it?

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GF2K(u64);

//// The basic utility functions are extremely straightforward now; we have
//// almost no work to do here

impl GF2K {
  pub fn new(val: u64) -> GF2K {
    // Sanity check!
    assert!((val as usize) < GF2K::number_of_elements());
    GF2K(val)
  }

  pub fn number_of_elements() -> usize {
    let k: u32 = K.try_into().unwrap(); // Abort if the number doesn't fit in 32-bits
    let p_k = (2 as u32).checked_pow(k).unwrap(); // Abort if the number doesn't fit in 32-bits
    p_k as usize
  }

  pub fn value(&self) -> u64 {
    self.0
  }
}

//// Addition in GF(2^k)
////
//// Our coefficients are all in GF(2). Let's revisit the addition table for
//// GF(2).
////
////  +  |  0    1
//// ----|----------
////  0  |  0    1
////  1  |  1    0
////
//// Addition is XOR of the two digits.
////
//// Therefore, polynomial addition, of two bitvectors is just the xor of the
//// two bitvectors
////
////          1x^2 + 1x + 1   =>    0b0111
//// + 1x^3 + 0x^2 + 1x + 1   => ^  0b1011
//// ----------------------      ---------
////   1x^3 + 1x^2 + 0x + 0  =>     0b1100

impl Add<GF2K> for GF2K {
  type Output = GF2K;
  fn add(self, rhs: GF2K) -> GF2K {
    GF2K::new(self.0 ^ rhs.0)
  }
}

//// Negation in GF(2)
////
//// Look carefully at the addition table for the GF(2) coefficients.
////
//// To negate, we want the value b, such that a+b=0
//// For 0, that is 0
//// For 1, that is 1
////
//// So negate(a) = a in GF(2)

impl GF2K {
  pub fn negate(self) -> GF2K {
    self
  }
}

//// Subtraction in GF(2)
////
//// Trivially, just add a + inverse(b) to subtract!
//// inverse(b)=b, so really a-b=a+b

impl Sub<GF2K> for GF2K {
  type Output = GF2K;
  fn sub(self, rhs: GF2K) -> GF2K {
    self + rhs.negate() // or just `self + rhs`
  }
}

//// Multiplication in GF(2^K)
////
//// This function performs polynomial multiplication in GF(2), where each bit in
//// the u64 represents a coefficient of a polynomial.
////
//// When we know the coefficients are all in GF(2), the convolution algorithm
//// simplifies significantly:
//// - Multiplication of bits is logical AND
//// - Addition of bits is logical XOR
////
//// For each position in the result polynomial, we compute the coefficient by:
//// 1. Finding all pairs of positions (j, i-j) from a and b that contribute to position i
//// 2. ANDing each pair of bits (multiplication in GF(2))
//// 3. XORing all products together (addition in GF(2))
////
//// Example in GF(2^3):
//// Let a = x^2 + 1 (0b101) and b = x^2 + x (0b110)
////
//// Computing some coefficients:
//// - Position 1: (a[0]&b[1]) ^ (a[1]&b[0]) = (1&1) ^ (0&0) = 1
//// - Position 2: (a[0]&b[2]) ^ (a[1]&b[1]) ^ (a[2]&b[0]) = (1&1) ^ (0&1) ^ (1&0) = 1
////
//// Final result: x^4 + x^3 + x^2 + x (0b11110)
////
//// The result has degree up to 2*K-2, which will later be reduced modulo Q using poly_mod.

fn poly_mul(a: u64, b: u64) -> u64 {
  // Result goes into c
  let mut c: u64 = 0;

  // For each possible result position (in a full multiplication)
  for i in 0..(2 * K - 1) {
    // Calculate the i-th coefficient of the result
    let mut coef = 0;
    for j in 0..K {
      // Skip if the j-th coefficient from a or (i-j)-th coefficient from b is out of bounds
      if i < j || (i - j) >= K {
        continue;
      }

      // Get the j-th bit of a and (i-j)-th bit of b
      let bit_a = (a >> j) & 1;
      let bit_b = (b >> (i - j)) & 1;
      coef ^= bit_a & bit_b; // coef += a * b
    }

    c |= coef << i; // set the bit high, if the sum is 1
  }

  c
}

//// Modulus in GF(2^K)
////
//// Use two assumptions:
//// 1. The polynomial `a` is a result of some multiplication; it is of order 2*K - 2
//// 2. q is some irreducible polynomial of order K
////
//// So we (again) only need to eliminate terms starting at high order and
//// working our way down to K
////
//// Following the same algorithm we used in gf_p_k, we can just apply the
//// shifts. Notice that, because the coefficients can only be in {0,1}, we'll
//// never need to perform any scaling.
////
//// Let's work a short example in GF(2^3):
////
//// Using q = x^3 + x + 1 (0b1011)
//// Let a = x^4 + x^2 + 1 (0b10101)
////
//// Start at i=K*2 - 1 = 3*2 - 2 = 6 - 2 = 4:
////   a[4] = 1
////
//// So we need to shift q (of order K=3) to be order 4
////   shift = 4-3 = 1
////
//// Which, in polynomial representation is:
////   (x^3 + x + 1) * x = x^4 + x^2 + x
////
//// Notice that we can do this simply by shifting the entire integer q with a
//// bitshift.
////   q << shift
////   == 0b1011 << 1
////   == 0b10110      (same thing as x^4 + x^2 + x)
////
////
//// We can now "add" together the coefficients all at once using an XOR:
////   a = a ^ q_shifted
////
//// The rest of the algorithm works the same way as before; notice that we've
//// effectively vectorized over the addition again by using the xor!

fn poly_mod(mut a: u64, q: u64) -> u64 {
  // We'll iterate from high-terms to low-terms, eliminating each:
  //   2K-2, 2K-3, ..., K+1, K
  // This will leave exactly K coefficients
  for i in (K..=(2 * K - 2)).rev() {
    let a_i = a & (1u64 << i); // extract `i`th coefficient

    // If the i-th bit is already 0, nothing to do
    if a_i == 0 {
      continue;
    };

    let shift = i - K; // coefficient shift
    let q_shifted = q << shift; // shifting with coefs in GF(2) is a bit shift

    // add (xor) q_shifted into a (elimination the leading term)
    a = a ^ q_shifted;

    // the i'th term should be zero now
    assert_eq!(a & (1u64 << i), 0);
  }

  // All bits at position K and higher should be zero now
  a
}

//// Just like before, multiplication is actually two steps:
//// 1. Multiply (produces a higher order polynomial)
//// 2. Take mod to get back into the field

impl Mul<GF2K> for GF2K {
  type Output = GF2K;
  fn mul(self, rhs: GF2K) -> GF2K {
    let res = poly_mod(poly_mul(self.0, rhs.0), Q);
    GF2K(res)
  }
}

//// Just like before, we will just brute force these
////
//// Is there a better way to do this in GF(2^K)?

impl GF2K {
  pub fn invert(self) -> Result<GF2K, String> {
    // Important: Zero has no inverse, it's invalid
    if self == GF2K::new(0) {
      return Err("Zero has no inverse".to_string());
    }
    // Scan the numbers {1, 2, ..., P-1} until we find the inverse
    for x in 1..GF2K::number_of_elements() {
      let candidate = GF2K::new(x as u64);
      if self * candidate == GF2K::new(1) {
        return Ok(candidate); // Found!
      }
    }
    unreachable!("Every non-zero number has an inverse");
  }
}

//// And division can be exactly the same also.
impl Div<GF2K> for GF2K {
  type Output = Result<GF2K, String>;
  fn div(self, rhs: Self) -> Result<GF2K, String> {
    // Important: Cannot divide by zero
    if rhs == GF2K::new(0) {
      return Err("Cannot divide by zero".to_string());
    }
    Ok(self * rhs.invert().unwrap())
  }
}

//// #### Some final things
////
//// Just as before, we need to teach Rust a few extra tricks.
//// These are all quite similar to `GF(p^k)` if you'd rather skim them.
////
//// Printing out numbers:

impl std::fmt::Display for GF2K {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    write!(f, "{}", self.value())
  }
}

//// And converting strings into our Field's numbers:
impl std::str::FromStr for GF2K {
  type Err = String;
  fn from_str(s: &str) -> Result<GF2K, String> {
    let num: u64 = s.parse().map_err(|_| format!("Not an 64-bit integer"))?;
    // Return an error if the number is too big for the field
    let limit = GF2K::number_of_elements() as u64;
    if num >= limit {
      return Err(format!(
        "Number too large, got {}, but limit is {}",
        num, limit
      ));
    }
    Ok(GF2K::new(num))
  }
}

//// And telling Rust that we built a Field type:
impl crate::field::Field for GF2K {
  fn number_of_elements() -> usize {
    GF2K::number_of_elements()
  }
}

//// #### Testing Time
////
//// Note that these tests assume `GF(2^3)`. If you change the field, they are not expected to pass.
////
//// The year is also 2025, so obviously these were written by AI (then fixed,
//// because claude wrote failing tests).

#[cfg(test)]
mod tests {
  use super::*;

  // Ensure we're using a specific field size for testing
  #[test]
  fn test_field_size() {
    // For K=3, we should have 2^3 = 8 elements
    assert_eq!(GF2K::number_of_elements(), 8);
  }

  // TEST: Verify basic field element construction
  #[test]
  fn test_element_construction() {
    for i in 0..8 {
      let element = GF2K::new(i);
      assert_eq!(element.value(), i);
    }
  }

  // TEST: We shouldn't be able to construct numbers out of the range
  #[should_panic]
  #[test]
  fn test_invalid_numbers() {
    GF2K::new(8); // With K=3, the largest valid value is 7
  }

  // TEST: Addition (which is XOR in GF(2^K))
  #[test]
  fn test_add() {
    assert_eq!(GF2K::new(0) + GF2K::new(0), GF2K::new(0));
    assert_eq!(GF2K::new(0) + GF2K::new(1), GF2K::new(1));
    assert_eq!(GF2K::new(1) + GF2K::new(1), GF2K::new(0));
    assert_eq!(GF2K::new(1) + GF2K::new(2), GF2K::new(3));
    assert_eq!(GF2K::new(3) + GF2K::new(4), GF2K::new(7));
    assert_eq!(GF2K::new(5) + GF2K::new(7), GF2K::new(2));
    assert_eq!(GF2K::new(6) + GF2K::new(3), GF2K::new(5));
  }

  // TEST: Negation (which is identity in GF(2^K))
  #[test]
  fn test_negation() {
    // In GF(2^K), negation is the identity function
    for i in 0..8 {
      assert_eq!(GF2K::new(i).negate(), GF2K::new(i));
    }
  }

  // TEST: Subtraction (which is identical to addition in GF(2^K))
  #[test]
  fn test_sub() {
    // In GF(2^K), addition and subtraction are the same operation
    // GF(2^3) | 0 - 0 = 0
    assert_eq!(GF2K::new(0) - GF2K::new(0), GF2K::new(0));

    // GF(2^3) | 0 - 1 = 1
    assert_eq!(GF2K::new(0) - GF2K::new(1), GF2K::new(1));

    // GF(2^3) | 1 - 1 = 0
    assert_eq!(GF2K::new(1) - GF2K::new(1), GF2K::new(0));

    // GF(2^3) | 3 - 2 = 1
    assert_eq!(GF2K::new(3) - GF2K::new(2), GF2K::new(1));

    // GF(2^3) | 7 - 5 = 2
    assert_eq!(GF2K::new(7) - GF2K::new(5), GF2K::new(2));
  }

  #[test]
  fn test_poly_operations() {
    assert_eq!(poly_mul(2, 2), 4);
    assert_eq!(poly_mul(3, 3), 5);
    assert_eq!(poly_mul(5, 3), 15);

    assert_eq!(poly_mod(8, 11), 3);
    assert_eq!(poly_mod(12, 11), 7);
    assert_eq!(poly_mod(15, 11), 4);
  }

  // TEST: Multiplication using the irreducible polynomial Q = x^3 + x + 1 (0b1011)
  #[test]
  fn test_mul() {
    assert_eq!(GF2K::new(0) * GF2K::new(0), GF2K::new(0));
    assert_eq!(GF2K::new(0) * GF2K::new(1), GF2K::new(0));
    assert_eq!(GF2K::new(1) * GF2K::new(1), GF2K::new(1));
    assert_eq!(GF2K::new(2) * GF2K::new(2), GF2K::new(4));
    assert_eq!(GF2K::new(3) * GF2K::new(3), GF2K::new(5));
    assert_eq!(GF2K::new(4) * GF2K::new(4), GF2K::new(6));
    assert_eq!(GF2K::new(5) * GF2K::new(6), GF2K::new(3));
    assert_eq!(GF2K::new(7) * GF2K::new(7), GF2K::new(3));
  }

  #[test]
  fn test_invert() {
    // 0 has no inverse
    assert!(matches!(GF2K::new(0).invert(), Err(_)));

    // 1 is its own inverse
    assert_eq!(GF2K::new(1).invert().unwrap(), GF2K::new(1));

    // Precomputed inverses for GF(2^3) with Q = 11
    let inverse_table = [
      0, // 0 has no inverse
      1, // inv(1) = 1
      5, // inv(2) = 5
      6, // inv(3) = 6
      7, // inv(4) = 7
      2, // inv(5) = 2
      3, // inv(6) = 3
      4, // inv(7) = 4
    ];

    for i in 1..8 {
      assert_eq!(
        GF2K::new(i).invert().unwrap(),
        GF2K::new(inverse_table[i as usize])
      );

      // Verify that x * x^(-1) = 1
      assert_eq!(
        GF2K::new(i) * GF2K::new(inverse_table[i as usize]),
        GF2K::new(1)
      );
    }
  }

  #[test]
  fn test_div() {
    // Division by zero is an error
    assert!(matches!(GF2K::new(1) / GF2K::new(0), Err(_)));

    // 0 divided by anything non-zero is 0
    for i in 1..8 {
      assert_eq!(GF2K::new(0) / GF2K::new(i), Ok(GF2K::new(0)));
    }

    // Division by 1 is identity
    for i in 0..8 {
      assert_eq!(GF2K::new(i) / GF2K::new(1), Ok(GF2K::new(i)));
    }

    // For all non-zero a and b, a/b = a * inv(b)
    for a in 1..8 {
      for b in 1..8 {
        let expected = GF2K::new(a) * GF2K::new(b).invert().unwrap();
        assert_eq!(GF2K::new(a) / GF2K::new(b), Ok(expected));
      }
    }

    // Some specific test cases
    assert_eq!(GF2K::new(2) / GF2K::new(3), Ok(GF2K::new(7)));
    assert_eq!(GF2K::new(5) / GF2K::new(6), Ok(GF2K::new(4)));
    assert_eq!(GF2K::new(7) / GF2K::new(4), Ok(GF2K::new(3)));
  }

  // TEST: Distributivity law (a+b)*c = a*c + b*c
  #[test]
  fn test_distributivity() {
    for a in 0..8 {
      for b in 0..8 {
        for c in 0..8 {
          let left = (GF2K::new(a) + GF2K::new(b)) * GF2K::new(c);
          let right = GF2K::new(a) * GF2K::new(c) + GF2K::new(b) * GF2K::new(c);
          assert_eq!(left, right);
        }
      }
    }
  }

  // TEST: Associativity of multiplication (a*b)*c = a*(b*c)
  #[test]
  fn test_associativity() {
    for a in 0..8 {
      for b in 0..8 {
        for c in 0..8 {
          let left = (GF2K::new(a) * GF2K::new(b)) * GF2K::new(c);
          let right = GF2K::new(a) * (GF2K::new(b) * GF2K::new(c));
          assert_eq!(left, right);
        }
      }
    }
  }
}
