//// #### Implementing Binary Fields GF(2^k)
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
//// Similarly, a 16-bit unsigned integer (`u16`) can be viewed as a *vector of 16 bits*. And so on!
////
//// Let's also review the coefficent field `GF(2)`. Here are the addition and multiplication tables.
//// You may wish to review [section 2](galois_fields_for_great_good_02.html) and [section 3](galois_fields_for_great_good_03.html):
////
//// Addition:
////
////  **+**  |  **0**  | **1**
//// --------|---------|--------
////  **0**  |    0    |   1
////  **1**  |    1    |   0
////
//// Multiplication:
////
////  *****  |  **0**  | **1**
//// --------|---------|--------
////  **0**  |    0    |   0
////  **1**  |    0    |   1
////
//// <i><u>Exercise:</u></i> These correspond to two very well-known common bitwise operations, which are they?
////
//// That's right! Addition is a bitwise XOR and Multiplication is a bitwise AND!
////
//// Nifty, huh?
////
//// All of these properties make `GF(2^k)` fields incredibly useful for computer science applications as
//// they map very well onto ordinary binary-state transistor computers.
////
//// #### Enough chatter, let's go code!
////
//// As before, we will implement addition, subtraction, multiplication, and division using
//// operator-overloading so we can use the normal operators: +, -, *, /
use std::ops::{Add, Div, Mul, Sub};

//// Similar to our other implementations, we'll define the constant parameters
//// here instead of using advanced Rust features.
////
//// In particular, we need:
////
////   - `K`: A polynomial order limit parameter
////   - `Q`: An irreducible polynomial for the modulus operation
////
//// We use some defaults below for a `GF(2^3)` Field.
////
//// Feel free change the parameters to get a different field. But, please do be careful to configure
//// correctly. Most notably: `Q` must be irreducible.

pub const K: usize = 3; // Must be less-or-equal 64, as we are using u64 for the bitvector
pub const Q: u64 = 11; // Default: x^3 + x + 1

//// We will store polynomials as a single u64 because all of our polynomial coefficients
//// are in {0,1}
////
//// For example `x^3 + x^2 (+ 0x) + 1` can be represented as `0b1101` or as `13` in decimal
////
//// This is slightly different than the previous integer representation
//// conversion scheme we used for `GF(p^k)`.

//// <i><u>Exercise:</u></i> Convert 32 to a polynomial in `GF(2)[x]`, then convert to binary</br>

//// <p><i><u>Exercise:</u></i> Convert 31 to a polynomial in `GF(2)[x]`, then convert to binary</p>

//// <p><i><u>Exercise:</u></i> Convert the polynomial x^4 + x in `GF(2)[x]` to binary, then convert to decimal</p>

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GF(u64);

//// The basic utility functions are very straightforward now:

impl GF {
  pub fn new(val: u64) -> GF {
    // Sanity check!
    assert!((val as usize) < GF::number_of_elements());
    GF(val)
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

//// #### Addition in `GF(2^k)`
////
//// Recall that polynomial addition is pointwise addition of coefficents. Since addition in `GF(2)` is an XOR, it is also in `GF(2^k)`:
////
//// ```text
////          1x^2 + 1x + 1   =>    0b0111
//// + 1x^3 + 0x^2 + 1x + 1   => ^  0b1011
//// ----------------------      ---------
////   1x^3 + 1x^2 + 0x + 0  =>     0b1100
//// ```

//// In code, this is simply:

impl Add<GF> for GF {
  type Output = GF;
  fn add(self, rhs: GF) -> GF {
    GF::new(self.0 ^ rhs.0)
  }
}

//// #### Negation and Subtraction in `GF(2^k)`
////
//// To negate, we want the value b, such that a+b=0
////
//// <i><u>Exercise:</u></i> Look at the addition table for `GF(2)`. For each number, what is it's negation?
////
//// That's right, the negation of a number in `GF(2)` is itself. Negation is an identity operation! And, in `GF(2^k)` it is the same also.
////
//// So, we have the simple implementation for negation:

impl GF {
  pub fn negate(self) -> GF {
    self
  }
}

//// A fascinating consequence of this is that subtraction and addition are the same operation. Just an XOR!
////
//// Here, we'll write the subtraction implementation as we have for other fields. But, in many optimized implementations, the same routine as addition is used!

impl Sub<GF> for GF {
    type Output = GF;
    fn sub(self, rhs: GF) -> GF {
        self + rhs.negate() // or `self + rhs` or `self.0 ^ rhs.0`
    }
}

//// #### Multiplication in `GF(2^k)`
////
//// Now we come to multiplication, the tricky operation in polynomial fields.
////
//// For `GF(p^k)` we implemented multiplication as two-steps:
////
////   1. Convolution of coefficents yielding `2k` new coefficents
////   2. Polynomial modulus reduction to `k` coefficents
////
//// For `GF(2^k)` we could do exactly the same thing but we can simplify by better utilizing bitwise operations.
////
//// First let's recall the two operations we used in reduction for `GF(p^k)`:
////
//// - `shift`: Multiplying by a monomial in `(1, x, x^2, x^3, ...etc...)` is a simple shift of the coefficient array
//// - `scale`: Multiplying by any scalar is a simple pointwise-coefficient multiply.
////
//// For `GF(2^k)` these operations simplify:
////
//// - `shift`: Shifting the coefficent array is just performing a "bit-shift" (i.e. `<<` or `>>`)
//// - `scale`: Scaling means multiplying by 0 or 1, this is a "selection operation" (i.e. `if (not_zero) use_it()`)
////
//// <i><u>Exercise:</u></i> Convince yourself that these are equivalent operations for `GF(2^k)`?
////
////
//// #### Multiplication by shift and select
////
//// Another way to perform polynomial multiplication is to decompose it into a several simpler products using the distribution law.
////
//// Consider polynomials `a` and `b` and we wish to obtain `a*b`. If we expand polynomial `b`, we have:
////
//// ```text
//// a * b = a * (b_n x^n + ... + b_1 x + b_0)
////       = b_n * a * x^n + ... + b_1 * a * x + b_0 * a
//// ```
////
//// As discussed in the previous section, multiplying by a monomial (e.g. `x^k`) is a simple coefficent shift operation.
////
//// Thus, we'll use the notation:
//// ```text
//// a * x^k = a << k
//// ```
////
//// This gives us:
//// ```text
//// a * b = b_n * (a << n) + ... + b_1 * (a << 1) + b_0 * (a << 0)
//// ```
////
//// Now observe that for coefficents in `GF(2)`, `b_i` is in the set `{0, 1}`. This means that in the above decomposition, there is no coefficent multiplication. Instead, we either add the term to the result or we skip it (selection).
////
//// As an example, consider:
//// ```text
//// b = x^5 + x^3 + 1  (or 0b101001)
//// ```
////
//// Then, our decomposition is:
//// ```text
//// a * b = (a << 5) + (a << 3) + (a << 0)
//// ```
////
//// In other words, we simply sum all shifts, but skip anywhere the coefficent is 0.
////
//// In pseudocode, we have:
////
//// ```text
//// c = 0
//// for i in 0..k:
////   if extract_bit(b, i) == 1:
////     c = poly_add(c, a)   ## c ^= a
////   a = a << 1
//// return c
//// ```
////
//// Here, we iterate through each possible term (`a << i`) one at a time. Then, we select the terms that we need and add them.
////
//// <i><u>Exercise:</u></i> Think about this algorithm a bit. Convince yourself it is correct. Why is this a
////
//// <i><u>Exercise:</u></i> This algorithm appears to only consider `n` terms, but a convolution involves `O(n^2)` terms. What's going on here?
////
//// <i><u>Exercise:</u></i> How would this scale to very large `k`? (e.g. consider k = 1,000,000,000)
////
//// #### Applying the polynomial reduction
////
//// To complete the multiplication, we need to apply the polynomial reduction. We could do this in the same way as our `GF(p^k)` implementation, but
//// there's a better approach that allows us to integrate this operation into the above algorithm.
////
//// The first observation is that we never need to apply a reduction after an addition in `GF(2^k)`. This means that in the above algorithm, we only have to
//// worry about the `a << 1` monomial shifting (`a * x`). Instead of applying the reduction at the end, we can instead apply this after every `a << 1` computation.
////
//// At first glance, this doesn't seem like an improvement. Now we're doing `k` reductions instead of a single one (at the end)! But, it turns out each of
//// these `k` reductions is a much simper and faster operation.
////
//// If `a` is in `GF(2^k)`, then it can be represented by at most `k` bits. Thus, `a << 1` can be represented in
//// at most `k + 1` bits. This means that we only have to consider the most-significant-bit. If it is set, we simply need to subtract off `Q` (the irreducible polynomial) to obtain just `k` bits.
////
//// Let's do an example using:
//// ```text
//// K = 3
//// Q = x^3 + x + 1  (0b1011)
//// A = x^2 + x      (0b0110)
//// ```
////
//// Step by step:
////
//// | Iter # | `A`      | `A << 1` | `(A << 1) % Q`  |
//// |--------|----------|----------|-----------------|
//// |   `1`  | `0b0110` | `0b1100` | `0b0111`        |
//// |   `2`  | `0b0111` | `0b1110` | `0b0101`        |
//// |   `3`  | `0b0101` | `0b1010` | `0b0001`        |
//// |   `4`  | `0b0001` | `0b0010` | `0b0010`        |
//// |   `5`  | `0b0010` | `0b0100` | `0b0100`        |
//// |   `6`  | `0b0100` | `0b1000` | `0b0011`        |
//// |   `7`  | `0b0011` | `0b0110` | `0b0110`        |
////
//// As you can see, when the most-significant-bit is set (iters 1,2,3,6), we subtract off `Q` (XOR).
//// And, when the most-significant-bit is clear (iters 4,5,7), we do nothing.
////
//// <i><u>Exercise:</u></i> Work through an iteration table using A = x^2 + 1 (0b1001) instead.
////
//// Let's add this reduction step to out multiplication algorithm:
////
//// ```text
//// c = 0
//// for i in 0..k:
////   if extract_bit(b, i) == 1:
////     c = poly_add(c, a)   ## c ^= a
////   a = a << 1
////   if extract_bit(a, k) == 1:
////     a = poly_sub(a, q)
//// return c
//// ```
////
//// <i><u>Exercise:</u></i> Make sure you fully understand the algorithm before proceeding.
////
//// Okay, lets get back to coding!
////

fn extract_bit(n: u64, i: usize) -> u64 {
    (n >> i) & 1
}

impl Mul<GF> for GF {
    type Output = GF;
    fn mul(self, rhs: GF) -> GF {
        // First we unpack to get the raw u64 and we implement the algorithm
        // directly over the bits, rather than using the field add / sub operators.
        let mut a: u64 = self.0;
        let b: u64 = rhs.0;
        let mut c: u64 = 0;

        // Loop over each possible term
        for i in 0..K {
            if extract_bit(b, i) == 1 {
                c ^= a; // c = poly_add(c, a)
            }
            a <<= 1;
            if extract_bit(a, K) == 1 {
                a ^= Q; // a = poly_sub(a, Q)
            }
        }
        GF::new(c)
    }
}

//// Just like before, we will just use brute force to find inverses.
////
//// There are much better approaches, but we'll need a dedicated article to focus on those.

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
//// These are all quite similar to `GF(p^k)` if you'd rather skim them.
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
      return Err(format!(
        "Number too large, got {}, but limit is {}",
        num, limit
      ));
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
//// Note that these tests assume `GF(2^3)`. If you change the field, they are not expected to pass.
#[cfg(test)]
mod tests {
  use super::*;

  // Ensure we're using a specific field size for testing
  #[test]
  fn test_field_size() {
    // For K=3, we should have 2^3 = 8 elements
    assert_eq!(GF::number_of_elements(), 8);
  }

  // TEST: Verify basic field element construction
  #[test]
  fn test_element_construction() {
    for i in 0..8 {
      let element = GF::new(i);
      assert_eq!(element.value(), i);
    }
  }

  // TEST: We shouldn't be able to construct numbers out of the range
  #[should_panic]
  #[test]
  fn test_invalid_numbers() {
    GF::new(8); // With K=3, the largest valid value is 7
  }

  // TEST: Addition (which is XOR in GF(2^K))
  #[test]
  fn test_add() {
    assert_eq!(GF::new(0) + GF::new(0), GF::new(0));
    assert_eq!(GF::new(0) + GF::new(1), GF::new(1));
    assert_eq!(GF::new(1) + GF::new(1), GF::new(0));
    assert_eq!(GF::new(1) + GF::new(2), GF::new(3));
    assert_eq!(GF::new(3) + GF::new(4), GF::new(7));
    assert_eq!(GF::new(5) + GF::new(7), GF::new(2));
    assert_eq!(GF::new(6) + GF::new(3), GF::new(5));
  }

  // TEST: Negation (which is identity in GF(2^K))
  #[test]
  fn test_negation() {
    // In GF(2^K), negation is the identity function
    for i in 0..8 {
      assert_eq!(GF::new(i).negate(), GF::new(i));
    }
  }

  // TEST: Subtraction (which is identical to addition in GF(2^K))
  #[test]
  fn test_sub() {
    // In GF(2^K), addition and subtraction are the same operation
    // GF(2^3) | 0 - 0 = 0
    assert_eq!(GF::new(0) - GF::new(0), GF::new(0));

    // GF(2^3) | 0 - 1 = 1
    assert_eq!(GF::new(0) - GF::new(1), GF::new(1));

    // GF(2^3) | 1 - 1 = 0
    assert_eq!(GF::new(1) - GF::new(1), GF::new(0));

    // GF(2^3) | 3 - 2 = 1
    assert_eq!(GF::new(3) - GF::new(2), GF::new(1));

    // GF(2^3) | 7 - 5 = 2
    assert_eq!(GF::new(7) - GF::new(5), GF::new(2));
  }

  // TEST: Multiplication using the irreducible polynomial Q = x^3 + x + 1 (0b1011)
  #[test]
  fn test_mul() {
    assert_eq!(GF::new(0) * GF::new(0), GF::new(0));
    assert_eq!(GF::new(0) * GF::new(1), GF::new(0));
    assert_eq!(GF::new(1) * GF::new(1), GF::new(1));
    assert_eq!(GF::new(2) * GF::new(2), GF::new(4));
    assert_eq!(GF::new(3) * GF::new(3), GF::new(5));
    assert_eq!(GF::new(4) * GF::new(4), GF::new(6));
    assert_eq!(GF::new(5) * GF::new(6), GF::new(3));
    assert_eq!(GF::new(7) * GF::new(7), GF::new(3));
  }

  #[test]
  fn test_invert() {
    // 0 has no inverse
    assert!(matches!(GF::new(0).invert(), Err(_)));

    // 1 is its own inverse
    assert_eq!(GF::new(1).invert().unwrap(), GF::new(1));

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
        GF::new(i).invert().unwrap(),
        GF::new(inverse_table[i as usize])
      );

      // Verify that x * x^(-1) = 1
      assert_eq!(
        GF::new(i) * GF::new(inverse_table[i as usize]),
        GF::new(1)
      );
    }
  }

  #[test]
  fn test_div() {
    // Division by zero is an error
    assert!(matches!(GF::new(1) / GF::new(0), Err(_)));

    // 0 divided by anything non-zero is 0
    for i in 1..8 {
      assert_eq!(GF::new(0) / GF::new(i), Ok(GF::new(0)));
    }

    // Division by 1 is identity
    for i in 0..8 {
      assert_eq!(GF::new(i) / GF::new(1), Ok(GF::new(i)));
    }

    // For all non-zero a and b, a/b = a * inv(b)
    for a in 1..8 {
      for b in 1..8 {
        let expected = GF::new(a) * GF::new(b).invert().unwrap();
        assert_eq!(GF::new(a) / GF::new(b), Ok(expected));
      }
    }

    // Some specific test cases
    assert_eq!(GF::new(2) / GF::new(3), Ok(GF::new(7)));
    assert_eq!(GF::new(5) / GF::new(6), Ok(GF::new(4)));
    assert_eq!(GF::new(7) / GF::new(4), Ok(GF::new(3)));
  }
}
