//// #### An implementation of GF(p) where p is a prime number
////
//// We will implement addition, subtraction, multiplication, and division using
//// operator-overloading so we can use the normal operators: +, -, *, /
use std::ops::{Add, Mul, Sub, Div};

//// Rust supports [Const Generics](https://practice.rs/generics-traits/const-generics.html),
//// but that adds a lot of extra syntax noise which may be distracting for readers without
//// a background in Rust or C++.
////
//// Instead, we'll just define a simple constant. If you'd like to explore a `GF(p)` field for a
//// different prime, just change the number below. You are responsible for ensuring that the number is actually
//// prime. By default, we'll use `GF(5)`:
pub const P: u8 = 5;

//// The type we use to represent a number in `GF(p)` is just an unsigned 8-bit integer (`u8`).
//// This design will allow for prime fields up to size `251` (the largest prime that fits
//// in 8-bits). We will also tell Rust that its okay to copy (`Clone` & `Copy`) and to compare
//// by two numbers via equality (`PartialEq`)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GF(u8);


//// Constructing a new number is `GF(p)` is easy: we just need to ensure it's within
//// the set `{0, 1, ..., p-1}`
impl GF {
  pub fn new(val: u8) -> GF {
    assert!(val < P);
    GF(val)
  }
  pub fn value(&self) -> u8 {
    self.0
  }
}

//// #### Addition and Subtraction
////
//// Now we implement addition via simple modular arithmetic: `(a + b) % p`.
//// But we need to take care to avoid 8-bit overflow, so we do the math in
//// 16-bit (`u16`) and then apply the modulus operation.
impl Add<GF> for GF {
  type Output = GF;
  fn add(self, rhs: GF) -> GF {
    let a = self.0 as u16;
    let b = rhs.0 as u16;
    let p = P as u16;
    GF::new(((a + b) % p) as u8)
  }
}

//// Negation (or additive-inverse) is fairly straight-forward also.
////
//// Given `a`, we want to solve for `neg(a)` in:
////
////   ```a + neg(a) == 0 (mod p)```
////
//// It seems like the answer would be:
////
////   ```neg(a) = p - a```
////
//// But if `a == 0`, then we'd compute `neg(a) = p` which is not a valid number.
////
//// This has an easy fix though:
////
////   ```neg(a) = (p - a) % p```
impl GF {
  pub fn negate(self) -> GF {
    GF::new((P - self.0) % P)
  }
}

//// Now, we can easily implement subtraction in terms of addition and negation
//// since: `a - b = a + (-b)`
impl Sub<GF> for GF {
  type Output = GF;
  fn sub(self, rhs: GF) -> GF {
    self + rhs.negate()
  }
}

//// #### Multiplication and Division
////
//// Now we implement multiplication via simple modular arithmetic: `(a * b) % p`.
//// But we need to take care to avoid 8-bit overflow, so we do the math in
//// 16-bit (`u16`) and then apply the modulus operation.
impl Mul<GF> for GF {
  type Output = GF;
  fn mul(self, rhs: GF) -> GF {
    let a = self.0 as u16;
    let b = rhs.0 as u16;
    let p = P as u16;
    GF::new(((a * b) % p) as u8)
  }
}

//// Multiplicative inverses are the trickiest operation in this field. We will
//// implement them by brute force. If we've constructed a proper field, then each
//// number will have an inverse. If we try all numbers, one will succeed.
////
//// Notice the `Result<>` type here. This is Rust's way of returning errors. We need to
//// use a `Result<>` here since the number 0 has no inverse and need to communicate that
//// the operation has failed.
////
//// Faster approaches exist (e.g. [The Extended Euclidean Algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm)), but for
//// a first implementation we will avoid adding that complexity. We'll discuss faster
//// methods in later sections.
impl GF {
  pub fn invert(self) -> Result<GF, String> {
    // Important: Zero has no inverse, it's invalid
    if self.0 == 0 {
      return Err("Zero has no inverse".to_string());
    }
    // Scan the numbers {1, 2, ..., P-1} until we find the inverse
    for x in 1..P {
      let candidate = GF::new(x);
      if self * candidate == GF::new(1) {
        return Ok(candidate); // Found!
      }
    }
    unreachable!("Every non-zero number has an inverse");
  }
}

//// Similairly to subtraction, we can implement division in terms of multiplication
//// and inversion since: `a / b = a * inv(b)`
////
//// Notice again that we use a `Result<>` to communicate that dividision by zero will fail
impl Div<GF> for GF {
  type Output = Result<GF, String>;
  fn div(self, rhs: Self) -> Result<GF, String> {
    // Important: Cannot divide by zero
    if rhs.0 == 0 {
      return Err("Cannot divide by zero".to_string());
    }
    Ok(self * rhs.invert().unwrap())
  }
}

//// #### Some final things
////
//// Rust as a language doesn't implicitly assume much about new types so we have to explicitly
//// tell it how to do a few more-or-less trivial things.
////
//// We need to tell Rust how to print these new numbers to the screen.
//// We will just print them out as ordinary integers.
impl std::fmt::Display for GF {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    write!(f, "{}", self.0)
  }
}

//// And, we need to tell Rust how to convert strings into our field's numbers.
impl std::str::FromStr for GF {
  type Err = String;
  fn from_str(s: &str) -> Result<GF, String> {
    let num: u8 = s.parse().map_err(|_| format!("Not an 8-bit integer"))?;
    // Return an error if the number is too big for the field
    if num >= P {
      return Err(format!("Number too large, got {}, but limit is {}", num, P-1));
    }
    Ok(GF::new(num))
  }
}

//// Finally, we'll tell rust that that our implementation can be treated as a Field in any
//// code that is written for Fields! This will come in handy for our calculator.
impl crate::field::Field for GF {
  fn number_of_elements() -> usize {
    P as usize
  }
}

//// #### Testing Time
////
//// Note that these tests assume `GF(5)`. If you change the size of the field, they are not expected to pass.
#[cfg(test)]
mod tests {
  use super::*;

  // TEST: We shouldn't be able to construct numbers out of the range
  #[should_panic]
  #[test]
  fn test_invalid_numbers() {
    GF::new(5);
  }

  // TEST: Addition
  #[test]
  fn test_add() {
    // GF(5) | 0 + 1 = 1
    assert_eq!(GF::new(0) + GF::new(1), GF::new(1));
    // GF(5) | 1 + 1 = 2
    assert_eq!(GF::new(1) + GF::new(1), GF::new(2));
    // GF(5) | 1 + 2 = 3
    assert_eq!(GF::new(1) + GF::new(2), GF::new(3));
    // GF(5) | 2 + 2 = 4
    assert_eq!(GF::new(2) + GF::new(2), GF::new(4));
    // GF(5) | 1 + 3 = 4
    assert_eq!(GF::new(1) + GF::new(3), GF::new(4));
    // GF(5) | 2 + 3 = 0
    assert_eq!(GF::new(2) + GF::new(3), GF::new(0));
    // GF(5) | 3 + 3 = 1
    assert_eq!(GF::new(3) + GF::new(3), GF::new(1));
    // GF(5) | 2 + 4 = 1
    assert_eq!(GF::new(2) + GF::new(4), GF::new(1));
    // GF(5) | 4 + 4 = 3
    assert_eq!(GF::new(4) + GF::new(4), GF::new(3));
  }

  // TEST: Subtraction
  #[test]
  fn test_sub() {
    // GF(5) | 1 - 0 = 1
    assert_eq!(GF::new(1) - GF::new(0), GF::new(1));
    // GF(5) | 0 - 1 = 4
    assert_eq!(GF::new(0) - GF::new(1), GF::new(4));
    // GF(5) | 1 - 1 = 0
    assert_eq!(GF::new(1) - GF::new(1), GF::new(0));
    // GF(5) | 1 - 2 = 4
    assert_eq!(GF::new(1) - GF::new(2), GF::new(4));
    // GF(5) | 2 - 2 = 0
    assert_eq!(GF::new(2) - GF::new(2), GF::new(0));
    // GF(5) | 1 - 3 = 3
    assert_eq!(GF::new(1) - GF::new(3), GF::new(3));
    // GF(5) | 2 - 3 = 4
    assert_eq!(GF::new(2) - GF::new(3), GF::new(4));
    // GF(5) | 3 - 3 = 0
    assert_eq!(GF::new(3) - GF::new(3), GF::new(0));
    // GF(5) | 2 - 4 = 3
    assert_eq!(GF::new(2) - GF::new(4), GF::new(3));
    // GF(5) | 4 - 4 = 0
    assert_eq!(GF::new(4) - GF::new(4), GF::new(0));
  }

  // TEST: Multiplication
  #[test]
  fn test_mul() {
    // GF(5) | 0 * 1 = 0
    assert_eq!(GF::new(0) * GF::new(1), GF::new(0));
    // GF(5) | 0 * 2 = 0
    assert_eq!(GF::new(0) * GF::new(2), GF::new(0));
    // GF(5) | 1 * 1 = 1
    assert_eq!(GF::new(1) * GF::new(1), GF::new(1));
    // GF(5) | 1 * 2 = 2
    assert_eq!(GF::new(1) * GF::new(2), GF::new(2));
    // GF(5) | 2 * 2 = 4
    assert_eq!(GF::new(2) * GF::new(2), GF::new(4));
    // GF(5) | 2 * 3 = 1
    assert_eq!(GF::new(2) * GF::new(3), GF::new(1));
    // GF(5) | 3 * 3 = 4
    assert_eq!(GF::new(3) * GF::new(3), GF::new(4));
    // GF(5) | 2 * 4 = 3
    assert_eq!(GF::new(2) * GF::new(4), GF::new(3));
    // GF(5) | 3 * 4 = 2
    assert_eq!(GF::new(3) * GF::new(4), GF::new(2));
    // GF(5) | 4 * 4 = 1
    assert_eq!(GF::new(4) * GF::new(4), GF::new(1));
  }


  // TEST: Division
  #[test]
  fn test_div() {
    // GF(5) | 0 / 1 = 0
    assert_eq!(GF::new(0) / GF::new(1), Ok(GF::new(0)));
    // GF(5) | 0 / 2 = 0
    assert_eq!(GF::new(0) / GF::new(2), Ok(GF::new(0)));
    // GF(5) | 1 / 0 = ERROR
    assert!(matches!(GF::new(1) / GF::new(0), Err(_)));
    // GF(5) | 2 / 0 = ERROR
    assert!(matches!(GF::new(2) / GF::new(0), Err(_)));
    // GF(5) | 1 / 1 = 1
    assert_eq!(GF::new(1) / GF::new(1), Ok(GF::new(1)));
    // GF(5) | 1 / 2 = 3
    assert_eq!(GF::new(1) / GF::new(2), Ok(GF::new(3)));
    // GF(5) | 2 / 2 = 1
    assert_eq!(GF::new(2) / GF::new(2), Ok(GF::new(1)));
    // GF(5) | 2 / 3 = 4
    assert_eq!(GF::new(2) / GF::new(3), Ok(GF::new(4)));
    // GF(5) | 3 / 3 = 1
    assert_eq!(GF::new(3) / GF::new(3), Ok(GF::new(1)));
    // GF(5) | 2 / 4 = 3
    assert_eq!(GF::new(2) / GF::new(4), Ok(GF::new(3)));
    // GF(5) | 3 / 4 = 2
    assert_eq!(GF::new(3) / GF::new(4), Ok(GF::new(2)));
    // GF(5) | 4 / 4 = 1
    assert_eq!(GF::new(4) / GF::new(4), Ok(GF::new(1)));
  }
}
