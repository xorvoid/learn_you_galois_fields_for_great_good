//// #### An implementation of GF(256) based on the generic GF(2^k)
//// To understand this implementation, you may wish to refer to `src/gf_2_k.rs`
//// There are a few minor changes here, but the core remains identical.

use std::ops::{Add, Div, Mul, Sub};

pub const K: usize = 8;
pub const Q: u64 = 283;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GF(u8);

impl GF {
  pub fn new(val: u8) -> GF {
      // Sanity check!
      assert!((val as usize) < GF::number_of_elements());
      GF(val)
  }

  pub fn number_of_elements() -> usize {
      (1 << K) as usize
  }

  pub fn value(&self) -> u8 {
    self.0
  }
}

impl Add<GF> for GF {
  type Output = GF;
  fn add(self, rhs: GF) -> GF {
    GF::new(self.0 ^ rhs.0)
  }
}

impl GF {
  pub fn negate(self) -> GF {
    self
  }
}

impl Sub<GF> for GF {
    type Output = GF;
    fn sub(self, rhs: GF) -> GF {
        GF::new(self.0 ^ rhs.0)
    }
}

fn extract_bit(n: u64, i: usize) -> u64 {
    (n >> i) & 1
}

impl Mul<GF> for GF {
    type Output = GF;
    fn mul(self, rhs: GF) -> GF {
        // First we unpack to get the raw u64 and we implement the algorithm
        // directly over the bits, rather than using the field's add/sub operators.
        let mut a: u64 = self.0 as u64;
        let b: u64 = rhs.0 as u64;
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
        GF::new(c as u8)
    }
}

static INVERSE_LUT: std::sync::OnceLock<Vec<GF>> = std::sync::OnceLock::new();

impl GF {
    fn get_inverse_lut() -> &'static [GF] {
        INVERSE_LUT.get_or_init(|| {
            // Build up the inverse table using brute force
            let mut lut = vec![];
            lut.resize(GF::number_of_elements(), GF::new(0));

            // Find the inverse for each of the numbers {1, 2, ..., N-1}
            for x in 1..GF::number_of_elements() {
                // Scan the numbers {1, 2, ..., N-1} until we find the inverse
                let x = GF::new(x as u8);
                let mut found = false;
                for y in 1..GF::number_of_elements() {
                    let y = GF::new(y as u8);
                    if x * y == GF::new(1) {
                        lut[x.0 as usize] = y;
                        found = true;
                        break;
                    }
                }
                if !found {
                    unreachable!("Every non-zero number has an inverse");
                }
            }

            lut
        })
    }

    pub fn invert(self) -> Result<GF, String> {
        // Important: Zero has no inverse, it's invalid
        if self == GF::new(0) {
            return Err("Zero has no inverse".to_string());
        }
        // Perform a lookup in the pre-computed table
        Ok(GF::get_inverse_lut()[self.0 as usize])
    }
}

impl GF {
    fn checked_div(self, rhs: Self) -> Option<GF> {
        // Important: Cannot divide by zero
        if rhs == GF::new(0) {
            return None;
        }
        Some(self * rhs.invert().unwrap())
    }
}

impl Div<GF> for GF {
  type Output = GF;
  fn div(self, rhs: Self) -> GF {
    self.checked_div(rhs).unwrap()
  }
}

impl std::fmt::Display for GF {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    write!(f, "{}", self.value())
  }
}

impl std::str::FromStr for GF {
  type Err = String;
  fn from_str(s: &str) -> Result<GF, String> {
    let num: u64 = s.parse().map_err(|_| format!("Not an 64-bit integer"))?;
    // Return an error if the number is too big for the field
    let limit = GF::number_of_elements() as u64;
    if num >= limit {
      return Err(format!(
        "Number too large, got {}, but limit is {}",
        num, limit-1
      ));
    }
    Ok(GF::new(num as u8))
  }
}

//// And telling Rust that we built a Field type:
impl crate::field::Field for GF {
  fn number_of_elements() -> usize {
    GF::number_of_elements()
  }
}

impl GF {
    pub fn initialize_all_lookup_tables() {
        GF::get_inverse_lut();
    }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_add() {
    assert_eq!(GF::new(0) + GF::new(0), GF::new(0));
    assert_eq!(GF::new(0) + GF::new(1), GF::new(1));
    assert_eq!(GF::new(1) + GF::new(0), GF::new(1));
    assert_eq!(GF::new(1) + GF::new(1), GF::new(0));
    assert_eq!(GF::new(45) + GF::new(67), GF::new(110));
  }

  #[test]
  fn test_sub() {
    assert_eq!(GF::new(0) - GF::new(0), GF::new(0));
    assert_eq!(GF::new(0) - GF::new(1), GF::new(1));
    assert_eq!(GF::new(1) - GF::new(0), GF::new(1));
    assert_eq!(GF::new(1) - GF::new(1), GF::new(0));
    assert_eq!(GF::new(45) - GF::new(67), GF::new(110));
  }

  #[test]
  fn test_mul() {
    assert_eq!(GF::new(0) * GF::new(0), GF::new(0));
    assert_eq!(GF::new(0) * GF::new(1), GF::new(0));
    assert_eq!(GF::new(1) * GF::new(0), GF::new(0));
    assert_eq!(GF::new(1) * GF::new(1), GF::new(1));
    assert_eq!(GF::new(5) * GF::new(1), GF::new(5));
    assert_eq!(GF::new(1) * GF::new(5), GF::new(5));
    assert_eq!(GF::new(1) * GF::new(5), GF::new(5));
    assert_eq!(GF::new(3) * GF::new(5), GF::new(15));
    assert_eq!(GF::new(50) * GF::new(5), GF::new(250));
    assert_eq!(GF::new(5) * GF::new(50), GF::new(250));
    assert_eq!(GF::new(4) * GF::new(64), GF::new(27));
    assert_eq!(GF::new(5) * GF::new(64), GF::new(91));
    assert_eq!(GF::new(6) * GF::new(50), GF::new(172));
    assert_eq!(GF::new(255) * GF::new(255), GF::new(19));
    assert_eq!(GF::new(128) * GF::new(2), GF::new(27));
    assert_eq!(GF::new(175) * GF::new(98), GF::new(1));
  }

  #[test]
  fn test_inv() {
    assert!(GF::new(0).invert().is_err());
    assert_eq!(GF::new(1).invert().unwrap(), GF::new(1));
    assert_eq!(GF::new(2).invert().unwrap() * GF::new(2), GF::new(1));
    assert_eq!(GF::new(83).invert().unwrap() * GF::new(83), GF::new(1));
    assert_eq!(GF::new(123) * GF::new(123).invert().unwrap(), GF::new(1));
  }

  #[test]
  fn test_div() {
    assert_eq!(GF::new(1) / GF::new(1), GF::new(1));
    assert_eq!(GF::new(175) / GF::new(1), GF::new(175));
    assert_eq!(GF::new(175) / GF::new(2), GF::new(218));
    assert_eq!(GF::new(1) / GF::new(175), GF::new(98));
    assert_eq!((GF::new(1) / GF::new(175)) * GF::new(175), GF::new(1));
    assert_eq!(GF::new(1) / GF::new(175), GF::new(175).invert().unwrap());
  }

  #[should_panic]
  #[test]
  fn test_div_zero() {
    let _ = GF::new(1) / GF::new(0);
  }
}
