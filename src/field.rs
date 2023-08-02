//// #### A type definition for a Field
////
//// A Field in math requires addition, subtraction, multiplication, and division.
//// But we will also require copying, equality comparisions, parsing from a string, and converting to a string.
use std::ops::{Add, Mul, Sub, Div};

pub trait Field:
  Add<Output=Self> +
  Mul<Output=Self> +
  Sub<Output=Self> +
  Div<Output=Result<Self, String>> +
  Copy + Clone + PartialEq +
  std::str::FromStr<Err=String> + std::fmt::Display
{
  fn number_of_elements() -> usize;
}
