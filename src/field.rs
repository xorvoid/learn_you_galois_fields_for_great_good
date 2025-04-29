//// #### A type definition for a Field
////
//// A Field in math requires addition, subtraction, multiplication, and division.
//// But we will also require copying, equality comparisons, parsing from a string, and converting to a string.
//// These additional properties will make it easier to write generic code that will work over any Field!
use std::ops::{Add, Mul, Sub, Div};

pub trait Field:
  Add<Output=Self> +
  Mul<Output=Self> +
  Sub<Output=Self> +
  Div<Output=Self> +
  Copy + Clone + PartialEq +
  std::str::FromStr<Err=String> + std::fmt::Display
{
  fn number_of_elements() -> usize;
}
