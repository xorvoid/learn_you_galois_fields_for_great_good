//// Here is a simple tool to show the various representations of some number in `GF(p)[x]`.
//// It accepts a coefficient field prime p and number in decimal-representation n and displays
//// it in all the represeantations we've discussed.
////
//// The most consistent and useful representation for doing calculations is the vector format.
////
//// The polynomial representation is more verbose and inconsistent.
//// Thus, we will generally avoid using it directly. But, it can be helpful to drop into polynomial
//// arithmetic to understand something. In such cases, this tool can be helpful.
////
//// <div style="display: none">
fn main() {
  let args: Vec<_> = std::env::args().collect();
  if args.len() != 3 {
    eprintln!("usage: {} <prime> <number>", args[0]);
    std::process::exit(1);
  }

  let p: u64 = args[1].parse().unwrap();
  let mut num: u64 = args[2].parse().unwrap();

  println!("Using GF({})[x]:", p);
  println!("  Number:     {}", num);

  // Decompose the number into a vector of coefficients
  let mut coeffs = vec![];
  if num == 0 { // Handle special case
    coeffs.push(0);
  }
  while num > 0 {
    coeffs.push(num % p);
    num /= p;
  }

  // Loop and print each coeff of the vector (in reverse order)
  let mut first_term = true;
  print!("  Vector:     [");
  for i in (0..coeffs.len()).rev() {
    // Special case: insert ',' punctuation only between present terms
    if !first_term {
      print!(", ");
    }
    first_term = false;
    // General case: print number
    print!("{}", coeffs[i]);
  }
  println!("]");

  // Loop and print each term of the polynomial (in reverse order)
  let mut first_term = true;
  print!("  Polynomial: ");
  for i in (0..coeffs.len()).rev() {
    let c = coeffs[i];

    // Special case: ignore terms that aren't present
    if c == 0 {
      continue;
    }

    // Special case: insert '+' punctuation only between present terms
    if !first_term {
      print!(" + ");
    }
    first_term = false;

    // Special case: handle 1's coeff
    if i == 0 {
      print!("{}", c);
      continue;
    }

    // General: handle x, x^2, x^3, etc

    // Special case: don't print the '1' coeff, it is implied
    if c != 1 {
      print!("{}", c);
    }

    print!("x");

    // Special case: don't print 'x^1', it is implied by 'x'
    if i != 1 {
      print!("^{}", i);
    }
  }
  println!("");
}
//// </div>
