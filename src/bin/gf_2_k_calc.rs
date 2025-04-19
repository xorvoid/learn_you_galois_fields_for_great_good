//// #### A Calculator for GF(p)
////
//// Print out the addition and multiplication tables and then provide an
//// interactive calculator for the user.
////
use learn_you_galois_fields_for_great_good::calc;
use learn_you_galois_fields_for_great_good::gf_2_k::{GF, K};

fn main() {
  // Simply call into the generic interactive calculator using our field
  let name = format!("GF({}^{})", 2, K);
  calc::interactive_calculator::<GF>(&name);
}
