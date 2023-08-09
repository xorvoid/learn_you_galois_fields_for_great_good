//// #### A Calculator for GF(p)
////
//// Print out the addition and multiplication tables and then provide an
//// interactive calculator for the user.
////
use learn_you_galois_fields_for_great_good::gf_p_k::{GF, P, K};
use learn_you_galois_fields_for_great_good::calc;

fn main() {
  // Simply call into the generic interactive calculator using our field
  let name = format!("GF({}^{})", P, K);
  calc::interactive_calculator::<GF>(&name);
}
