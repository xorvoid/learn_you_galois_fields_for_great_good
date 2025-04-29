//// #### A Calculator for GF(256)
////
//// Provide an interactive calculator for the user.
////
use learn_you_galois_fields_for_great_good::gf_256::GF;
use learn_you_galois_fields_for_great_good::calc;

fn main() {
  // Simply call into the generic interactive calculator using our field
  calc::interactive_calculator::<GF>("GF(256)", false);
}
