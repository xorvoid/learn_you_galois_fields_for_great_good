use std::io::{self, Read};
use learn_you_galois_fields_for_great_good::crc32;

fn main() -> io::Result<()> {
    let mut buffer = String::new();
    io::stdin().read_to_string(&mut buffer)?;

    let cksum = crc32::cksum(buffer.as_bytes());
    println!("{} {}", cksum, buffer.len());

    Ok(())
}
