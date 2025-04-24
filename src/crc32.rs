//// #### Implementing CRC-32 Checksum
////
//// CRC-32 is a checksum that reduces n-bits into 32-bits. Viewed as an operation in an abstract algebra, we are reducing a polynomial in `GF(2^n)` to a polynomial in `GF(2^32)`.
////
//// The polynomial used for this reduction is:
//// ```text
//// x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x + 1
//// ```
//// In hexadecimal this polynomial is `0x104c11db7`, but its typical to drop the leading term (`x^32`) so the polynomial number fits
//// in a 32-bit unsigned integer. After applying that adjustment, we get:
const Q: u32 = 0x04c11db7;

//// Our approach to implementing CRC-32 will be similar to a hardware implementation. This is not the fastest approach
//// in software, but it is the simplest and easiest to understand. It also helps us understand why CRC codes are so
//// nice to implement in hardware.
////
//// To start, we'll define a type for the current checksum state, which should be exactly 32-bits
#[derive(Debug, Clone, Copy, PartialEq)]
struct CRC(u32);

impl CRC {
    fn new() -> CRC {
        CRC(0)
    }
}

//// #### Decomposing into shifts
////
//// Consider the polynomial:
//// ```text
//// x^4 + x^3 + x + 1
//// ```
////
//// We can decompose this polynomial as follows:
//// ```text
//// (x^3 + x^2 + 1)*x + 1
//// ```
////
//// As we have discussed in the previous section, multiplying by `x` in `GF(2^k)` is the same as shifting bits to the left by one. Rewriting, we have:
//// ```text
//// (x^3 + x^2 + 1)<<1 + 1
//// ```
////
//// Decomposing repeatedly:
//// ```text
//// ((x^2 + x)<<1 + 1)<<1 + 1
//// (((x + 1)<<1 + 0)<<1 + 1)<<1 + 1
//// ((((1)<<1 + 1)<<1 + 0)<<1 + 1)<<1 + 1
//// ```

//// Thinking about this as raw bit operations, notice that we can build up any polynomial by left-shifting and adding (xor-ing) new bits one-at-a-time. In this particular case, we perform:
//// ```text
//// xor 1
//// left-shift
//// xor 1
//// left-shift
//// xor 0
//// left-shift
//// xor 1
//// ```

//// <i><u>Exercise:</u></i> Decompose `x^5 + x^2 + x` into shifts and bit insertions.

//// #### Reducing the left-shift (`a << 1`)

//// In the implementation of `GF(2^k)` we found that we needed to reduce `a << 1` from a `k+1` bit number to a `k` bit number. For CRC, we need to do the same thing. Whenever, the kth bit is set, we need to subtract (XOR) away `Q`.

//// Let's modify the sequence above to include these reduce operations:

//// ```text
//// xor 1
//// left-shift
//// if most-significant-bit set, xor Q
//// xor 1
//// left-shift
//// if most-significant-bit set, xor Q
//// xor 0
//// left-shift
//// if most-significant-bit set, xor Q
//// xor 1
//// ```

//// <i><u>Exercise:</u></i> Decompose `x^5 + x^2 + x` into shifts and bit insertions and reductions in `GF(2^3)`

//// As you can see, a CRC code can be implemented very simply by just shifting the bits and applying xors. By reducing
//// every time the MSB is set, we ensure that we always have a number in the checksum field `GF(2^p)`. When the input
//// stream runs out of bits, we have our checksum!

//// Let's implement it:

impl CRC {
    fn append_bit(&mut self, bit: u8) {
        // Extract the most-significant-bit (MSB), we will need to check it after shifting
        let msb = (self.0 >> 31) & 1;

        // Shift the bits (notice: this overwrites the old MSB)
        let mut shifted = self.0 << 1;

        // Check if we need to reduce. If the MSB was set, then we have a 33-bit number
        // and we need to reduce it back to 32-bits using Q (notice: Q doesn't include
        // the MSB. Since the MSB was shifted off, it has been implicitly cleared already)
        if msb == 1 { shifted ^= Q; }

        // Add the next bit and save the result
        self.0 = shifted ^ (bit as u32);
    }
}

//// And that's the core operation of CRC-32, the rest is just supporting details!
////
//// <i><u>Exercise:</u></i> Really make sure you follow this algorithm. It can be helpful to draw diagrams
//// of bits shifting, new bits entering, and Q getting XOR'd. You can use the following parameters to construct
//// a simpler example: `p = 3, Q = x^3 + x + 1 (decimal: 11, binary: 1011)`
////
//// Next up, we need a method to append an entire byte, rather than one bit. This is pretty simple. We will
//// just iterate each bit from MSB to LSB, inserting each bit:
impl CRC {
    fn append_byte(&mut self, byte: u8) {
        // For each byte, just append each bit from MSB to LSB
        for i in (0..8).rev() {
            let bit = (byte >> i)&1;
            self.append_bit(bit);
        }
    }
}

//// #### A routine for the `cksum` tool
//// We'll now use the above routines to construct the same checksum as the Unix `cksum` tool. This tool is specified
//// by the POSIX Standard, IEEE Std 1003.1-2017 ([here](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/cksum.html)).
////
//// For the most part, we'll simply be calling our `append_byte()` routine, but there are a few additional quirks that we'll explain as we go:

pub fn cksum(dat: &[u8]) -> u32 {
    // Initialize our CRC state
    let mut crc = CRC::new();

    // Insert all of the data bytes in order
    for b in dat {
        crc.append_byte(*b);
    }

    // Quirk #1: The `cksum` tool appends the data length also. There is no fixed size
    // for the length field. Instead, we append data bytes from LSB to MSB. We stop
    // appending bytes when the MSB becomes zero.
    let mut len: u64 = dat.len() as u64;
    loop {
        // Append the lowest byte
        crc.append_byte(len as u8);

        // Shift the byte off
        len >>= 8;

        // When the remaining bytes are zero, we're done!
        if len == 0 {
            break;
        }
    }

    // Quirk #2: The `cksum` tool always appends an empty 32-bits after everything.
    // The main reason is to force each bit to cause a reduction (if required).
    // Each reduction causes bit-mixing/scrambling throughout the checksum.
    for _ in 0..32 {
        crc.append_bit(0);
    }

    // Quirk #3: The `cksum` tool inverts all the bits in the final result. The author
    // doesn't have an immediate justification for this operation. (Notice: In Rust
    // code, '!crc' is bitwise inversion. In C code and many related languages, this
    // would be '~crc' instead)
    !crc.0
}

//// #### Testing Time
////
//// Here are some tests that verify our implementation matches the values computed by the `cksum` tool. We also have
//// a test here demonstrating a desired property of a good checksum. Notably, if one bit changes, we want many
//// checksum bits to change.
////
#[cfg(test)]
mod test {
    use super::*;

    // Check some well-known strings for the right CRC-32 value
    #[test]
    fn test_strings() {
        assert_eq!(cksum(b""), 0xffffffff);
        assert_eq!(cksum(b"a"), 0x48c279fe);
        assert_eq!(cksum(b"123456789"), 0x377a6011);
        assert_eq!(cksum(b"finite fields are super fun when you really understand them!"), 0xa695ef0f);
    }

    // Test that a single flipped bit causes the checksum to flip many bits!
    #[test]
    fn test_one_flipped_bit() {
        let a: u64 = 0x42d5151330d94a84;
        let b = a ^ (1u64 << 30);  // flip the 30th bit

        // One bit flip causes many bits to flip in the checksum, very good!
        assert_eq!(cksum(&a.to_le_bytes()), 0xd146283d);
        assert_eq!(cksum(&b.to_le_bytes()), 0x01c33b8f);
    }
}

//// <i><u>Exercise:</u></i> Why does `cksum(b"")` result in `0xffffffff` when no set bits were appended ?

//// <i><u>Exercise:</u></i> Why does `cksum(b"a")` result in many set bits even though "a" is only one byte?

//// <i><u>Exercise:</u></i> Can you explain why altering a single input bit can significantly change the checksum?
