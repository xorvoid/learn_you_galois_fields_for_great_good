//// #### An implementation of Computer Science Fields GF(2^k)
////
//// After implementing the general `GF(p^k)` Fields, why would we implement a special case?
//// Can't we just use that implemention, configured for `p = 2`?
////
//// Fair question. You could indeed use it. But, often it's a bit overkill and significantly
//// better implementations exist for `GF(2^k)`
////
//// Let us explain.
////
//// Recall that we used a vector representation for `GF(p^k)`. In computer memory, this was an
//// array of bytes: one byte per coefficient. What happens when we use `p = 2`?  Well, it becomes
//// a vector of bits (0 or 1). And, we have a very efficient way to store that vector in a computer:
////
//// *An Ordinary Binary Number*
////
//// Consider an 8-bit unsigned integer (`u8`), we can view this as a *vector of 8 bits*.
//// Similarly, a 16-bit unsigned integer (`u8`) can be viewed as a *vector of 16 bits*. And so on!
////
//// Nifty, huh?
////
//// This is why I call these fields the *Computer Science Fields*. These map extremely well onto
//// ordinary binary-state transister computers. And, it turns out that common bitwise operations
//// are useful as well! We'll discuss these as we implement `GF(2^k)`
////
