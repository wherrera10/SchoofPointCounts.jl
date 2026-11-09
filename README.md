# SchoofPointCounts

Schoof's algorithm for counting prime field elliptic curve points. Uses Nemo's finite fields.

    schoof(A, B, q) -> Integer
    schoof(A::FinFieldElem, B::FinFieldElem) -> Integer

Return `#E(Fq)` for `y^2 = x^3 + A*x + B`, the usual prime order elliptical curve in cryptography.

The return value will be an integer, generally Int or BigInt depending on the size of the result.

The integer interface for `schoof` accepts a prime-power order `q = p^k` and interprets integer
coefficients in the prime subfield. The finite-field interface accepts coefficients in an
existing Nemo finite field.
