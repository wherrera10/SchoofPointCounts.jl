# SchoofPointCounts

Schoof's algorithm for counting prime field elliptic curve points. Uses Nemo's finite fields.

    schoof(A, B, q) -> Integer
    schoof(A::FinFieldElem, B::FinFieldElem) -> Integer

Return `#E(Fq)` for `y^2 = x^3 + A*x + B`, the usual prime order elliptical curve in cryptography.

The return value will be an integer, generally Int or BigInt depending on the size of the result.

The integer interface for `schoof` accepts a prime-power order `q = p^k` and interprets integer
coefficients in the prime subfield. The finite-field interface accepts coefficients in an
existing Nemo finite field.

## Ecosystem notes 

Nemo.jl does include elliptic-curve research code, as does the separate EllipticCurves.jl project. Nemo's reasearch has included
Frobenius computations used in Schoof–Atkin–Elkies point counting, but that is not the same as the general-purpose Schoof
implementation of `schoof` implemented here.

## Installation

You may be able to install the package from the Julia REPL in the usual way with  `]add SchoofPointCounts`,  
or to install from the current master copy:

```julia
using Pkg
Pkg.add(url = "https://github.com/wherrera10/SchoofPointCounts.jl")
```

Requires Julia 1.10 or later, plus [Nemo.jl](https://github.com/Nemocas/Nemo.jl), which should be added automatically.
