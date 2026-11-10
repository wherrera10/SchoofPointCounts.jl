# SchoofPointCounts

SchoofPointCounts implements Schoof's algorithm for counting points on elliptic
curves over finite fields, using Nemo.jl.

```julia
schoof(A, B, q) -> Integer
schoof(A::FinFieldElem, B::FinFieldElem) -> Integer
```

Both methods return the number of points on the curve
`y^2 = x^3 + A*x + B` over a finite field. The integer method takes a
prime-power field order `q = p^k` and integer coefficients `A` and `B`, which
are interpreted in the prime subfield. The two-argument method takes
coefficients from the same existing Nemo finite field. The field
characteristic must be greater than 3, and the curve must be nonsingular.

For example, count the points on `y^2 = x^3 + x + 1` over `F_5`:

```julia
julia> using SchoofPointCounts

julia> schoof(1, 1, 5)
9
```

## Ecosystem notes 

Nemo.jl does include elliptic-curve research code, as does the separate EllipticCurves.jl project. Nemo's research 
code included Frobenius computations used in Schoof–Atkin–Elkies point counting, but that is not the same as the
general-purpose Schoof implementation of `schoof` implemented here. OSCAR.jl uses calls to the PARI library to
provide a call to a Schoof–Atkin–Elkies algorithm for its order() function.

## Installation

You may be able to install the package from the Julia REPL in the usual way with  `]add SchoofPointCounts`,  
or to install from the current master copy:

```julia
using Pkg
Pkg.add(url = "https://github.com/wherrera10/SchoofPointCounts.jl")
```

Requires Julia 1.10 or later, plus Primes.jl and [Nemo.jl](https://github.com/Nemocas/Nemo.jl), which should be added automatically.
