# SchoofPointCounts

[![Documentation Status](https://img.shields.io/badge/docs-latest-blue.svg)](https://wherrera10.github.io/SchoofPointCounts.jl/)
[![CI](https://github.com/wherrera10/SchoofPointCounts.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/wherrera10/SchoofPointCounts.jl/actions/workflows/CI.yml)


SchoofPointCounts implements Schoof's algorithm for counting points on elliptic
curves over finite fields, using Nemo.jl.


This package exports a single function: `schoof`, with two methods:

The first method takes three integer arguments, and returns an Int or a BigInt, depending on the size of the return value:
```julia
schoof(A, B, q)::Integer # All arguments are handled as integers
```

The second method assumes the user has defined a finite field in Nemo.jl and takes two Nemo.FinFieldElem values in that same field:

```julia
schoof(A::Nemo.FinFieldElem, B::Nemo.FinFieldElem)::Integer # Both arguments are Nemo.FinFieldElem type values
```

Both methods return an integer (either Int or if too large for Int a BigInt) which is the number of points on the curve
`y^2 = x^3 + A*x + B` 
over a finite field. The integer method takes a prime-power field order integer 
`q = p^k` 
and integer coefficients `A` and `B`, which are interpreted in the prime subfield. The two-argument method takes coefficients
from the same existing Nemo finite field. The field characteristic must be greater than 3, and the curve must be nonsingular,
or an argument error will be thrown.

For example, to count the points on elliptic curve

```julia
y^2 = x^3 + x + 1
```

with prime p = 5 and k = 1, so that q is 5:

```julia
julia> using SchoofPointCounts

julia> schoof(1, 1, 5)
9
```

Note that the return value is the number of points on the projective elliptic curve; its affine points are supplemented by the single point at infinity, so the point at infinity is included in the returned count.

## Ecosystem notes

Nemo.jl does include elliptic-curve research code, as does the separate EllipticCurves.jl project. Nemo's research code includes Frobenius computations used in Schoof–Atkin–Elkies (SAE) point counting, but that is not the same as the general-purpose Schoof implementation of `schoof` implemented here. OSCAR.jl may use calls to the SAE code in the PARI library or to Hecke.jl for its order() function. Hecke.jl, part of the OSCAR project, uses a linear algebra-based general algorithm for its order() function. Benchmarking suggests Hecke.jl is 10 to 40% faster than schoof() for integers from 16 to 96 bits, but that as big integers expand past 100 bits in size, Hecke.jl seems to slow more than schoof(), so that, with cryptographically significant random p curves, schoof() is often, but not always, faster than Hecke.jl.

## Installation

You may be able to install the package from the Julia REPL in the usual way with  `]add SchoofPointCounts`, or to install from the current GitHub archive:

```julia
using Pkg
Pkg.add(url = "https://github.com/wherrera10/SchoofPointCounts.jl")
```

Requires Julia 1.10 or later, plus Primes.jl and [Nemo.jl](https://github.com/Nemocas/Nemo.jl), which should be added automatically.
