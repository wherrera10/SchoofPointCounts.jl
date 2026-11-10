"""
    SchoofPointCounts

Schoof's algorithm for counting points on `y^2 = x^3 + A*x + B` over a
finite field of characteristic greater than 3. The integer interface accepts
a prime-power order `q = p^k` and interprets integer coefficients in the
prime subfield. The finite-field interface accepts coefficients in an
existing Nemo finite field.
"""
module SchoofPointCounts

using Nemo
using Primes

export schoof


"""
    schoof(A, B, q) -> Integer
    schoof(A::FinFieldElem, B::FinFieldElem) -> Integer

Return `#E(Fq)` for `y^2 = x^3 + A*x + B`.
"""
function schoof(A::Integer, B::Integer, q::Integer)::Integer
    q > 3 || throw(ArgumentError("field order must exceed 3"))
    fac = Primes.factor(q)
    length(fac) == 1 || throw(ArgumentError("q = $q must be a prime power"))
    p, k = first(fac)
    p > 3 || throw(ArgumentError(
        "this implementation requires field characteristic greater than 3"
    ))
    Fq = finite_field(p, k, "a")[1]
    return schoof(Fq(A), Fq(B))
end

function schoof(A::FinFieldElem, B::FinFieldElem)::Integer
    Fq = parent(A)
    parent(B) == Fq || throw(ArgumentError("A and B must lie in the same field"))
    p = BigInt(characteristic(Fq))
    q = BigInt(order(Fq))
    p > 3 || throw(ArgumentError(
        "this implementation requires field characteristic greater than 3"
    ))
    iszero(4A^3 + 27B^2) &&
        throw(ArgumentError("the curve is singular (zero discriminant)"))

    R, X = polynomial_ring(Fq, "X")
    f = X^3 + A*X + B
    S, _ = auxiliaryprimes(p, q)
    g = gcd(powermod(X, q, f) - X, f)
    congruences = Tuple{Int,Int}[(2, isone(g) ? 1 : 0)]

    for ℓ in S
        ψ = divisionpolynomial(ℓ, A, B, R)
        tl = 0
        while true
            try
                tl, _ = tracemodℓ(ℓ, q, R(A), X, f, ψ)
                break
            catch e
                e isa NotInvertible || rethrow()
                # A nonunit identifies a factor on which the attempted
                # quotient-ring operation is not valid. Refine the modulus.
                degree(e.g) < degree(ψ) ||
                    error("unexpected zero divisor modulo ψ_$ℓ")
                degree(e.g) > 0 ||
                    error("noninvertibility did not produce a proper factor modulo ψ_$ℓ")
                ψ = e.g * inv(leading_coefficient(e.g))
            end
        end
        push!(congruences, (ℓ, mod(tl, ℓ)))
    end

    t = reconstruct_trace(congruences)
    result = q + 1 - t
    return typemin(Int) <= result <= typemax(Int) ? Int(result) : result
end

"""
    auxiliaryprimes(p, q) -> (S, M)

Return increasing odd primes `S`, excluding `p`, and `M = 2*prod(S)`,
stopping once `M^2 > 16q`, so `M > 4sqrt(q)` exactly.
"""
function auxiliaryprimes(p, q)
    S = Int[]
    M = big(2)
    ℓ = 3
    while M^2 <= 16 * big(q)
        ℓ = Primes.nextprime(ℓ)
        if ℓ != p
            push!(S, ℓ)
            M *= ℓ
        end
        ℓ += 1
    end
    return S, M
end

"""
    reconstruct_trace(congruences) -> BigInt

Combine pairwise-coprime congruences `(t ≡ r mod m)` by CRT and return
the representative in `(-M/2, M/2]`, where `M` is the product of moduli.
"""
function reconstruct_trace(congruences)
    isempty(congruences) &&
        throw(ArgumentError("at least one congruence is required"))
    moduli = BigInt[c[1] for c in congruences]
    all(m -> m > 1, moduli) ||
        throw(ArgumentError("CRT moduli must be integers greater than one"))
    for i in eachindex(moduli), j in firstindex(moduli):(i-1)
        gcd(moduli[i], moduli[j]) == 1 ||
            throw(ArgumentError("CRT moduli must be pairwise coprime"))
    end

    M = prod(moduli)
    t = BigInt(0)
    for (m0, r0) in congruences
        m, r = BigInt(m0), BigInt(r0)
        Mi = M ÷ m
        t += mod(r, m) * Mi * invmod(Mi, m)
    end
    t = mod(t, M)
    return t > M ÷ 2 ? t - M : t
end

"""
    NotInvertible

Raised when an element is not invertible modulo `ψ`. `g` stores the gcd
that identifies a nontrivial factor, or `ψ` when the residue is zero.
"""
struct NotInvertible <: Exception
    g
end

function Base.showerror(io::IO, e::NotInvertible)
    print(io, "element is not invertible modulo ψ; gcd = ", e.g)
end

"""
    inv_mod(a, ψ)

Return the inverse of `a` in `Fq[X]/(ψ)`. Throw `NotInvertible` if `a`
is not a unit. This is distinct from Julia's integer `invmod`.
"""
function inv_mod(a, ψ)
    a = mod(a, ψ)
    iszero(a) && throw(NotInvertible(ψ))
    g, s, _ = gcdx(a, ψ)
    degree(g) > 0 && throw(NotInvertible(g))
    return mod(s * inv(coeff(g, 0)), ψ)
end

"""
    comparepoints(P, Q, ψ) -> :equal | :opposite | :different

Compare points represented by `(a, b)` as `(a(X), Y*b(X))`, with `Y^2=f`.
`nothing` represents infinity. For odd-index division polynomials on a
nonsingular curve in characteristic different from the index, `gcd(f, ψ)=1`:
an odd-order torsion point cannot be a nontrivial point of order two.
Thus `f` is a unit modulo `ψ` and every factor of `ψ`. If the comparison
differs on only part of the quotient, throw `NotInvertible` with a gcd.
"""
function comparepoints(P, Q, ψ)
    (P === nothing || Q === nothing) &&
        return P === Q ? :equal : :different

    xd = mod(P[1] - Q[1], ψ)
    if !iszero(xd)
        g = gcd(xd, ψ)
        degree(g) > 0 && throw(NotInvertible(g))
        return :different
    end

    d = mod(P[2] - Q[2], ψ)
    iszero(d) && return :equal
    s = mod(P[2] + Q[2], ψ)
    iszero(s) && return :opposite

    # Equal x gives (b1-b2)*(b1+b2) = 0 because f is a unit.
    g = gcd(d, ψ)
    degree(g) == 0 && (g = gcd(s, ψ))
    degree(g) > 0 ||
        error("inconsistent point comparison: neither y-difference has a proper gcd")
    throw(NotInvertible(g))
end

"""
    ecdouble(P, A, f, ψ)

Double a point on the curve over `Fq[X]/(ψ)`.
"""
function ecdouble(P, A, f, ψ)
    P === nothing && return nothing
    x, b = P
    iszero(mod(b, ψ)) && return nothing
    m = mod((3x^2 + A) * inv_mod(2f*b, ψ), ψ)
    x3 = mod(f*m^2 - 2x, ψ)
    return (x3, mod(m*(x - x3) - b, ψ))
end

"""
    ecadd(P, Q, A, f, ψ)

Add two points on the curve over `Fq[X]/(ψ)`.
"""
function ecadd(P, Q, A, f, ψ)
    P === nothing && return Q
    Q === nothing && return P
    x1, b1 = P
    x2, b2 = Q

    if iszero(mod(x2 - x1, ψ))
        rel = comparepoints(P, Q, ψ)
        rel === :opposite && return nothing
        rel === :equal && return ecdouble(P, A, f, ψ)
        error("points have equal x-coordinates but are neither equal nor opposite")
    end

    m = mod((b2 - b1) * inv_mod(x2 - x1, ψ), ψ)
    x3 = mod(f*m^2 - x1 - x2, ψ)
    return (x3, mod(m*(x1 - x3) - b1, ψ))
end

"""
    ecneg(P) -> -P

Negate a point. Infinity is returned unchanged.
"""
ecneg(P) = P === nothing ? nothing : (P[1], -P[2])

"""
    scalarmul(P, k, A, f, ψ)

Compute `[k]P` by iterative double-and-add, reducing modulo `ψ`.
"""
function scalarmul(P, k, A, f, ψ)
    k < 0 && return scalarmul(ecneg(P), -k, A, f, ψ)
    (P === nothing || iszero(k)) && return nothing

    R = nothing
    while k > 0
        isodd(k) && (R = ecadd(R, P, A, f, ψ))
        k >>= 1
        k > 0 && (P = ecdouble(P, A, f, ψ))
    end
    return R
end

"""
    frobeniuspower(X, e, f, ψ)

Compute the image of `(X,Y)` under exponentiation by odd `e` in
`Fq[X,Y]/(Y^2-f, ψ)`, represented as `(X^e, f^((e-1)/2))`.
"""
function frobeniuspower(X, e, f, ψ)
    isodd(e) || throw(ArgumentError("Frobenius exponent must be odd"))
    return (powermod(X, e, ψ), powermod(f, (e - 1) ÷ 2, ψ))
end

"""
    tracemodℓ(ℓ, q, A, X, f, ψ) -> (t mod ℓ, case)

Compute the trace of Frobenius modulo the odd prime `ℓ`. The second return
value identifies the branch: `:pi2_is_minus_q`, `:nonresidue`,
`:residue`, or `:generic`.
"""
function tracemodℓ(ℓ, q, A, X, f, ψ)
    P = (X, one(parent(X)))
    qbar = Int(mod(q, ℓ))
    π1 = frobeniuspower(X, q, f, ψ)
    π2 = frobeniuspower(X, q^2, f, ψ)
    Q = scalarmul(P, qbar, A, f, ψ)

    rel = comparepoints(π2, Q, ψ)
    rel === :opposite && return (0, :pi2_is_minus_q)

    if rel === :equal
        powermod(qbar, (ℓ - 1) ÷ 2, ℓ) == 1 ||
            return (0, :nonresidue)
        w = sqrtmod(qbar, ℓ)
        W = scalarmul(P, w, A, f, ψ)
        r = comparepoints(π1, W, ψ)
        r === :equal && return (2w, :residue)
        r === :opposite && return (-2w, :residue)
        error("inconsistent Frobenius relation in residue case modulo $ℓ")
    end

    # π² + [qbar] = [τ]π: find τ in 1:(ℓ-1)/2, allowing either sign.
    S = ecadd(π2, Q, A, f, ψ)
    τπ = π1
    for τ in 1:((ℓ-1)÷2)
        r = comparepoints(S, τπ, ψ)
        r === :equal && return (τ, :generic)
        r === :opposite && return (-τ, :generic)
        τπ = ecadd(τπ, π1, A, f, ψ)
    end
    error("no trace candidate found modulo $ℓ")
end

"""
    divisionpolynomial(ℓ, A, B, R)

Return the odd-index division polynomial `ψ_ℓ` in `R = Fq[X]`.
Even-index terms are stored divided by `Y`, using `Y^2=f`.
"""
function divisionpolynomial(ℓ, A, B, R)
    isodd(ℓ) && ℓ >= 3 ||
        throw(ArgumentError("ℓ must be an odd integer ≥ 3"))
    F = base_ring(R)
    ℓ != characteristic(F) ||
        throw(ArgumentError("division-polynomial index must differ from the field characteristic"))

    X = gen(R)
    a, b = F(A), F(B)
    f2 = (X^3 + a*X + b)^2
    half = inv(F(2))
    g = Vector{elem_type(R)}(undef, max(ℓ, 4) + 1)
    g[1] = zero(R)
    g[2] = one(R)
    g[3] = R(2)
    g[4] = 3X^4 + 6a*X^2 + 12b*X - a^2
    g[5] = 4*(X^6 + 5a*X^4 + 20b*X^3 - 5a^2*X^2 -
        4a*b*X - 8b^2 - a^3)
    G(n) = g[n+1]

    for n in 5:ℓ
        m = n ÷ 2
        if isodd(n)
            t1 = G(m+2)*G(m)^3
            t2 = G(m-1)*G(m+1)^3
            g[n+1] = iseven(m) ? f2*t1 - t2 : t1 - f2*t2
        else
            g[n+1] = G(m) *
                (G(m+2)*G(m-1)^2 - G(m-2)*G(m+1)^2) * half
        end
    end
    return G(ℓ)
end

"""
    sqrtmod(a, ℓ) -> w

Return a square root of `a` modulo the small odd prime `ℓ` by trial division.
"""
function sqrtmod(a, ℓ)
    for w in 1:(ℓ-1)
        w*w % ℓ == a % ℓ && return w
    end
    error("$a is not a square modulo $ℓ")
end


end # module SchoofPointCounts
