"""
    This test file runs the code against precomputed known values. For further testing you may want to
    run the test_pureemp_compare,jl file, which requires Rust to ber installed.
"""

using Test
using Nemo
using SchoofPointCounts
const S = SchoofPointCounts # shorthand to reference nonexported functions

""" Independent direct point count over any finite field (no Schoof code involved). """
function bruteforce(A, B)
    F = parent(A)
    q = BigInt(order(F))
    p = Int(characteristic(F))
    d = round(Int, log(p, q))
    g = gen(F)
    n = 1 # point at infinity
    for coeffs in Iterators.product(ntuple(_ -> 0:p-1, d)...)
        x = sum(F(c) * g^(i - 1) for (i, c) in enumerate(coeffs))
        r = x^3 + A*x + B
        n += iszero(r) ? 1 : (r^((q - 1) ÷ 2) == 1 ? 2 : 0)
    end
    return n
end

primefield(p) = finite_field(p, 1, "a")[1]

@testset "Schoof point counting" begin
    @testset "direct comparison, prime fields" begin
        for (p, A, B) in [(5,1,1), (7,2,3), (11,1,6), (13,1,2), (97,2,3), (101,1,1),
                          (1009,2,3), (10007,-3,4), (65537,3,8)]
            F = primefield(p)
            @test S.schoof(A, B, p) == bruteforce(F(A), F(B))
        end
    end

    @testset "direct comparison, extension fields" begin
        for (p, k) in [(5, 2), (7, 2), (5, 3)]
            F, a = finite_field(p, k, "a")
            for (A, B) in [(a, F(1)), (a + 1, a^2), (F(1), a), (F(2), F(3))]
                iszero(4A^3 + 27B^2) && continue
                @test S.schoof(A, B) == bruteforce(A, B)
            end
        end
        @test S.schoof(1, 1, 25) == bruteforce(finite_field(5, 2, "a")[1](1), finite_field(5, 2, "a")[1](1))
        @test S.schoof(2, 3, 49) isa Integer
    end

    @testset "input validation" begin
        @test_throws ArgumentError S.schoof(1, 1, 8)    # characteristic 2
        @test_throws ArgumentError S.schoof(1, 1, 9)    # characteristic 3
        @test_throws ArgumentError S.schoof(1, 1, 27)
        @test_throws ArgumentError S.schoof(0, 0, 7)    # singular
        @test_throws ArgumentError S.schoof(-3, 2, 101) # x³-3x+2 = (x-1)²(x+2), singular
        @test_throws ArgumentError S.schoof(1, 1, 15)   # not a prime power
    end

    @testset "CRT bound and reconstruction" begin
        for q in (5, 7, 101, 10007, big(10)^30 + 57), p in (5, 7, 11)
            Sset, M = S.auxiliaryprimes(p, q)
            @test !(p in Sset)
            @test M^2 > 16 * big(q)
            @test M == 2 * prod(big.(Sset))
        end
        # every t in the Hasse interval is recovered, including the boundaries
        q = 101
        _, M = S.auxiliaryprimes(11, q)
        Sset, _ = S.auxiliaryprimes(11, q)
        for t in -2*isqrt(q)-1:2*isqrt(q)+1
            cong = [(2, mod(t, 2))]
            append!(cong, [(ℓ, mod(t, ℓ)) for ℓ in Sset])
            @test S.reconstruct_trace(cong) == t
        end
    end

    @testset "division polynomials" begin
        F = primefield(10007)
        R, X = polynomial_ring(F, "X")
        for (A, B) in [(2, 3), (5, 7), (-3, 4)]
            a, b = F(A), F(B)
            @test S.divisionpolynomial(5, a, b, R) == 5X^12 + 62a*X^10 + 380b*X^9 - 105a^2*X^8 +
                240a*b*X^7 + (-300a^3 - 240b^2)*X^6 - 696a^2*b*X^5 + (-125a^4 - 1920a*b^2)*X^4 +
                (-80a^3*b - 1600b^3)*X^3 + (-50a^5 - 240a^2*b^2)*X^2 +
                (-100a^4*b - 640a*b^3)*X + (a^6 - 32a^3*b^2 - 256b^4)
            for ℓ in (3, 5, 7, 11, 13)
                ψ = S.divisionpolynomial(ℓ, a, b, R)
                @test degree(ψ) == (ℓ^2 - 1) ÷ 2
                @test coeff(ψ, degree(ψ)) == ℓ
            end
        end
        # ψ_ℓ vanishes at the x-coordinate of every F_p-rational point of order dividing ℓ
        p = 103
        Fp = primefield(p)
        Rp, Xp = polynomial_ring(Fp, "X")
        A, B = 2, 3
        affine_add(P, Q) = begin
            P === nothing && return Q
            Q === nothing && return P
            x1, y1 = P; x2, y2 = Q
            x1 == x2 && (y1 + y2) % p == 0 && return nothing
            m = x1 == x2 ? (3x1^2 + A) * invmod(2y1, p) : (y2 - y1) * invmod(x2 - x1, p)
            m = mod(m, p)
            x3 = mod(m^2 - x1 - x2, p)
            (x3, mod(m * (x1 - x3) - y1, p))
        end
        pts = [(x, y) for x in 0:p-1 for y in 0:p-1 if mod(y^2 - x^3 - A*x - B, p) == 0]
        for ℓ in (3, 5, 7, 11, 13)
            ψ = S.divisionpolynomial(ℓ, Fp(A), Fp(B), Rp)
            for P in pts
                Q = nothing
                for _ in 1:ℓ
                    Q = affine_add(Q, P)
                end
                if Q === nothing
                    @test iszero(evaluate(ψ, Fp(P[1])))
                end
            end
        end
    end

    @testset "polynomial inversion and zero divisors" begin
        F = primefield(101)
        R, X = polynomial_ring(F, "X")
        ψ = X * (X + 1)
        err = try S.inv_mod(X, ψ); nothing catch e; e end
        @test err isa S.NotInvertible
        @test err.g == X
        u = X + 2
        @test mod(u * S.inv_mod(u, ψ), ψ) == 1
    end

    @testset "point arithmetic" begin
        F = primefield(101)
        R, X = polynomial_ring(F, "X")
        A, B = R(2), R(3)
        f = X^3 + 2X + 3
        ψ = S.divisionpolynomial(3, F(2), F(3), R)
        P = (X, one(R))
        @test S.ecadd(P, S.ecneg(P), A, f, ψ) === nothing
        @test S.scalarmul(P, 0, A, f, ψ) === nothing
        @test S.comparepoints(S.scalarmul(P, 2, A, f, ψ), P, ψ) == :opposite
        @test S.scalarmul(P, 3, A, f, ψ) === nothing
        @test S.scalarmul(P, -2, A, f, ψ) == S.ecneg(S.scalarmul(P, 2, A, f, ψ))
    end

    @testset "trace_mod_l cases" begin
        seen = Set{Symbol}()
        for p in (11, 13, 17, 19, 23), A in 0:p-1, B in 0:p-1
            F = primefield(p)
            iszero(4F(A)^3 + 27F(B)^2) && continue
            t = p + 1 - bruteforce(F(A), F(B))
            R, X = polynomial_ring(F, "X")
            f = X^3 + F(A)*X + F(B)
            for ℓ in (3, 5, 7)
                ℓ == p && continue
                ψ = S.divisionpolynomial(ℓ, F(A), F(B), R)
                result = nothing
                while result === nothing
                    try
                        result = S.tracemodℓ(ℓ, p, R(A), X, f, ψ)
                    catch e
                        e isa S.NotInvertible || rethrow()
                        ψ = e.g * inv(leading_coefficient(e.g))
                    end
                end
                tl, kind = result
                @test mod(tl, ℓ) == mod(t, ℓ)
                push!(seen, kind)
            end
        end
        # π²P = +[q̄]P forces t² ≡ 4q̄, so :nonresidue is only a defensive guard
        @test seen == Set([:pi2_is_minus_q, :residue, :generic])
    end

    @testset "sqrtmod" begin
        for ℓ in (3, 5, 7, 11, 13), a in 1:ℓ-1
            if powermod(a, (ℓ - 1) ÷ 2, ℓ) == 1
                @test S.sqrtmod(a, ℓ)^2 % ℓ == a
            else
                @test_throws ErrorException S.sqrtmod(a, ℓ)
            end
        end
    end
end

true
