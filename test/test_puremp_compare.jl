""" This file tests a Rust Schoof implementation against the Julia code's implementation of Schoof. """

using Test
using Random
using Primes

# Compares SchoofPointCounts.schoof against the Rust crate `puremp`
# Rust must be installed and reachable to build the testing harness below.
# (EllipticCurve::schoof_point_count and ::sea_point_count) through the small CLI in
# puremp_compare/. puremp only supports prime fields GF(p), so only those are compared.

const MODULE_FILE = get(ENV, "SCHOOF_MODULE_FILE") do
    local_copy = joinpath(@__DIR__, "SchoofPointCounts.jl")
    isfile(local_copy) ? local_copy : joinpath(homedir(), "Downloads", "SchoofPointCounts.jl")
end
include(MODULE_FILE)
using .SchoofPointCounts

const RUST_DIR = joinpath(@__DIR__, "puremp_compare")
const RUST_EXE = joinpath(RUST_DIR, "target", "release",
                          Sys.iswindows() ? "puremp_compare.exe" : "puremp_compare")

if !isfile(RUST_EXE)
    @info "Building the puremp comparison binary"
    run(Cmd(`cargo build --release`; dir = RUST_DIR))
end

"""
Run the Rust tool on `(p, A, B)` triples; returns a Dict mapping each triple to
`:singular` or `(schoof, sea)` as BigInts. Coefficients are reduced mod p first.
"""
function puremp_counts(cases)
    input = join(("$p $(mod(A, p)) $(mod(B, p))" for (p, A, B) in cases), "\n") * "\n"
    out = open(`$RUST_EXE`, "r+") do io
        write(io, input)
        close(io.in)
        read(io, String)
    end
    results = Dict{Tuple{BigInt,BigInt,BigInt},Any}()
    for line in split(out, '\n'; keepempty = false)
        f = split(line)
        key = (parse(BigInt, f[1]), parse(BigInt, f[2]), parse(BigInt, f[3]))
        results[key] = f[4] == "singular" ? :singular : (parse(BigInt, f[4]), parse(BigInt, f[5]))
    end
    return results
end

rng = MersenneTwister(20261109)

fixed = [(5, 1, 1), (7, 2, 3), (11, 1, 6), (13, 1, 2), (97, 2, 3), (101, 1, 1),
         (1009, 2, 3), (10007, 5, 7), (10007, -3, 4), (65537, 3, 8), (100003, 2, 3),
         (1000003, 2, 3), (1000033, 1, 1)]

function random_cases(rng, ndigits, n)
    cases = Tuple{BigInt,BigInt,BigInt}[]
    while length(cases) < n
        p = nextprime(rand(rng, big(10)^(ndigits - 1):big(10)^ndigits - 1))
        p > 3 || continue
        push!(cases, (p, rand(rng, 0:p-1), rand(rng, 0:p-1)))
    end
    return cases
end

# small primes: every curve, so singular, special-case and ordinary curves are all hit
exhaustive = [(p, A, B) for p in (5, 7, 11, 13) for A in 0:p-1 for B in 0:p-1]

groups = [
    "fixed curves"            => fixed,
    "every curve over F_p, p ≤ 13" => exhaustive,
    "random 3-digit primes"   => random_cases(rng, 3, 25),
    "random 6-digit primes"   => random_cases(rng, 6, 15),
    "random 9-digit primes"   => random_cases(rng, 9, 6),
]

@testset "SchoofPointCounts vs puremp" begin
    for (name, cases) in groups
        cases = [(BigInt(p), BigInt(A), BigInt(B)) for (p, A, B) in cases]
        rust = puremp_counts(cases)
        @testset "$name" begin
            for (p, A, B) in cases
                expected = rust[(p, mod(A, p), mod(B, p))]
                if expected === :singular
                    @test_throws ArgumentError schoof(Int(mod(A, p)), Int(mod(B, p)), Int(p))
                else
                    got = schoof(Int(mod(A, p)), Int(mod(B, p)), Int(p))
                    @test got == expected[1] # puremp Schoof
                    @test got == expected[2] # puremp SEA
                end
            end
        end
    end

    @testset "invalid orders are rejected by the Julia module" begin
        @test_throws ArgumentError schoof(1, 1, 8)
        @test_throws ArgumentError schoof(1, 1, 9)
        @test_throws ArgumentError schoof(1, 1, 15)
    end
end
