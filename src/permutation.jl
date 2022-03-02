export ExactPermutationTest, ApproximatePermutationTest

function ptstats(x,y)
    xy = vcat(x,y)
    rx = 1:length(x)
    ry = (length(xy)-length(y)+1):length(xy)
    (xy, rx, ry)
end

struct PermutationTest{T<:Real} <: HypothesisTest
    observation::T
    samples::Vector{T}
end

"""
    ExactPermutationTest(x::Vector, y::Vector, f::Function)

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(x)` is equal to `f(y)`.  All possible permutations are sampled.
"""
function ExactPermutationTest(x::AbstractVector{R}, y::AbstractVector{S},
                              f::Function) where {R<:Real,S<:Real}
    xy, rx, ry = ptstats(x,y)
    P = permutations(xy)
    samples = [f(view(p,rx)) - f(view(p,ry)) for p in P]
    PermutationTest(f(x) - f(y), samples)
end

"""
    ApproximatePermutationTest([rng::AbstractRNG,] x::Vector, y::Vector, f::Function, n::Int)

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(x)` is equal to `f(y)`.  `n` of the `factorial(length(x)+length(y))`
permutations are sampled at random. A random number generator can optionally
be passed as the first argument. The default generator is `Random.default_rng()`.
"""
function ApproximatePermutationTest(rng::AbstractRNG, x::AbstractVector{R}, y::AbstractVector{S},
                                    f::Function, n::Int) where {R<:Real,S<:Real}
    xy, rx, ry = ptstats(x,y)
    samples = [(shuffle!(rng, xy); f(view(xy,rx)) - f(view(xy,ry))) for i = 1:n]
    PermutationTest(f(x) - f(y), samples)
end
ApproximatePermutationTest(x::AbstractVector{R}, y::AbstractVector{S},
                           f::Function, n::Int) where {R<:Real,S<:Real} =
    ApproximatePermutationTest(Random.default_rng(), x, y, f, n)

function pvalue(apt::PermutationTest; tail=:both)
    if tail == :both
        count = sum(abs(apt.observation) <= abs(x) for x in apt.samples)
    elseif tail == :left
        count = sum(apt.observation >= x for x in apt.samples)
    elseif tail == :right
        count = sum(apt.observation <= x for x in apt.samples)
    end
    return count / length(apt.samples)
end

testname(::PermutationTest) = "Permutation Test"

function show_params(io::IO, apt::PermutationTest, ident)
    println(io, ident, "observation: ", apt.observation)
    print(io, ident, "samples: ")
    show(io, apt.samples)
    println(io)
end
