export PermutationTest

function ptstats(x,y)
    xy = vcat(x,y)
    rx = 1:length(x)
    ry = (length(xy) - length(y) + 1):length(xy)
    return (xy, rx, ry)
end

"""
    PermutationTest(f, x::AbstractVector{<:Real}, y::AbstractVector{<:Real} [, n::Integer])

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(x)` is equal to `f(y)`. If `n` is provided, an approximate permutation test is
performed, where `n` of the `factorial(length(x) + length(y))` permutations are sampled at
random. Otherwise, an exact permutation test is performed, where all possible permutations
are sampled.
"""
struct PermutationTest{T<:Real} <: HypothesisTest
    observation::T
    samples::Vector{T}
end

function PermutationTest(
    f::Function,
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
)
    xy, rx, ry = ptstats(x, y)
    P = permutations(xy)
    samples = [f(view(p, rx)) - f(view(p, ry)) for p in P]
    PermutationTest(f(x) - f(y), samples)
end

function PermutationTest(
    f::Function,
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    n::Integer,
)
    xy, rx, ry = ptstats(x, y)
    samples = [(shuffle!(xy); f(view(xy, rx)) - f(view(xy, ry))) for i = 1:n]
    PermutationTest(f(x) - f(y), samples)
end

function pvalue(pt::PermutationTest; tail=:both)
    c = if tail == :both
        count(abs(pt.observation) ≤ abs(x) for x in pt.samples)
    elseif tail == :left
        count(pt.observation ≥ x for x in pt.samples)
    elseif tail == :right
        count(pt.observation ≤ x for x in pt.samples)
    end
    return c / length(pt.samples)
end

testname(::PermutationTest) = "Permutation Test"

function show_params(io::IO, pt::PermutationTest, ident)
    println(io, ident, "observation: ", pt.observation)
    print(io, ident, "samples: ")
    show(io, pt.samples)
    println(io)
end
