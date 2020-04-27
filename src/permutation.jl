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
    ApproximatePermutationTest(x::Vector, y::Vector, f::Function, n::Int)

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(x)` is equal to `f(y)`.  `n` of the `factorial(length(x)+length(y))`
permutations are sampled at random.
"""
function ApproximatePermutationTest(x::AbstractVector{R}, y::AbstractVector{S},
                                    f::Function, n::Int) where {R<:Real,S<:Real}
    xy, rx, ry = ptstats(x,y)
    samples = [(shuffle!(xy); f(view(xy,rx)) - f(view(xy,ry))) for i = 1:n]
    PermutationTest(f(x) - f(y), samples)
end

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

function buildind(xs)
    i1 = 1
    y = Vector{UnitRange{Int64}}(undef, length(xs))
    for (i,x) in enumerate(xs)
        i2 = i1 + length(x) - 1
        y[i] = i1:i2
        i1 = i2 + 1
    end
    return y
end

function ExactPermutationTest(data::AbstractVector{T}, 
                              f::Function) where {T<:AbstractVector}
    xy = vcat(data...)
    r = buildind(data)
    P = permutations(xy)
    samples = [mapreduce(i -> f(view(p, i)), -, r) for p in P]
    PermutationTest(mapreduce(i -> f(view(xy, i)), -, r), samples)
end

function ApproximatePermutationTest(data::AbstractVector{T}, f::Function,
                                    n::Int) where {T<:AbstractVector}
    xy = vcat(data...)
    r = buildind(data)
    observation = mapreduce(i -> f(view(xy, i)), -, r)
    samples = [(shuffle!(xy); mapreduce(i -> f(view(xy, i)), -, r)) for _ in 1:n]
    PermutationTest(observation, samples)
end
