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

"""
    ExactPermutationTest(data::Vector, f::Function)

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(data)` is ??? across all groups (typically more than 2 groups).  All
possible permutations are sampled.

The function `f` should reduce on the indices that comprise the groups in
the data using the statistics of interest.

# Examples
Say you want to test the difference in means between 3 groups (1,1), (2,2), 
and (9,9,9,9), you'll first need to concatenate your data into a vector.
Then you'll need to create a function that accumalativly subtracts the 
mean of the 1:2, 3:4, and 5:8 values of a given vector. 

```julia-repl

julia> indices = [1:2, 3:4, 5:8] # indices to the vector
3-element Array{UnitRange{Int64},1}:
 1:2
 3:4
 5:8

julia> data = [1,1,2,2,9,9,9,9] # the data vector
8-element Array{Int64,1}:
 1
 1
 2
 2
 9
 9
 9
 9

julia> f(x) = mapreduce(i -> mean(view(x, i)), -, indices) # calculate the cumalative
# difference between the means of the 1:2, 3:4, and 5:8 elements of vector `x`
f (generic function with 1 method)

julia> ExactPermutationTest(data, f)
Permutation Test
----------------
Population details:
    parameter of interest:   not implemented yet
    value under h_0:         NaN
    point estimate:          NaN

Test summary:
    outcome with 95% confidence: fail to reject h_0
    p-value:                     0.2024

Details:
    observation: -10.0
    samples: [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0  …  -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5]
```
"""
function ExactPermutationTest(data::AbstractVector{T}, 
                              f::Function) where {T<:Real}
    P = permutations(data)
    PermutationTest(f(data), f.(P))
end

"""
    ApproximatePermutationTest(data::Vector, f::Function, n::Int)

Perform a permutation test (a.k.a. randomization test) of the null hypothesis
that `f(data)` is ??? across all groups (typically more than 2 groups).
`n` of the `factorial(length(x)+length(y))` permutations are sampled at
random.

The function `f` should reduce on the indices that comprise the groups in
the data using the statistics of interest.

# Examples
Say you want to test the difference in means between 3 groups (1,1), (2,2), 
and (9,9,9,9), you'll first need to concatenate your data into a vector.
Then you'll need to create a function that accumalativly subtracts the 
mean of the 1:2, 3:4, and 5:8 values of a given vector. 

```julia-repl
julia> indices = [1:20, 21:40, 41:80]
3-element Array{UnitRange{Int64},1}:
 1:20
 21:40
 41:80

julia> data = [fill(1,20); fill(2, 20); fill(9, 40)];

julia> f(x) = mapreduce(i -> mean(view(x, i)), -, indices) # calculate the cumalative
# difference between the means of the 1:2, 3:4, and 5:8 elements of vector `x`

f (generic function with 1 method)

julia> ApproximatePermutationTest(data, f, 10^5)
Permutation Test
----------------
Population details:
    parameter of interest:   not implemented yet
    value under h_0:         NaN
    point estimate:          NaN

Test summary:
    outcome with 95% confidence: reject h_0
    p-value:                     <1e-4

Details:
    observation: -10.0
    samples: [-4.15, -5.85, -5.975, -5.949999999999999, -2.0749999999999993, -4.6, -4.6499999999999995, -6.05, -7.449999999999999, -5.0  …  -6.475, -4.45, -8.575, -5.6000000000000005, -4.075, -5.425, -7.5, -2.825, -5.175, -6.625]
```
"""
function ApproximatePermutationTest(data::AbstractVector{T}, f::Function,
                                    n::Int) where {T<:Real}
    observation = f(data)
    samples = [(shuffle!(data); f(data)) for _ in 1:n]
    PermutationTest(observation, samples)
end
