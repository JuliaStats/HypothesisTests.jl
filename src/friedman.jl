# Friedman.jl
# Friedman two-way ANOVA by ranks
#
# Copyright (C) 2012   Simon Kornblith
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

export FriedmanTest

struct FriedmanTest <: HypothesisTest
    n::Int                    # number of observations
    k::Int                    # number of treatments
    df::Int                   # degrees of freedom
    rank_sums::Vector{<:Real} # ranks sums vector
    chi_sq::Real              # chi-squared statistic
end

"""
    FriedmanTest(groups::AbstractVector{<:Real}...)

Perform the Friedman two-way ANOVA by ranks, a rank sum test to test the difference of ``k``
treatments across ``N`` repeated tests. This is a non-parametric test similar to the
Kruskall-Wallis one-way ANOVA by ranks. It is a special case of the Durbin test.

The p-value is computed using a ``χ^2`` statistic, as follows:
```math
    χ^2 & = \\frac{12}{Nk(N+1)} \\sum_{j = 1}^{k}R_j^2 - 3N(k+1)
```
where ``N`` is the number of tests, ``k`` is the number of treatments, and ``R`` is the rank
sum vector.

Implements: [`pvalue`](@ref)

# References

  * Daniel, W.W., Friedman two-way analysis of variance by ranks. Applied Nonparametric Statistics
    (2nd ed.). Boston: PWS-Kent. pp. 262–74, 1990.

# External links

  * [Friedman test on Wikipedia
    ](https://en.wikipedia.org/wiki/Friedman_test)
"""
function FriedmanTest(groups::AbstractVector{T}...) where {T<:Real}
    x = mapreduce(permutedims, vcat, groups)
    n = size(x, 2)
    k = size(x, 1)
    df = k - 1
    rank_sums = sum.(convert_rank(x) |> eachcol)
    squared_rank_sums_sum = sum(rank_sum^2 for rank_sum in rank_sums)
    chi_sq = (12 / (n*k*(k+1))) * squared_rank_sums_sum - 3*n*(k+1)
    FriedmanTest(n, k, df, rank_sums, chi_sq)
end

testname(::FriedmanTest) = "Friedman two-way ANOVA by ranks"
population_param_of_interest(x::FriedmanTest) = ("location parameter", "all equal", NaN) # parameter of interest: name, value under h0, point estimate
default_tail(test::FriedmanTest) = :both

function show_params(io::IO, x::FriedmanTest, ident)
    println(io, ident, "number of observations: ", x.n)
    println(io, ident, "number of treatments:   ", x.k)
    println(io, ident, "degrees of freedom:     ", x.df)
    println(io, ident, "rank sums vector:       ", x.rank_sums)
    println(io, ident, "χ² statistic:           ", x.chi_sq)
end


StatsAPI.pvalue(x::FriedmanTest; tail=:both) = pvalue(Chisq(x.df), x.chi_sq, tail=tail)

#helper functions

#find the average index of an item in a vector
function avg_index(item::T, x::AbstractVector{T}) where {T<:Real}
    equal = item .== x
    sum([i for (i, value) in enumerate(equal) if value]) / sum(equal) # get average rank (for ties)
end

#find the ranks for a single test/row
function row_rank(row::AbstractVector{T}) where {T<:Real}
    sorted_row = sort(row)
    map(x -> avg_index(x, sorted_row), row)
end

#convert an input matrix into ranks
function convert_rank(x::AbstractMatrix{T}) where {T<:Real}
    x = map(row_rank, eachrow(x))
    mapreduce(permutedims, vcat, x)
end