# kruskal_wallis.jl
# Kruskal-Wallis rank sum test
#
# Copyright (C) 2014   Christoph Sawade
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

export KruskalWallisTest

struct KruskalWallisTest <: HypothesisTest
    n_i::Vector{Int}         # number of observations in each group
    df::Int                  # degrees of freedom
    R_i::Vector{Float64}     # rank sums
    H::Float64               # test statistic: chi-square statistic
    tie_adjustment::Float64  # adjustment for ties
end

"""
    KruskalWallisTest(groups::AbstractVector{<:Real}...)

Perform Kruskal-Wallis rank sum test of the null hypothesis that the `groups`
``\\mathcal{G}`` come from the same distribution against the alternative hypothesis that
that at least one group stochastically dominates one other group.

The Kruskal-Wallis test is an extension of the Mann-Whitney U test to more than two groups.

The p-value is computed using a ``χ^2`` approximation to the distribution of the test
statistic ``H_c=\\frac{H}{C}``:
```math
    \\begin{align*}
    H & = \\frac{12}{n(n+1)} \\sum_{g ∈ \\mathcal{G}} \\frac{R_g^2}{n_g} - 3(n+1)\\\\
    C & = 1-\\frac{1}{n^3-n}\\sum_{t ∈ \\mathcal{T}} (t^3-t),
    \\end{align*}
```
where ``\\mathcal{T}`` is the set of the counts of tied values at each tied position, ``n``
is the total number of observations across all groups, and ``n_g`` and ``R_g`` are the number of
observations and the rank sum in group ``g``, respectively. See references for further
details.

Implements: [`pvalue`](@ref)

# References

  * Meyer, J.P, Seaman, M.A., Expanded tables of critical values for the Kruskal-Wallis
    H statistic. Paper presented at the annual meeting of the American Educational Research
    Association, San Francisco, April 2006.

# External links

  * [Kruskal-Wallis test on Wikipedia
    ](https://en.wikipedia.org/wiki/Kruskal-Wallis_one-way_analysis_of_variance)
"""
function KruskalWallisTest(groups::AbstractVector{T}...) where T<:Real
    (H, R_i, tieadj, n_i) = kwstats(groups...)
    if length(groups)<=3 && any(n_i .< 6)
        @warn("This test is only asymptotically correct and might be inaccurate for the given group size")
    end
    df = length(groups) - 1
    KruskalWallisTest(n_i, df, R_i, H, tieadj)
end

testname(::KruskalWallisTest) = "Kruskal-Wallis rank sum test (chi-square approximation)"
population_param_of_interest(x::KruskalWallisTest) = ("Location parameters", "all equal", NaN) # parameter of interest: name, value under h0, point estimate
default_tail(test::KruskalWallisTest) = :right

function show_params(io::IO, x::KruskalWallisTest, ident)
    print(io, ident, "number of observation in each group: ")
    show(io, x.n_i)
    println(io)
    println(io, ident, "χ²-statistic:                        ", x.H)
    print(io, ident, "rank sums:                           ")
    show(io, x.R_i)
    println(io)
    println(io, ident, "degrees of freedom:                  ", x.df)
    println(io, ident, "adjustment for ties:                 ", x.tie_adjustment)
end

pvalue(x::KruskalWallisTest) = pvalue(Chisq(x.df), x.H; tail=:right)


## helper

# Get H, rank sums, and tie adjustment for Kruskal-Wallis test
function kwstats(groups::AbstractVector{T}...) where T<:Real
    n_i = [length(g) for g in groups]
    n = sum(n_i)

    # get ranks and adjustment for ties
    (ranks, tieadj) = tiedrank_adj([groups...;])
    C = 1-tieadj/(n^3 - n)

    # compute rank sums
    R_i = Vector{Float64}(undef, length(groups))
    n_end = 0
    for i in 1:length(groups)
        R_i[i] = sum(ranks[n_end+1:n_end+n_i[i]])
        n_end += n_i[i]
    end

    # compute test statistic and correct for ties
    H = 12 * sum(R_i .^ 2 ./ n_i) / (n * (n + 1)) - 3 * (n + 1)
    H /= C

    (H, R_i, C, n_i)
end
