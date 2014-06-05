# kruskal_wallis.jl
# Kruskal-Wallis rank sum test
#
# Copyright (C) 2013   Simon Kornblith
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

# TODO(cs): Add post-hoc tests

immutable KruskalWallisTest <: HypothesisTest
    g::Int             # number of groups
    n::Int             # number of observations
    n_i::Vector{Int}   # number of observation in each group
    df::Int            # degrees of freedom
    R_i::Vector{Int}   # rank sums
    H::Float64         # test statistic: chi-square statistic
    C::Float64         # adjustment for ties
end

function KruskalWallisTest{T<:Real}(groups::AbstractVector{T}...) 
    g = length(groups)
    n_i = [length(g) for g in groups]
    n = sum(n_i)
    df = g - 1
    obs = [groups...]
    ranks = sortperm(sortperm(obs))
    
    # adjustment for ties

    # create dictionary from observation to ranks
    obs2ranks=Dict{Int, Vector{Int}}()
    for i=1:length(obs)
        old = get(obs2ranks, obs[i], [])
        obs2ranks[obs[i]] = [old, ranks[i]]
    end
    # map each observation to average rank (if ties)
    correction = [(k=>mean(v)) for (k,v) in zip(keys(obs2ranks), values(obs2ranks))]
    ranks = [correction[o] for o in obs]

    T_i = [length(d)^3-length(d) for d in values(obs2ranks)]
    C = 1-sum(T_i)/(n^3 - n)

    # compute rank sums and test statistic
    R_i = [sum(R) for R in split(ranks, n_i[1:end-1])]
    H = 12 * sum(R_i.^2./n_i) / (n * (n + 1)) - 3 * (n + 1) 
    H /= C

    KruskalWallisTest(g, n, n_i, df, R_i, H, C)
end
    
testname(::KruskalWallisTest) = "Kruskal-Wallis rank sum test"

pvalue(x::KruskalWallisTest) = pvalue(Chisq(x.df), x.H; tail=:right)
    

## helpers

# Return a tuple of vectors by splitting the given vector xs in multiple 
# vectors each with ns[i] components
function split(xs, ns::Vector{Int})
    if length(xs) == 0 || length(ns) == 0
        {xs}
    else
        rest = split(xs[ns[1]+1:end], ns[2:end])
        {xs[1:ns[1]], rest...}
    end
end
