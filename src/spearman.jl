# spearman.jl
# Spearman's rank correlation test in Julia
#
# Copyright (C) 2016   Diego Javier Zea
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

export CorrelationTest, SpearmanCorrelationTest

abstract CorrelationTest <: HypothesisTest

"Sum squared difference of ranks (midranks for ties)"
spearman_S(xrank, yrank) = sumabs2(xrank .- yrank)

immutable SpearmanCorrelationTest <: CorrelationTest
    # Tied ranking for x and y
    xrank::Vector{Float64}
    yrank::Vector{Float64}
    # Adjustment for ties
    xtiesadj::Float64
    ytiesadj::Float64
    # Sum squared difference of ranks
    S::Float64
    # Number of points
    n::Int
    # Spearman's ρ
    ρ::Float64

    function SpearmanCorrelationTest(x, y)

        n = length(x)
        (n != length(y)) && throw(ErrorException("x and y must have the same length"))

        xrank, xtiesadj = HypothesisTests.tiedrank_adj(x)
        yrank, ytiesadj = HypothesisTests.tiedrank_adj(y)

        S = spearman_S(xrank, yrank)

        ρ = corspearman(x, y)

        new(xrank, yrank, xtiesadj, ytiesadj, S, n, ρ)
    end

end

testname(::SpearmanCorrelationTest) = "Spearman's rank correlation test"

# parameter of interest: name, value under h0, point estimate
population_param_of_interest(x::SpearmanCorrelationTest) = ("Spearman's ρ", 0.0, x.ρ)

function show_params(io::IO, x::SpearmanCorrelationTest, ident)
    println(io, ident, "Number of points:                    ", x.n)
    println(io, ident, "Spearman's ρ:                        ", x.ρ)
    println(io, ident, "S (Sum squared difference of ranks): ", x.S)
    println(io, ident, "adjustment for ties in x:            ", x.xtiesadj)
    println(io, ident, "adjustment for ties in y:            ", x.ytiesadj)
end

function P_from_null_S_values(S_null, x::SpearmanCorrelationTest, tail)
    S_null_mean = mean(S_null)
    # S is approximately normally distributed
    # S and ρ are inversely proportional
    S_null[:] = S_null .- S_null_mean # center
    S_centered = x.S - S_null_mean    # center
    if tail == :both
        modS = abs(S_centered)
        mean(S_null .<= -modS) + mean(S_null .>= modS)
    elseif tail == :right
        mean(S_null .<= S_centered)
    elseif tail == :left
        mean(S_null .>= S_centered)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

function spearman_P_exact(x::SpearmanCorrelationTest, tail)
    S_null = Float64[ spearman_S(perm, x.yrank) for perm in permutations(x.xrank) ]
    P_from_null_S_values(S_null, x, tail)
end

function spearman_P_sampling(x::SpearmanCorrelationTest, tail)
    # 360000 samples gives an se(P) < 0.0005 for P < 0.1
    X = copy(x.xrank)
    S_null = Float64[ spearman_S(shuffle!(X), x.yrank) for sample in 1:360000 ]
    P_from_null_S_values(S_null, x, tail)
end

# Use estimated mean and std for the S null distribution as in:
#
# Press WH, Teukolsky SA, Vetterling WT, Flannery BP.
# Numerical recipes in C.
# Cambridge: Cambridge university press; 1996.
function spearman_P_estimated(x::SpearmanCorrelationTest, tail)
    N = float(x.n)
    a = (N^3 - N)
    # Numerical Recipes (14.6.6)
    S_null_mean = (a/6.) - (x.xtiesadj/12.) - (x.ytiesadj/12.)
    # Numerical Recipes (14.6.7)
    S_null_std  = sqrt( ((N - 1.)/36.) * (1. - (x.xtiesadj/a) ) * (1. - (x.ytiesadj/a) ) ) * N * (N + 1.)
    zscore = (x.S - S_null_mean)/S_null_std
    # S is approximately normally distributed
    # S and ρ are inversely proportional
    if tail == :both
        cdf(Normal(), -abs(zscore)) + ccdf(Normal(), abs(zscore))
    elseif tail == :right
        cdf(Normal(),  zscore)
    elseif tail == :left
        ccdf(Normal(), zscore)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# Using T test for n > 10 as in:
#
# McDonald JH.
# Handbook of biological statistics.
# Baltimore, MD: Sparky House Publishing; 2009 Aug.
function spearman_P_ttest(x::SpearmanCorrelationTest, tail)
    ρ2 = x.ρ^2
    df = x.n-2
    t = sqrt((df*ρ2)/(1-ρ2))
    if tail == :both
        cdf(TDist(df), -t) + ccdf(Normal(), t)
    elseif tail == :right
        ccdf(TDist(df), t)
    elseif tail == :left
        cdf(TDist(df),  t)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

function pvalue(x::SpearmanCorrelationTest; tail=:both, method=:estimated)
    if x.n <= 10
        # Exact P value using permutations
        return( spearman_P_exact(x, tail) )
    end
    if method == :sampling
        return( spearman_P_sampling(x, tail) )
    elseif method == :exact
        return( spearman_P_exact(x, tail) )
    elseif method == :estimated
        return( spearman_P_estimated(x, tail) )
    elseif method == :ttest
        return( spearman_P_ttest(x, tail) )
    else
        throw(ArgumentError("method=$(method) is invalid"))
    end
end

