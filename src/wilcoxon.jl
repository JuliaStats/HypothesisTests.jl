# Wilcoxon.jl
# Wilcoxon rank sum (Mann-Whitney U) and signed rank tests in Julia
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

export MannWhitneyUTest, ExactMannWhitneyUTest, ApproximateMannWhitneyUTest,
    SignedRankTest, ExactSignedRankTest, ApproximateSignedRankTest

# RMATH WRAPPERS
macro rmath_deferred_free(base)
    libcall = symbol(string(base, "_free"))
    func = symbol(string(base, "_deferred_free"))
    quote
        let gc_tracking_obj = []
            global $func
            function $libcall(x::Vector{None})
                gc_tracking_obj = []
                ccall(($(string(libcall)),:libRmath), Void, ())
            end
            function $func()
                if !isa(gc_tracking_obj, Bool)
                    finalizer(gc_tracking_obj, $libcall)
                    gc_tracking_obj = false
                end
            end
        end
    end
end

@rmath_deferred_free(signrank)
function psignrank(q::Number, p1::Number, lower_tail::Bool,
                   log_p::Bool=false)
    signrank_deferred_free()
    ccall((:psignrank,:libRmath), Float64, (Float64,Float64,Int32,Int32), q, p1,
          lower_tail, log_p)
end
@rmath_deferred_free(wilcox)
function pwilcox(q::Number, p1::Number, p2::Number, lower_tail::Bool,
                 log_p::Bool=false)
    wilcox_deferred_free()
    ccall((:pwilcox,:libRmath), Float64, (Float64,Float64,Float64,Int32,Int32),
          q, p1, p2, lower_tail, log_p)
end

# TYPES

immutable ApproximateSignedRankTest <: HypothesisTest
    W::Float64
    tie_adjustment::Float64
    mu::Float64
    sigma::Float64
end

# COMMON TO ALL

# Tied rank from Base, modified to compute the adjustment for ties
function tiedrank_adj(v::AbstractArray)
    n     = length(v)
    place = sortperm(v)
    ord   = Array(Float64, n)
    tieadj = 0.0

    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[place[i]] == v[place[j + 1]]
            j += 1
        end

        if j > i
            t = j - i + 1
            m = sum(i:j) / t
            tieadj += t^3 - t
            for k = i:j
                ord[place[k]] = m
            end
        else
            ord[place[i]] = i
        end

        i = j + 1
    end

    (ord, tieadj)
end

## COMMON MANN-WHITNEY U

# Get U, ranks, and tie adjustment for Mann-Whitney U test
function mwustats{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T})
    nx = length(x)
    ny = length(y)
    if nx <= ny
        (ranks, tieadj) = tiedrank_adj([x, y])
        U = sum(ranks[1:nx]) - nx*(nx+1)/2
    else
        (ranks, tieadj) = tiedrank_adj([y, x])
        U = nx*ny - sum(ranks[1:ny]) + ny*(ny+1)/2
    end
    (U, ranks, tieadj, nx, ny)
end

# Automatic exact/normal selection
function MannWhitneyUTest{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T})
    (U, ranks, tieadj, nx, ny) = mwustats(x, y)
    if nx + ny <= 10 || (nx + ny <= 50 && tieadj == 0)
        ExactMannWhitneyUTest(U, ranks, tieadj, nx, ny)
    else
        ApproximateMannWhitneyUTest(U, ranks, tieadj, nx, ny)
    end
end

## EXACT MANN-WHITNEY U TEST

immutable ExactMannWhitneyUTest{T <: Real} <: HypothesisTest
    U::Float64
    ranks::Vector{T}
    tie_adjustment::Float64
    nx::Int
    ny::Int
end
ExactMannWhitneyUTest{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) =
    ExactMannWhitneyUTest(mwustats(x, y)...)

testname(::ExactMannWhitneyUTest) = "Exact Mann-Whitney U test"

# Enumerate all possible Mann-Whitney U results for a given vector,
# determining left-and right-tailed p values
function mwuenumerate(x::ExactMannWhitneyUTest)
    # Get the other U if inverted by mwu_stats
    n = min(x.nx, x.ny)
    if x.ny > x.nx
        U = x.nx*x.ny - U
    end
    le = 0
    gr = 0
    tot = 0
    k = n*(n+1)/2
    for comb in @task combinations(x.ranks, n)
        Up = sum(comb) - k
        tot += 1
        le += Up <= x.U
        gr += Up >= x.U
    end
    (le/tot, gr/tot)
end

pvalue(x::ExactMannWhitneyUTest; tail=:both) =
    if x.tie_adjustment == 0
        # Compute exact p-value using method from Rmath, which is fast but
        # cannot account for ties
        if tail == :both
            if x.U < x.nx * x.ny / 2
                2 * pwilcox(x.U, x.nx, x.ny, true)
            else
                2 * pwilcox(x.U - 1, x.nx, x.ny, false)
            end
        elseif tail == :left
            pwilcox(x.U, x.nx, x.ny, true)
        elseif tail == :right
            pwilcox(x.U - 1, x.nx, x.ny, false)
        else
            error("tail=$(tail) is invalid")
        end
    else
        # Compute exact p-value by enumerating possible ranks in the tied data
        if tail == :both
            min(1, 2 * min(mwuenumerate(x)))
        elseif tail == :left
            mwuenumerate(x)[1]
        elseif tail == :right
            mwuenumerate(x)[2]
        else
            error("tail=$(tail) is invalid")
        end
    end

## APPROXIMATE MANN-WHITNEY U TEST

immutable ApproximateMannWhitneyUTest <: HypothesisTest
    U::Float64
    tie_adjustment::Float64
    mu::Float64
    sigma::Float64
end
function ApproximateMannWhitneyUTest(U::Real, ::Vector, tie_adjustment::Real,
                                     nx::Int, ny::Int)
    mu = U - nx * ny / 2
    sigma = sqrt((nx * ny * (nx + ny + 1 - tie_adjustment /
        ((nx + ny) * (nx + ny - 1)))) / 12)
    ApproximateMannWhitneyUTest(U, tie_adjustment, mu, sigma)
end
ApproximateMannWhitneyUTest{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) =
    ApproximateMannWhitneyUTest(mwustats(x, y)...)

testname(::Type{ApproximateMannWhitneyUTest}) =
    "Approximate Mann-Whitney U test"

pvalue(x::Union(ApproximateMannWhitneyUTest, ApproximateSignedRankTest); tail=:both) =
    if x.mu == x.sigma == 0
        1
    else
        if tail == :both
            2 * ccdf(Normal(), abs(x.mu - 0.5 * sign(x.mu))/x.sigma)
        elseif tail == :left
            cdf(Normal(), (x.mu + 0.5)/x.sigma)
        elseif tail == :right
            ccdf(Normal(), (x.mu - 0.5)/x.sigma)
        end
    end

## COMMON SIGNED RANK

# Get W and absolute ranks for signed rank test
function signedrankstats{S <: Real}(x::Vector{S})
   nonzero_x = x[x .!= 0]
   (ranks, tieadj) = tiedrank_adj(abs(nonzero_x))
   W = 0.0
   for i = 1:length(nonzero_x)
           if nonzero_x[i] > 0
                   W += ranks[i]
           end
   end
   (W, ranks, tieadj)
end


# Automatic exact/normal selection
function SignedRankTest{T <: Real}(x::Vector{T})
    (W, ranks, tie_adjustment) = signedrankstats(x)
    n = length(ranks)
    if n <= 15 || (n <= 50 && tieadj == 0)
        ExactSignedRankTest(W, ranks, tie_adjustment)
    else
        ApproximateSignedRankTest(W, tie_adjustment, n)
    end
end
SignedRankTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S}) =
    SignedRankTest(x - y)

## EXACT WILCOXON SIGNED RANK TEST

immutable ExactSignedRankTest{T <: Real} <: HypothesisTest
    W::Float64
    ranks::Vector{T}
    tie_adjustment::Float64
end
ExactSignedRankTest{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) =
    ExactSignedRankTest(x - y)
ExactSignedRankTest{T <: Real}(x::Vector{T}) =
    ExactSignedRankTest(signedrankstats(x)...)

testname(::ExactSignedRankTest) = "Exact Wilcoxon signed rank test"

# Enumerate all possible Mann-Whitney U results for a given vector, determining left-
# and right-tailed p values
function signedrankenumerate(x::ExactSignedRankTest)
    le = 0
    gr = 0
    n = length(x.ranks)
    tot = 2^n
    for i = 0:tot-1
        # Interpret bits of i as signs to generate wp for all possible sign combinations
        Wp = 0
        b = i
        j = 1
        while b != 0
            Wp += (b & 1)*x.ranks[j]
            j += 1
            b >>= 1
        end
        le += Wp <= x.W
        gr += Wp >= x.W
    end
    (le/tot, gr/tot)
end

pvalue(x::ExactSignedRankTest; tail=:both) = 
    if length(x.ranks) == 0
        1
    elseif x.tie_adjustment == 0
        # Compute exact p-value using method from Rmath, which is fast but cannot account for ties
        n = length(x.ranks)
        if tail == :both
            if x.W <= n * (n + 1)/4
                p = 2 * psignrank(x.W, n, true)
            else
                p = 2 * psignrank(x.W - 1, n, false)
            end
        elseif tail == :left
            psignrank(x.W, length(x.ranks), true)
        elseif tail == :right
            psignrank(x.W - 1, length(x.ranks), false)
        end
    else
        # Compute exact p-value by enumerating all possible ranks in the tied data
        if tail == :both
            min(1, 2 * min(signedrankenumerate(x)))
        elseif tail == :left
            singedrankenumerate(x.W, x.ranks)[1]
        elseif tail == :right
            signedrankenumerate(x.W, x.ranks)[2]
        end
    end
    

## APPROXIMATE SIGNED RANK TEST

ApproximateSignedRankTest(W::Float64, tie_adjustment::Float64, n::Int) =
    ApproximateSignedRankTest(W, tie_adjustment, W - n * (n + 1)/4,
        sqrt(n * (n + 1) * (2 * n + 1) / 24 - tie_adjustment / 48))
ApproximateSignedRankTest{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) =
    ApproximateSignedRankTest(x - y)
function ApproximateSignedRankTest{T <: Real}(x::Vector{T})
    (W, ranks, tie_adjustment) = signedrankstats(x)
    ApproximateSignedRankTest(W, tie_adjustment, length(ranks))
end

testname(::Type{ApproximateSignedRankTest}) = "Approximate Wilcoxon signed rank test"