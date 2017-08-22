# fisher.jl
# Fisher's exact test
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

export FisherExactTest

immutable FisherExactTest <: HypothesisTest
    # Format:
    # X1  X2
    # Y1  a  b
    # Y2  c  d
    a::Int
    b::Int
    c::Int
    d::Int

    # conditional maximum likehood estimate of odd ratio
    ω::Float64

    function FisherExactTest(a::Int, b::Int, c::Int, d::Int)
        ω = cond_mle_odds_ratio(a, b, c, d)
        new(a, b, c, d, ω)
    end
end

testname(::FisherExactTest) = "Fisher's exact test"
population_param_of_interest(x::FisherExactTest) = ("Odds ratio", 1.0, x.ω) # parameter of interest: name, value under h0, point estimate
tail(test::FisherExactTest) = :both

# The sizing argument to print_matrix was removed during the 0.5 dev period
function _print_matrix(io::IO, X::AbstractVecOrMat, pre::AbstractString)
    Base.print_matrix(io, X, pre)
end

function show_params(io::IO, x::FisherExactTest, ident="")
    println(io, ident, "contingency table:")
    _print_matrix(io, [x.a x.b; x.c x.d], repeat(ident, 2))
    println(io)
end

# DOC: for tail=:both there exist multiple ``method``s for computing a pvalue and the corresponding ci.
function pvalue(x::FisherExactTest; tail=:both, method=:central)
    if tail == :both && method != :central
        if method == :minlike
            p = pvalue_both_minlike(x)
        else
            throw(ArgumentError("method=$(method) is not implemented yet"))
        end
    else
        p = pvalue(Hypergeometric(x.a + x.b, x.c + x.d, x.a + x.c), x.a, tail=tail)
    end
    p = max(min(p, 1.0), 0.0)

    return p
end

function pvalue_both_minlike(x::FisherExactTest, ω::Float64=1.0)
    a, b, c, d = reorder(x.a, x.b, x.c, x.d)

    dist = FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)

    p = pdf(dist, a)
    v = nextfloat(p)
    if a != 0
        p += cdf(dist, a-1)
    end

    # Add p-values of all tables in other tail equally or less probable
    for i = a+c:-1:a+1
        curp = pdf(dist, i)
        if curp > v
            break
        end
        p += curp
    end
    p
end

# confidence interval by inversion of p-value
function StatsBase.confint(x::FisherExactTest, alpha::Float64=0.05; tail=:both, method=:central)
    check_alpha(alpha)
    dist(ω) = FisherNoncentralHypergeometric(x.a+x.b, x.c+x.d, x.a+x.c, ω)
    obj(ω) = pvalue(dist(ω), x.a, tail=tail) - alpha

    if tail == :left # upper bound
        if (x.a == maximum(dist(1.0)))
            (0.0, Inf)
        else
            (0.0, fzero(obj, find_brackets(obj)...))
        end
    elseif tail == :right # lower bound
        if (x.a == minimum(dist(1.0)))
            (0.0, Inf)
        else
            (fzero(obj, find_brackets(obj)...), Inf)
        end
    elseif tail == :both
        if method == :central
            (StatsBase.confint(x, alpha/2; tail=:right)[1], StatsBase.confint(x, alpha/2; tail=:left)[2])
        else
            throw(ArgumentError("method=$(method) is not implemented yet"))
        end
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

## helpers

function reorder(a,b,c,d)
    if a + c > b + d
        a, b, c, d = b, a, d, c
    end
    if a/c > b/d
        a, b, c, d = c, d, a, b
    end
    (a, b, c, d)
end

# find values x_lower and x_upper, s.t. f(x_lower) < 0 and f(x_upper) > 0 or vice versa
function find_brackets(f::Function, x_init::Float64=1.0)
    f_init = f(x_init)

    if f_init > f(x_init + 1.0)
        find_brackets(x -> -f(x), x_init)
    else
        x_upper = x_lower = x_init
        if f_init > 0.0
            while x_lower > eps(0.0) && f(x_lower) > 0.0
                x_lower /= 2
            end
        else
            while f(x_upper) < 0.0
                x_upper *= 2
            end
        end
        (x_lower, x_upper)
    end
end

# find odds ratio ω that maximizes the (conditional) likelihood;
# since the mode and mean of Fisher's Noncentral Hypergeometric distribution
# coincide, this is equivalent to find ω, s.t., mean(dist(ω)) = a
function cond_mle_odds_ratio(a::Int, b::Int, c::Int, d::Int)
    dist(ω) = FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)

    if (a == minimum(dist(1.0)))
        0.0
    elseif a == maximum(dist(1.0))
        Inf
    else
        obj(ω) = mean(dist(ω))-a
        fzero(obj, find_brackets(obj)...)
    end
end
