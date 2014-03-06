# binomial.jl
# Binomial tests
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

export BinomialTest, SignTest

## BINOMIAL TEST

immutable BinomialTest <: HypothesisTest
    p::Float64
    x::Int
    n::Int

    BinomialTest(x::Real, n::Real, p::Real=0.5) = new(p, x, n)
end

BinomialTest(x::AbstractVector{Bool}, p=0.5) =
	BinomialTest(sum(x), length(x), float64(p))

testname(::BinomialTest) = "Binomial test"
pvalue(x::BinomialTest; tail=:both) = pvalue(Binomial(x.n, x.p), x.x; tail=tail)

# Clopper-Pearson interval
function ci(x::BinomialTest, alpha::Float64=0.05; tail=:both)
    check_alpha(alpha)

    if tail == :both
        (quantile(Beta(x.x, x.n - x.x + 1), alpha/2), quantile(Beta(x.x + 1, x.n - x.x), 1-alpha/2))
    elseif tail == :left
        (0.0, quantile(Beta(x.x + 1, x.n - x.x), 1-alpha))
    elseif tail == :right
        (quantile(Beta(x.x, x.n - x.x + 1), alpha), 1.0)
    else
        error("tail=$(tail) is invalid")
    end
end

## SIGN TEST

immutable SignTest <: HypothesisTest
	median::Float64
	x::Int
	n::Int
end

SignTest{T<:Real}(x::AbstractVector{T}, median::Real=0) =
	SignTest(median, sum(x .> median), sum(x .!= median))
SignTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}) = SignTest(x - y, 0)

pvalue(x::SignTest; tail=:both) = pvalue(Binomial(x.n, 0.5), x.x; tail=tail)

testname(::SignTest) = "Sign test"
