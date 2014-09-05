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
population_param_of_interest(x::BinomialTest) = ("Probability of success", x.p, x.x/x.n) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::BinomialTest, ident="")
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "number of successes:    $(x.x)")
end

pvalue(x::BinomialTest; tail=:both) = pvalue(Binomial(x.n, x.p), x.x; tail=tail)

# Confidence interval

function ci(x::BinomialTest, alpha::Float64=0.05; tail=:both, method=:clopper_pearson)
    check_alpha(alpha)

    if tail == :left
        (0.0, ci(x, alpha*2, method=method)[2])
    elseif tail == :right
        (ci(x, alpha*2, method=method)[1], 1.0)
    elseif tail == :both
        if method == :clopper_pearson
            ci_clopper_pearson(x, alpha)
        elseif method == :wald
            ci_wald(x, alpha)
        elseif method == :wilson
            ci_wilson(x, alpha)
        elseif method == :jeffrey
            ci_jeffrey(x, alpha)
        elseif method == :agresti_coull
            ci_agresti_coull(x, alpha)
        else
            error("method=$(method) is not implemented yet")
        end
    else
        error("tail=$(tail) is invalid")
    end
end

# Clopper-Pearson interval (confidence interval by inversion)
function ci_clopper_pearson(x::BinomialTest, alpha::Float64=0.05)
    (quantile(Beta(x.x, x.n - x.x + 1), alpha/2), quantile(Beta(x.x + 1, x.n - x.x), 1-alpha/2))
end

# Wald interval / normal approximation interval
function ci_wald(x::BinomialTest, alpha::Float64=0.05)
    μ = x.x / x.n
    σ = sqrt(μ*(1-μ)/x.n)
    (quantile(Normal(μ, σ), alpha/2), quantile(Normal(μ, σ), 1-alpha/2))
end

# Jeffreys interval
function ci_jeffrey(x::BinomialTest, alpha::Float64=0.05)
    (quantile(Beta(x.x + 1/2, x.n - x.x + 1/2), alpha/2), quantile(Beta(x.x + 1/2, x.n - x.x + 1/2), 1-alpha/2))
end

# Agresti-Coull interval
function ci_agresti_coull(x::BinomialTest, alpha::Float64=0.05)
    q = quantile(Normal(), 1-alpha/2)
    n = x.n + q^2
    μ = (x.x + q^2/2)/n
    σ = sqrt(μ*(1-μ)/n)
    (μ-q*σ, μ+q*σ)
end

# Wilson score interval
function ci_wilson(x::BinomialTest, alpha::Float64=0.05)
    q = quantile(Normal(), 1-alpha/2)
    p = x.x / x.n
    denominator = 1 + q^2/x.n
    μ = p + q^2/(2*x.n)
    μ /= denominator
    σ = sqrt(p*(1-p)/x.n + q^2/(4x.n^2))
    σ /= denominator
    (μ-q*σ, μ+q*σ)
end


## SIGN TEST

immutable SignTest <: HypothesisTest
	median::Float64
	x::Int
	n::Int
	data
end

SignTest{T<:Real}(x::AbstractVector{T}, median::Real=0) =
	SignTest(float64(median), sum(x .> median), sum(x .!= median), sort(x))
SignTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}) = SignTest(x - y, 0.0)

testname(::SignTest) = "Sign Test"
population_param_of_interest(x::SignTest) = ("Median", x.median, median(x.data)) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::SignTest, ident="")
    text1 = "number of observations:"
    text2 = "observations larger than $(x.median): "
    maxlen = length(text2)
    
    println(io, ident, text1, repeat(" ", maxlen-length(text1)), x.n)
    println(io, ident, text2, x.x)
end

pvalue(x::SignTest; tail=:both) = pvalue(Binomial(x.n, 0.5), x.x; tail=tail)

# confidence interval by inversion
function ci(x::SignTest, alpha::Float64=0.05; tail=:both)
	check_alpha(alpha)

	if tail == :left
		q = quantile(Binomial(x.n, 0.5), alpha)
		(x.data[q+1], x.median)
	elseif tail == :right
		q = quantile(Binomial(x.n, 0.5), alpha)
		(x.median, x.data[end-q])
	elseif tail == :both
		q = quantile(Binomial(x.n, 0.5), alpha/2)
		(x.data[q+1], x.data[end-q])
	else
		error("tail=$(tail) is invalid")
	end
end
