# Rayleigh.jl
# Rayleigh test of randomness against a unimodal alternative
#
# For reference see:
# Statistical Analysis of Circular Data. Cambridge: Cambridge University
# Press, 1995.
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

export RayleighTest

# Rayleigh test of uniformity against a unimodal alternative
type RayleighTest <: HypothesisTest
	Rbar::Float64
	Z::Float64
	p_value::Float64
end
test_name(::Type{RayleighTest}) = "Rayleigh test"

# Complex numbers
test_statistic{S <: Complex}(::Type{RayleighTest}, samples::Vector{S}) =
	abs(sum(samples./abs(samples)))^2/length(samples)
# Angles (in radians)
test_statistic{S <: Real}(::Type{RayleighTest}, theta::Vector{S}) =
	abs(sum(exp(im*theta)))^2/length(theta)

# Z given
function p_value(::Type{RayleighTest}, Z::Float64, n::Int)
	if n < 1e6
		p = exp(-Z)*(1+(2*Z-Z^2)/(4*n)-(24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*n^2))
	else	# Avoid overflow
		p = exp(-Z)
	end
	p
end
p_value{S <: Number}(::Type{RayleighTest}, samples::Vector{S}) =
	p_value(RayleighTest, test_statistic(RayleighTest, samples), length(samples))

function RayleighTest{S <: Number}(samples::Vector{S})
	if isa(eltype(samples), Complex)
		s = abs(sum(samples./abs(samples)))
	else
		s = abs(sum(exp(im*samples)))
	end
	n = length(samples)
	Z = s^2/n
	p = p_value(RayleighTest, Z, length(samples))
	RayleighTest(s/n, Z, p)
end