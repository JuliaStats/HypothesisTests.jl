# Rayleigh.jl
# Rayleigh test of randomness against a unimodal alternative
#
# For reference see:
# Fisher, N. I., & Lee, A. J. (1983). A correlation coefficient for
#     circular data. Biometrika, 70(2), 327–332. doi:10.2307/2335547
# Fisher, N. I., Statistical Analysis of Circular Data. Cambridge:
#     Cambridge University Press, 1995.
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

export RayleighTest, TLinearAssociation

## RAYLEIGH TEST OF UNIFORMITY AGAINST AN UNSPECIFIED UNIMODAL ALTERNATIVE

type RayleighTest <: HypothesisTest
	Rbar::Float64
	Z::Float64
	p_value::Float64
end
test_name(::Type{RayleighTest}) = "Rayleigh test"

# Complex numbers
test_statistic{S <: Complex}(::Type{RayleighTest}, samples::Vector{S}) =
	float64(abs(sum(samples./abs(samples))))^2/length(samples)
# Angles (in radians)
test_statistic{S <: Real}(::Type{RayleighTest}, theta::Vector{S}) =
	float64(abs(sum(exp(im*theta))))^2/length(theta)

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
		s = float64(abs(sum(samples./abs(samples))))
	else
		s = float64(abs(sum(exp(im*samples))))
	end
	n = length(samples)
	Z = s^2/n
	p = p_value(RayleighTest, Z, length(samples))
	RayleighTest(s/n, Z, p)
end

## TEST OF T-LINEAR CIRCULAR-CIRCULAR ASSOCIATION

type TLinearAssociation <: HypothesisTest
	rho_t::Float64
	p_value::Float64
end
test_name(::Type{TLinearAssociation}) = "T-linear test of circular-circular association"

check_same_length(x::Vector, y::Vector) = if length(x) != length(y)
		error("Vectors must be the same length")
	end

# Calculate T from Fisher, 1993
function tlinear_T{S <: Real, T <: Real}(theta::Vector{S}, phi::Vector{T})
	A = sum(cos(theta).*cos(phi))
	B = sum(sin(theta).*sin(phi))
	C = sum(cos(theta).*sin(phi))
	D = sum(sin(theta).*cos(phi))
	A*B-C*D
end

# For large samples, compute the distribution and statistic of T
function tlinear_large_sample_dist_and_stat{S <: Real, T <: Real}(theta::Vector{S}, phi::Vector{T})
	rho_t = test_statistic(TLinearAssociation, theta, phi)
	theta_resultant = sum(exp(im*theta))
	phi_resultant = sum(exp(im*phi))
	if abs(theta_resultant) <= eps() || abs(phi_resultant) <= eps()
		# "If either distribution has a mean resultant length 0...the statistic has
		# approximately a double exponential distribution with density 1/2*exp(-abs(x))
		return (Laplace(), rho_t)
	end

	alpha_2_theta = mean(cos(2*(theta-angle(theta_resultant))))
	beta_2_theta = mean(sin(2*(theta-angle(theta_resultant))))
	alpha_2_phi = mean(cos(2*(phi-angle(phi_resultant))))
	beta_2_phi = mean(sin(2*(phi-angle(phi_resultant))))
	U_theta = (1-alpha_2_theta^2-beta_2_theta^2)
	U_phi = (1-alpha_2_phi^2-beta_2_phi^2)
	V_theta = abs(theta_resultant/n)^2*(1-alpha_2_theta)
	V_phi = abs(phi_resultant/n)^2*(1-alpha_2_phi)
	Z = sqrt(n)*U_theta*U_phi*rho_t/sqrt(V_theta*V_phi)
	(Normal(), Z)
end

# Angles (in radians)
function test_statistic{S <: Real, T <: Real}(::Type{TLinearAssociation}, theta::Vector{S}, phi::Vector{T})
	check_same_length(theta, phi)
	# Notation drawn from Fisher, 1993
	n = length(theta)
	E = sum(cos(2*theta))
	F = sum(sin(2*theta))
	G = sum(cos(2*phi))
	H = sum(sin(2*phi))
	4*tlinear_T(theta, phi)/sqrt((n^2 - E^2 - F^2)*(n^2-G^2-H^2))
end

# p-values
for (fn, transform, comparison, distfn) in ((:p_value, :abs, :>, :ccdf),
	                                         (:left_p_value, :+, :<, :cdf),
	                                         (:right_p_value, :+, :>, :ccdf))
	@eval begin
		function $(fn){S <: Real, T <: Real}(::Type{TLinearAssociation}, theta::Vector{S}, phi::Vector{T})
			check_same_length(theta, phi)
			n = length(theta)
			if n < 25
				# "For n < 25, use a randomisation test based on the quantity T = AB - CD"
				T = $(transform)(tlinear_T(theta, phi))
				greater = 0

				if n <= 8
					# Exact permutation test
					nperms = factorial(n)
					for i = 1:maxn
						greater += $(comparison)($(transform)(tlinear_T(nthperm(theta, i), phi)), T)
					end
					return greater / nperms
				end

				# Approximate permutation test
				const nperms = 10000
				thetap = copy(theta)
				for i = 1:nperms
					greater += $(comparison)($(transform)(tlinear_T(shuffle!(thetap), phi)), T)
				end
				return greater / nperms
			end

			(dist, stat) = tlinear_large_sample_dist_and_stat(theta, phi)
			p = 2 * $(distfn)(dist, $(transform)(stat))
		end
	end
end

# Complex numbers
for fn in (:test_statistic, :p_value, :left_p_value, :right_p_value)
	@eval begin
		$(fn){S <: Complex, T <: Complex}(::Type{TLinearAssociation}, x::Vector{S}, y::Vector{T}) =
			$(fn)(angle(x), angle(y))
	end
end

TLinearAssociation{S <: Number, T <: Number}(theta::Vector{S}, phi::Vector{T}) =
	TLinearAssociation(test_statistic(TLinearAssociation, theta, phi), p_value(TLinearAssociation, theta, phi))