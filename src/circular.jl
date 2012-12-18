# Rayleigh.jl
# Rayleigh test of randomness against a unimodal alternative
#
# For reference see:
# Fisher, N. I., & Lee, A. J. (1983). A correlation coefficient for
#     circular data. Biometrika, 70(2), 327–332. doi:10.2307/2335547
# Fisher, N. I. Statistical Analysis of Circular Data. Cambridge:
#     Cambridge University Press, 1995.
# Jammmaladak, S. R. & Sengupta, A. Topics in circular statistics
#     vol. 5.  World Scientific, 2001.
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

export RayleighTest, FisherTLinearAssociation, JammalamadakaCircularCorrelation

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

## N.I. FISHER'S TEST OF T-LINEAR CIRCULAR-CIRCULAR ASSOCIATION

type FisherTLinearAssociation <: HypothesisTest
	rho_t::Float64
	p_value::Float64
end
test_name(::Type{FisherTLinearAssociation}) = "T-linear test of circular-circular association"

# Calculate T from Fisher, 1993
function tlinear_T{S <: Real, T <: Real}(theta::Vector{S}, phi::Vector{T})
	A = sum(cos(theta).*cos(phi))
	B = sum(sin(theta).*sin(phi))
	C = sum(cos(theta).*sin(phi))
	D = sum(sin(theta).*cos(phi))
	A*B-C*D
end

# For large samples, compute the distribution and statistic of T
function tlinear_Z{S <: Real, T <: Real}(rho_t::FloatingPoint, theta::Vector{S}, phi::Vector{T})
	n = length(theta)
	theta_resultant = sum(exp(im*theta))
	phi_resultant = sum(exp(im*phi))
	theta_resultant_angle = angle(theta_resultant)
	phi_resultant_angle = angle(phi_resultant)
	alpha_2_theta = mean(cos(2*(theta-theta_resultant_angle)))
	beta_2_theta = mean(sin(2*(theta-theta_resultant_angle)))
	alpha_2_phi = mean(cos(2*(phi-phi_resultant_angle)))
	beta_2_phi = mean(sin(2*(phi-phi_resultant_angle)))
	U_theta = (1-alpha_2_theta^2-beta_2_theta^2)/2
	U_phi = (1-alpha_2_phi^2-beta_2_phi^2)/2
	V_theta = (abs2(theta_resultant)/n^2)*(1-alpha_2_theta)
	V_phi = (abs2(phi_resultant)/n^2)*(1-alpha_2_phi)
	sqrt(n)*U_theta*U_phi*rho_t/sqrt(V_theta*V_phi)

	# Alternative computational strategy from Fisher and Lee (1983)
	# a1 = [mean(cos(theta)), mean(cos(phi))]
	# b1 = [mean(sin(theta)), mean(sin(phi))]
	# a2 = [mean(cos(2*(theta))), mean(cos(2*(phi)))]
	# b2 = [mean(sin(2*(theta))), mean(sin(2*(phi)))]
	
	# mu = (1-a2.^2-b2.^2)/2
	# A = a1.^2 + b1.^2 + a2.*b1.^2 - a1.^2.*a2 - 2a1.*b1.*b2

	# Z = sqrt(n)*prod(mu)*rho_t/sqrt(prod(A))
end

# Angles (in radians)
function test_statistic{S <: Real, T <: Real}(::Type{FisherTLinearAssociation}, theta::Vector{S}, phi::Vector{T})
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
		function $(fn){S <: Real, T <: Real}(::Type{FisherTLinearAssociation}, theta::Vector{S}, phi::Vector{T}, uniformly_distributed::Bool...)
			check_same_length(theta, phi)
			n = length(theta)
			if n == 0
				return NaN
			elseif n < 25 || length(uniformly_distributed) == 0
				# If the number of samples is small, or if we don't know whether the 
				# distribution is uniform, use a permutation test.
				
				# "For n < 25, use a randomisation test based on the quantity T = AB - CD"
				ct = cos(theta)
				st = sin(theta)
				cp = cos(phi)
				sp = sin(phi)
				T = sum(ct.*cp)*sum(st.*sp)-sum(ct.*sp)*sum(st.*cp)
				greater = 0

				if n <= 8
					# Exact permutation test
					nperms = factorial(n)
					indices = [1:n]
					for i = 1:nperms
						tp = nthperm(indices, i)
						ctp = ct[tp]
						stp = st[tp]
						greater += $(comparison)($(transform)(sum(ctp.*cp)*sum(stp.*sp)-sum(ctp.*sp)*sum(stp.*cp)), T)
					end
					return greater / nperms
				end

				# Approximate permutation test
				const nperms = 10000
				thetap = copy(theta)
				for i = 1:nperms
					tp = randperm(n)
					ctp = ct[tp]
					stp = st[tp]
					greater += $(comparison)($(transform)(sum(ctp.*cp)*sum(stp.*sp)-sum(ctp.*sp)*sum(stp.*cp)), T)
				end
				return greater / nperms
			else
				# Use approximate distribution of test statistic. This only works if we
				# know whether the distributions of theta and phi are uniform. Otherwise,
				# the provided p-values are not conservative.

				rho_t = test_statistic(FisherTLinearAssociation, theta, phi)

				# "If either distribution has a mean resultant length 0...the statistic has
				# approximately a double exponential distribution with density 1/2*exp(-abs(x))
				(dist, stat) = uniformly_distributed[1] ? (Laplace(), n*rho_t) :
					(Normal(), tlinear_Z(rho_t, theta, phi))
				p = 2 * $(distfn)(dist, $(transform)(stat))
			end
		end
	end
end

FisherTLinearAssociation{S <: Number, T <: Number}(theta::Vector{S}, phi::Vector{T}) =
	FisherTLinearAssociation(test_statistic(FisherTLinearAssociation, theta, phi), p_value(FisherTLinearAssociation, theta, phi))

## JAMMALADAK's CIRCULAR CORRELATION

type JammalamadakaCircularCorrelation
	r::Float64
	p_value::Float64
end

test_statistic{S <: Real, T <: Real}(::Type{JammalamadakaCircularCorrelation}, alpha::Vector{S}, beta::Vector{T}) = 
	test_statistic(JammalamadakaCircularCorrelation, alpha, beta, angle(sum(exp(im*alpha))), angle(sum(exp(im*beta))))

function test_statistic{S <: Real, T <: Real}(::Type{JammalamadakaCircularCorrelation}, alpha::Vector{S}, beta::Vector{T}, alpha_bar::Real, beta_bar::Real)
	check_same_length(alpha, beta)
	sum(sin(alpha - alpha_bar).*sin(beta - beta_bar))/sqrt(sum(sin(alpha - alpha_bar).^2)*sum(sin(beta - beta_bar).^2))
end

function jammalamadaka_Z{S <: Real, T <: Real}(alpha::Vector{S}, beta::Vector{T})
	check_same_length(alpha, beta)
	alpha_bar = angle(sum(exp(im*alpha)))
	beta_bar = angle(sum(exp(im*beta)))
	sin2_alpha = sin(alpha - alpha_bar).^2
	sin2_beta = sin(beta - beta_bar).^2
	lambda_20 = mean(sin2_alpha)
	lambda_02 = mean(sin2_beta)
	lambda_22 = mean(sin2_alpha.*sin2_beta)
	r = test_statistic(JammalamadakaCircularCorrelation, alpha, beta, alpha_bar, beta_bar)
	(sqrt(length(alpha)*lambda_20*lambda_02/lambda_22)*r, r)
end

let dist = Normal()
	global p_value, left_p_value, right_p_value
	p_value(::Type{JammalamadakaCircularCorrelation}, args...) =
		2*ccdf(dist, abs(jammalamadaka_Z(args...)[1]))
	left_p_value(::Type{JammalamadakaCircularCorrelation}, args...) =
		cdf(dist, jammalamadaka_Z(args...)[1])
	right_p_value(::Type{JammalamadakaCircularCorrelation}, args...) =
		ccdf(dist, jammalamadaka_Z(args...)[1])

	function JammalamadakaCircularCorrelation{S <: Real, T <: Real}(alpha::Vector{S}, beta::Vector{T})
		(Z, r) = jammalamadaka_Z(alpha, beta)
		JammalamadakaCircularCorrelation(r, 2*ccdf(dist, abs(Z)))
	end
end
JammalamadakaCircularCorrelation{S <: Complex, T <: Complex}(alpha::Vector{S}, beta::Vector{T}) =
	JammalamadakaCircularCorrelation(angle(alpha), angle(beta))

## GENERAL

check_same_length(x::Vector, y::Vector) = if length(x) != length(y)
		error("Vectors must be the same length")
	end

# Complex numbers
for fn in (:test_statistic, :p_value, :left_p_value, :right_p_value)
	@eval begin
		$(fn){S <: Complex, T <: Complex}(T::Union(Type{FisherTLinearAssociation}, Type{JammalamadakaCircularCorrelation}), x::Vector{S}, y::Vector{T}) =
			$(fn)(T, angle(x), angle(y))
	end
end
