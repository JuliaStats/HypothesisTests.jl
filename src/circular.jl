# Rayleigh.jl
# Rayleigh test of randomness against a unimodal alternative
#
# For reference see:
# Fisher, N. I., & Lee, A. J. (1983). A correlation coefficient for
#     circular data. Biometrika, 70(2), 327–332. doi:10.2307/2335547
# Fisher, N. I. Statistical Analysis of Circular Data. Cambridge:
#     Cambridge University Press, 1995.
# Jammalamadaka, S. R. & SenGupta, A. Topics in circular statistics
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

struct RayleighTest <: HypothesisTest
    Rbar::Float64 # mean resultant length
    n::Int        # number of observations
end
function RayleighTest(samples::Vector{S}) where S <: Complex
    s = Float64(abs(sum(samples./abs(samples))))
    n = length(samples)
    Rbar = s/n
    RayleighTest(Rbar, n)
end
function RayleighTest(samples::Vector{S}) where S <: Real
    s = Float64(abs(sum(exp, im * samples)))
    n = length(samples)
    Rbar = s/n
    RayleighTest(Rbar, n)
end

testname(::RayleighTest) = "Rayleigh test"
population_param_of_interest(x::RayleighTest) = ("Mean resultant length", 0, x.Rbar) # parameter of interest: name, value under h0, point estimate
default_tail(test::RayleighTest) = :both

function show_params(io::IO, x::RayleighTest, ident="")
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "test statistic:         $(x.Rbar^2 * x.n)")
end

function pvalue(x::RayleighTest)
    Z = x.Rbar^2 * x.n
    x.n > 1e6 ? exp(-Z) :
        exp(-Z)*(1+(2*Z-Z^2)/(4*x.n)-(24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*x.n^2))
end

## N.I. FISHER'S TEST OF T-LINEAR CIRCULAR-CIRCULAR ASSOCIATION

struct FisherTLinearAssociation{S <: Real, T <: Real} <: HypothesisTest
    rho_t::Float64                              # circular correlation coefficient
    theta::Vector{S}                            # radians of group 1
    phi::Vector{T}                              # radians of group 2
    uniformly_distributed::Union{Bool,Nothing}     # is distribution of theta and phi uniform?
end
function FisherTLinearAssociation(theta::Vector{Stheta},
        phi::Vector{Sphi}, uniformly_distributed::Union{Bool,Nothing}) where {Stheta <: Real, Sphi <: Real}
    check_same_length(theta, phi)

    A = sum(cos.(theta).*cos.(phi))
    B = sum(sin.(theta).*sin.(phi))
    C = sum(cos.(theta).*sin.(phi))
    D = sum(sin.(theta).*cos.(phi))
    T = A*B-C*D

    # Notation drawn from Fisher, 1993
    n = length(theta)
    E = sum(cos, 2*theta)
    F = sum(sin, 2*theta)
    G = sum(cos, 2*phi)
    H = sum(sin, 2*phi)
    rho_t = 4*T/sqrt((n^2 - E^2 - F^2)*(n^2-G^2-H^2))
    FisherTLinearAssociation(rho_t, theta, phi, uniformly_distributed)
end
FisherTLinearAssociation(theta::Vector{S},
phi::Vector{T}) where {S <: Real, T <: Real} = FisherTLinearAssociation(theta, phi, nothing)

testname(::FisherTLinearAssociation) =
    "T-linear test of circular-circular association"
population_param_of_interest(x::FisherTLinearAssociation) = ("Circular correlation coefficient", 0, x.rho_t) # parameter of interest: name, value under h0, point estimate
default_tail(test::FisherTLinearAssociation) = :both

function show_params(io::IO, x::FisherTLinearAssociation, ident="")
    println(io, ident, "number of observations: [$(length(x.theta)),$(length(x.phi))]")
end

# For large samples, compute the distribution and statistic of T
function tlinear_Z(x::FisherTLinearAssociation)
    n = length(x.theta)
    theta_resultant       = sum(exp, im * x.theta)
    phi_resultant         = sum(exp, im * x.phi)
    theta_resultant_angle = angle(theta_resultant)
    phi_resultant_angle   = angle(phi_resultant)
    alpha_2_theta         = mean(cos, 2 * (x.theta - theta_resultant_angle))
    beta_2_theta          = mean(sin, 2 * (x.theta - theta_resultant_angle))
    alpha_2_phi           = mean(cos, 2 * (x.phi - phi_resultant_angle))
    beta_2_phi            = mean(sin, 2 * (x.phi - phi_resultant_angle))
    U_theta               = (1 - alpha_2_theta^2 - beta_2_theta^2) / 2
    U_phi                 = (1-alpha_2_phi^2-beta_2_phi^2)/2
    V_theta               = (abs2(theta_resultant) / n^2) * (1 - alpha_2_theta)
    V_phi                 = (abs2(phi_resultant) / n^2) * (1 - alpha_2_phi)
    return sqrt(n) * U_theta * U_phi * rho_t / sqrt(V_theta * V_phi)

    # Alternative computational strategy from Fisher and Lee (1983)
    # a1 = [mean(cos(theta)), mean(cos(phi))]
    # b1 = [mean(sin(theta)), mean(sin(phi))]
    # a2 = [mean(cos(2*(theta))), mean(cos(2*(phi)))]
    # b2 = [mean(sin(2*(theta))), mean(sin(2*(phi)))]

    # mu = (1-a2.^2-b2.^2)/2
    # A = a1.^2 + b1.^2 + a2.*b1.^2 - a1.^2.*a2 - 2a1.*b1.*b2

    # Z = sqrt(n)*prod(mu)*rho_t/sqrt(prod(A))
end

# p-values
function pvalue(x::FisherTLinearAssociation; tail=:both)
    n = length(x.theta)
    if n == 0
        return NaN
    elseif n < 25 || x.uniformly_distributed == nothing
        # If the number of samples is small, or if we don't know whether the
        # distribution is uniform, use a permutation test.

        # "For n < 25, use a randomisation test based on the quantity T = AB - CD"
        ct = cos.(x.theta)
        st = sin.(x.theta)
        cp = cos.(x.phi)
        sp = sin.(x.phi)
        T = sum(ct.*cp)*sum(st.*sp)-sum(ct.*sp)*sum(st.*cp)
        greater = 0

        exact = n <= 8
        nperms = exact ? factorial(n) : 100000
        indices = [1:n;]
        for i = 1:nperms
            tp = exact ? nthperm(indices, i) : shuffle!(indices)
            a = 0.0
            b = 0.0
            c = 0.0
            d = 0.0
            for j = 1:n
                a += ct[tp[j]]*cp[j]
                b += st[tp[j]]*sp[j]
                c += ct[tp[j]]*sp[j]
                d += st[tp[j]]*cp[j]
            end
            Tp = a*b-c*d
            greater += tail == :both ? abs(Tp) > T :
                       tail == :left ? Tp < T :
                       tail == :right ? Tp > T :
                       throw(ArgumentError("tail=$(tail) is invalid"))
        end
        greater / nperms
    else
        # Use approximate distribution of test statistic. This only works if we
        # know whether the distributions of theta and phi are uniform. Otherwise,
        # the provided p-values are not conservative.

        # "If either distribution has a mean resultant length 0...the statistic has
        # approximately a double exponential distribution with density 1/2*exp(-abs(x))""
        (dist, stat) = x.uniformly_distributed ? (Laplace(), n*x.rho_t) :
            (Normal(), tlinear_Z(x.rho_t, x.theta, x.phi))
        return pvalue(dist, stat; tail=tail)
    end
end

## JAMMALAMADAKA'S CIRCULAR CORRELATION

struct JammalamadakaCircularCorrelation <: HypothesisTest
    r::Float64  # circular-circular correlation coefficient
    Z::Float64  # test statistic
end
function JammalamadakaCircularCorrelation(alpha::Vector{S}, beta::Vector{T}) where {S <: Real, T <: Real}
    check_same_length(alpha, beta)
    # calculate sample mean directions
    alpha_bar = angle(sum(exp, im * alpha))
    beta_bar = angle(sum(exp, im * beta))
    r = sum(sin.(alpha .- alpha_bar) .* sin.(beta .- beta_bar)) /
        sqrt(sum(sin.(alpha .- alpha_bar).^2) * sum(sin.(beta .- beta_bar).^2))

    sin2_alpha = sin.(alpha .- alpha_bar).^2
    sin2_beta = sin.(beta .- beta_bar).^2
    lambda_20 = mean(sin2_alpha)
    lambda_02 = mean(sin2_beta)
    lambda_22 = mean(sin2_alpha .* sin2_beta)
    Z = sqrt(length(alpha) * lambda_20 * lambda_02 / lambda_22) * r

    JammalamadakaCircularCorrelation(r, Z)
end

testname(::JammalamadakaCircularCorrelation) = "Jammalamadaka circular correlation"
population_param_of_interest(x::JammalamadakaCircularCorrelation) = ("Circular-circular correlation coefficient", 0, x.r) # parameter of interest: name, value under h0, point estimate
default_tail(test::JammalamadakaCircularCorrelation) = :both

function show_params(io::IO, x::JammalamadakaCircularCorrelation, ident="")
    println(io, ident, "test statistic: $(x.Z)")
end

pvalue(x::JammalamadakaCircularCorrelation; tail=:both) = pvalue(Normal(), x.Z; tail=tail)

## GENERAL

# Complex numbers
for fn in (:JammalamadakaCircularCorrelation, :FisherTLinearAssociation)
    @eval begin
        $(fn)(x::Vector{S}, y::Vector{T}) where {S <: Complex, T <: Complex} =
            $(fn)(angle(x), angle(y))
    end
end
