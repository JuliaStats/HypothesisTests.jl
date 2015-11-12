# Anderson-Darling test

export OneSampleADTest, KSampleADTest

abstract ADTest <: HypothesisTest

## ONE SAMPLE AD-TEST
### http://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm

function adstats{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)
    n = length(x)
    y = sort(x)
    μ = mean(y)
    σ = std(y)
    z = cdf(d, (y-μ)/σ)
    i = 1:n
    S = ((2*i-1.0)/n) .* (log(z[i])+log(1-z[n+1-i]))
    S[isinf(S)] = 0. # remove infinity
    A² = -n-sum(S)
    (n, μ, σ, A²)
end

immutable OneSampleADTest <: ADTest
    n::Int      # number of observations
    μ::Float64  # sample mean
    σ::Float64  # sample std
    A²::Float64 # Anderson-Darling test statistic
end

function OneSampleADTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)
    OneSampleADTest(adstats(x, d)...)
end

testname(::OneSampleADTest) = "One sample Anderson-Darling test"

function show_params(io::IO, x::OneSampleADTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "sample mean:              $(x.μ)")
    println(io, ident, "sample SD:                $(x.σ)")
    println(io, ident, "A² statistic:             $(x.A²)")
end

function pvalue(x::OneSampleADTest)
    z = x.A²*(1.0 + .75/x.n + 2.25/x.n/x.n)

    # q = [.05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, .975, .99, .995]
    # b0 = [-.512, -.552, -.608, -.643, -.707, -.735, -.772, -.770, -.778, -.779, -.803, -.818, -.818, -.801, -.800, -.756, -.749, -.750, -.795, -.881, -1.013, -1.063]
    # b1 =[2.10, 1.25, 1.07, .93, 1.03, 1.02, 1.04, .90, .80, .67, .70, .58, .42, .12, -.09, -.39, -.59, -.8, -.89, -.94, -.93, -1.34]
    # a₊ = [.1674, .1938, .2147, .2333, .2509, .2681, .2853, .3030, .3213, .3405, .3612, .3836, .4085, .4367, .4695, .5091, .5597, .6305, .7514, .8728, 1.0348, 1.1578]
    # a = a₊.*(1.+ b0/x.n +b1/x.n/x.n)
    # qi = searchsortedlast(a, x.A²)
    # Z(i) = a₊[i]*(1+b0[i]/x.n +b1[i]/x.n/x.n)

    Z = [.200, .340, .600]
    if z < Z[1]
        1.0 - exp(-13.436+101.14z-223.73z^2)
    elseif Z[1] < z < Z[2]
        1.0 - exp(-8.318+42.796z-59.938z^2)
    elseif Z[2] < z < Z[3]
        exp(0.9177-4.279z-1.38z^2)
    else
        exp(1.2937-5.709z+0.0186z^2)
    end
end


# ### K-SAMPLE ANDERSON DARLING TEST
# function a2_ksample(samples, ties=true)
#     k = length(samples)
#     k < 2 && error("Need at least two samples")

#     n = map(length, samples)
#     Z = sort(vcat(samples...))
#     N = length(Z)
#     Zstar = unique(Z)

#     length(Zstar) < 2 && error("Need more then 1 observation")
#     any(n .== 0) && error("One of the samples is empty")

#     A2kN = 0.
#     if ties
#         Zssl = map(i->searchsortedfirst(Z, i), Zstar) -1
#         lj = N == length(Zstar) ? 1. : map(i->searchsortedlast(Z, i), Zstar) - Zssl
#         Bj = Zssl + lj / 2.
#         for i in 1:k
#             smpl = sort(samples[i])
#             ssr = map(i->searchsortedlast(smpl, i), Zstar)
#             Mij = map(Float64, ssr)
#             fij = ssr - (map(i->searchsortedfirst(smpl, i), Zstar)-1)
#             Mij -= fij / 2.
#             inner = (lj / N) .* (N * Mij - Bj * n[i]).^2 ./ (Bj .* (N - Bj) - N * lj / 4.)
#             A2kN += sum(inner) / n[i]
#         end
#         A2kN *= (N - 1.) / N
#     else
#         0.0
#         # lj = Z.searchsorted(Zstar[:-1], 'right') - Z.searchsorted(Zstar[:-1],
#         #                                                       'left')
#         # Bj = lj.cumsum()
#         # for i in arange(0, k):
#         #     s = np.sort(samples[i])
#         #     Mij = s.searchsorted(Zstar[:-1], side='right')
#         #     inner = lj / float(N) * (N * Mij - Bj * n[i])**2 / (Bj * (N - Bj))
#         #     A2kN += inner.sum() / n[i]
#     end

#     h = sum(1./(1:N-1))
#     H = sum(1./n)
#     g = 0
#     for l in 1:N-2
#         inner = [1. / ((N - l) * m) for m in (l+1):(N-1)]
#         g += sum(inner)
#     end

#     a = (4*g - 6) * (k - 1) + (10 - 6*g)*H
#     b = (2*g - 4)*k^2 + 8*h*k + (2*g - 14*h - 4)*H - 8*h + 4*g - 6
#     c = (6*h + 2*g - 2)*k^2 + (4*h - 4*g + 6)*k + (2*h - 6)*H + 4*h
#     d = (2*h + 6)*k^2 - 4*h*k
#     var_kn = (a*N^3 + b*N^2 + c*N + d) / ((N - 1.) * (N - 2.) * (N - 3.))
#     m = k - 1
#     A2 = (A2kN - m) / sqrt(var_kn)

#     b0 = [0.675, 1.281, 1.645, 1.96, 2.326]
#     b1 = [-0.245, 0.25, 0.678, 1.149, 1.822]
#     b2 = [-0.105, -0.305, -0.362, -0.391, -0.396]
#     critical = b0 + b1 / sqrt(m) + b2 / m
#     f = poly_fit(critical, log([0.25, 0.1, 0.05, 0.025, 0.01]), 2)
#     pf(x) = exp(f[1] + f[2]*x + f[3]*x^2)

#     return A2, critical, pf(A2)
# end
