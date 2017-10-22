# Anderson-Darling test

export OneSampleADTest, KSampleADTest

abstract type ADTest <: HypothesisTest end

## ONE SAMPLE AD-TEST
### http://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm

function adstats(x::AbstractVector{T}, d::UnivariateDistribution) where T<:Real
    n = length(x)
    y = sort(x)
    μ = mean(y)
    σ = std(y)
    A² = convert(typeof(μ), -n)
    for i = 1:n
        zi = (y[i] - μ)/σ
        zni1 = (y[n - i + 1] - μ)/σ
        lcdfz = logcdf(d, zi)
        lccdfz = logccdf(d, zni1)
        A² -= (2*i - 1)/n * (lcdfz + lccdfz)
    end
    (n, μ, σ, A²)
end

struct OneSampleADTest <: ADTest
    n::Int      # number of observations
    μ::Float64  # sample mean
    σ::Float64  # sample std
    A²::Float64 # Anderson-Darling test statistic
end

"""
    OneSampleADTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)

Perform a one-sample Anderson–Darling test of the null hypothesis that the data in vector
`x` come from the distribution `d` against the alternative hypothesis that the sample
is not drawn from `d`.

Implements: [`pvalue`](@ref)
"""
function OneSampleADTest(x::AbstractVector{T}, d::UnivariateDistribution) where T<:Real
    OneSampleADTest(adstats(x, d)...)
end

testname(::OneSampleADTest) = "One sample Anderson-Darling test"
default_tail(test::OneSampleADTest) = :right

function show_params(io::IO, x::OneSampleADTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "sample mean:              $(x.μ)")
    println(io, ident, "sample SD:                $(x.σ)")
    println(io, ident, "A² statistic:             $(x.A²)")
end

### Ralph B. D'Agostino, Goodness-of-Fit-Techniques, CRC Press, Jun 2, 1986
### https://books.google.com/books?id=1BSEaGVBj5QC&pg=PA97, p.127
function pvalue(x::OneSampleADTest)
    z = x.A²*(1.0 + .75/x.n + 2.25/x.n/x.n)

    if z < .200
        1.0 - exp(-13.436+101.14z-223.73z^2)
    elseif .200 < z < .340
        1.0 - exp(-8.318+42.796z-59.938z^2)
    elseif .340 < z < .600
        exp(0.9177-4.279z-1.38z^2)
    elseif z < 153.467
        exp(1.2937-5.709z+0.0186z^2)
    else
        0.0
    end
end

## K-SAMPLE ANDERSON DARLING TEST
struct KSampleADTest{T<:Real} <: ADTest
    k::Int             # number of samples
    n::Int             # number of observations
    σ::Float64         # variance A²k
    A²k::Float64       # Anderson-Darling test statistic
    modified::Bool     # Modified test statistic
    nsim::Int          # Number of simulations for P-value calculation (0 - for asymptotic calculation)
    samples::Vector{T} # Pooled samples
    sizes::Vector{Int} # sizes of samples
end

"""
    KSampleADTest(xs::AbstractVector{<:Real}...; modified = true, nsim=0)

Perform a ``k``-sample Anderson–Darling test of the null hypothesis that the data in the
``k`` vectors `xs` come from the same distribution against the alternative hypothesis that
the samples come from different distributions.

`modified` parameter enables a modified test calculation for samples whose observations
do not all coincide.

`nsim` parameter if equals to 0, specifies an asymptotic calculation of P-value.
If greater than zero, enables an estimation of P-values by generating `nsim` random splits
of the pooled data on ``k`` samples, evaluating the AD statistics for each split, and
computing proportion of simulated values >= observed AD values, which is reported as
a P-value estimate.

Implements: [`pvalue`](@ref)

# References

  * F. W. Scholz and M. A. Stephens, K-Sample Anderson-Darling Tests, Journal of the
    American Statistical Association, Vol. 82, No. 399. (Sep., 1987), pp. 918-924.
"""
KSampleADTest(xs::AbstractVector{T}...; modified=true, nsim=0) where T<:Real =
    a2_ksample(xs, modified, nsim)

testname(::KSampleADTest) = "k-sample Anderson-Darling test"
default_tail(test::KSampleADTest) = :right

function show_params(io::IO, x::KSampleADTest, ident="")
    println(io, ident, "number of samples:        $(x.k)")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "SD of A²k:                $(x.σ)")
    println(io, ident, "A²k statistic:            $(x.A²k)")
    println(io, ident, "Standardized statistic:   $((x.A²k - x.k + 1) / x.σ)")
    println(io, ident, "Modified:                 $(x.modified)")
    println(io, ident, "P-value calculation:      $(x.nsim == 0 ? "asymptotic" : "simulation" )")
    x.nsim != 0 && println(io, ident, "Number of simulations:    $(x.nsim)")
end

"""Monte-Carlo simulation of the p-value for AD test"""
function pvaluesim(x::KSampleADTest)
    Z = sort(x.samples)
    Z⁺ = unique(Z)

    cn = [0, cumsum(x.sizes)...]
    Xr = [cn[i]+1:cn[i+1] for i in 1:x.k]
    idxs = collect(1:x.n)
    IV = [view(idxs, Xr[i]) for i in 1:x.k]
    Xr = [view(x.samples, IV[i]) for i in 1:x.k]

    pv = 0
    for j in 1:x.nsim
        shuffle!(idxs)
        A²k, A²km = adkvals(Z⁺, x.n, Xr)
        adv = x.modified ? A²km : A²k
        adv >= x.A²k && (pv += 1)
    end
    return pv/x.nsim
end

"""Asymptotic evaluation of the p-value for AD test"""
function pvalueasym(x::KSampleADTest)
    m = x.k - 1
    Tk = (x.A²k - m) / x.σ
    sig = [0.25, 0.1, 0.05, 0.025, 0.01]
    b0 = [0.675, 1.281, 1.645, 1.96, 2.326]
    b1 = [-0.245, 0.25, 0.678, 1.149, 1.822]
    b2 = [-0.105, -0.305, -0.362, -0.391, -0.396]
    tm = b0 + b1 / sqrt(m) + b2 / m

    # Fit a quadratic polynomial to tm and log(sig)
    # Adapted from code by Paulo José Salz Jabardo
    A = zeros(eltype(tm), 5, 3)
    A[:,1] = 1.0
    @inbounds for i = 1:2, k = 1:5
        A[k,i+1] = A[k,i] * tm[k]
    end
    f = A \ log.(sig)

    pv = exp(f[1] + f[2]*Tk + f[3]*Tk^2)
    return pv > 1 ? 1.0 : pv # cap p-value
end

pvalue(x::KSampleADTest) = x.nsim == 0 ? pvalueasym(x) : pvaluesim(x)

function adkvals(Z⁺, N, samples)
    k = length(samples)
    n = map(length, samples)
    L = length(Z⁺)

    fij = zeros(Int, k, L)
    for i in 1:k
        for s in samples[i]
            fij[i, searchsortedfirst(Z⁺, s)] += 1
        end
    end
    ljs = sum(fij, 1)

    A²k = A²km = 0.
    for i in 1:k
        innerm = 0.
        inner = 0.
        Mij = 0.
        Bj = 0.
        for j = 1:L
            lj = ljs[j]
            Mij += fij[i, j]
            Bj += lj
            Maij = Mij - fij[i, j]/2.
            Baj = Bj - lj/2.
            innerm += lj/N * (N*Maij-n[i]*Baj)^2 / (Baj*(N-Baj) - N*lj/4.)
            if j < L
                inner += lj/N * (N*Mij-n[i]*Bj)^2 / (Bj*(N-Bj))
            end
        end
        A²km += innerm / n[i]
        A²k += inner / n[i]
    end
    return A²k, A²km * (N - 1.) / N
end

function a2_ksample(samples, modified, method)
    k = length(samples)
    k < 2 && error("Need at least two samples")

    n = map(length, samples)
    pooled = vcat(samples...)
    Z = sort(pooled)
    N = length(Z)
    Z⁺ = unique(Z)
    L = length(Z⁺)

    L < 2 && error("Need more then 1 observation")
    minimum(n) == 0 && error("One of the samples is empty")

    A²k, A²km = adkvals(Z⁺, N, samples)

    H = sum(map(i->1./i, n))
    h = sum(1./(1:N-1))
    g = 0.
    for i in 1:N-2
        for j in i+1:N-1
            g += 1. / ((N - i) * j)
        end
    end

    a = (4*g - 6)*(k - 1) + (10 - 6*g)*H
    b = (2*g - 4)*k^2 + 8*h*k + (2*g - 14*h - 4)*H - 8*h + 4*g - 6
    c = (6*h + 2*g - 2)*k^2 + (4*h - 4*g + 6)*k + (2*h - 6)*H + 4*h
    d = (2*h + 6)*k^2 - 4*h*k
    σ² = (a*N^3 + b*N^2 + c*N + d) / ((N - 1.) * (N - 2.) * (N - 3.))

    KSampleADTest(k, N, sqrt(σ²), (modified ? A²km : A²k), modified, method, pooled, [n...])
end
