# Anderson-Darling test

export OneSampleADTest, KSampleADTest

@compat abstract type ADTest <: HypothesisTest end

## ONE SAMPLE AD-TEST
### http://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm

function adstats{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)
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
### k-Sample Anderson-Darling Tests, F. W. Scholz; M. A. Stephens, Journal of the American Statistical Association, Vol. 82, No. 399. (Sep., 1987), pp. 918-924.

immutable KSampleADTest <: ADTest
    k::Int        # number of samples
    n::Int        # number of observations
    σ::Float64   # variance A²k
    A²k::Float64 # Anderson-Darling test statistic
end

function KSampleADTest{T<:Real}(xs::AbstractVector{T}...; modified=true)
    KSampleADTest(a2_ksample(xs, modified)...)
end

testname(::KSampleADTest) = "k-sample Anderson-Darling test"

function show_params(io::IO, x::KSampleADTest, ident="")
    println(io, ident, "number of samples:        $(x.k)")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "SD of A²k:               $(x.σ)")
    println(io, ident, "A²k statistic:           $(x.A²k)")
end

function pvalue(x::KSampleADTest)
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

    exp(f[1] + f[2]*Tk + f[3]*Tk^2)
end

function a2_ksample(samples, modified=true)
    k = length(samples)
    k < 2 && error("Need at least two samples")

    n = map(length, samples)
    Z = sort(vcat(samples...))
    N = length(Z)
    Z⁺ = unique(Z)
    L = length(Z⁺)

    L < 2 && error("Need more then 1 observation")
    minimum(n) == 0 && error("One of the samples is empty")

    fij = zeros(Int, k, L)
    for i in 1:k
        for s in samples[i]
            fij[i, searchsortedfirst(Z⁺, s)] += 1
        end
    end

    A²k = 0.
    if modified
        for i in 1:k
            inner = 0.
            Mij = 0.
            Bj = 0.
            for j = 1:L
                lj = sum(fij[:,j])
                Mij += fij[i, j]
                Bj += lj
                Maij = Mij - fij[i, j]/2.
                Baj = Bj - lj/2.
                inner += lj/N * (N*Maij-n[i]*Baj)^2 / (Baj*(N-Baj) - N*lj/4.)
            end
            A²k += inner / n[i]
        end
        A²k *= (N - 1.) / N
    else
        for i in 1:k
            inner = 0.
            Mij = 0.
            Bj = 0.
            for j = 1:L-1
                lj = sum(fij[:,j])
                Mij += fij[i, j]
                Bj += lj
                inner += lj/N * (N*Mij-n[i]*Bj)^2 / (Bj*(N-Bj))
            end
            A²k += inner / n[i]
        end
    end

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
    σ²n = (a*N^3 + b*N^2 + c*N + d) / ((N - 1.) * (N - 2.) * (N - 3.))

    (k, N, sqrt(σ²n), A²k)
end
