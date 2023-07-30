export ShapiroWilkTest

#=
From:
PATRICK ROYSTON
Approximating the Shapiro-Wilk W-test for non-normality
*Statistics and Computing* (1992) **2**, 117-119
DOI: [10.1007/BF01891203](https://doi.org/10.1007/BF01891203)
=#

# Coefficients from Royston (1992)
for (s, c) in [(:C1, [0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056]),
               (:C2, [0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633]),
               (:C3, [0.5440, -0.39978, 0.025054, -0.0006714]),
               (:C4, [1.3822, -0.77857, 0.062767, -0.0020322]),
               (:C5, [-1.5861, -0.31082, -0.083751, 0.0038915]),
               (:C6, [-0.4803, -0.082676, 0.0030302]),
               (:C7, [0.164, 0.533]),
               (:C8, [0.1736, 0.315]),
               (:C9, [0.256, -0.00635]),
               (:G, [-2.273, 0.459])]
    @eval $(Symbol(:__RS92_, s))(x) = Base.Math.@horner(x, $(c...))
end

struct ShapiroWilkCoefs <: AbstractVector{Float64}
    N::Int
    A::Vector{Float64}
end

Base.size(SWc::ShapiroWilkCoefs) = (SWc.N,)
Base.IndexStyle(::Type{ShapiroWilkCoefs}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(SWc::ShapiroWilkCoefs, i::Integer)
    @boundscheck checkbounds(SWc, i)
    @inbounds if checkbounds(Bool, SWc.A, i)
        return SWc.A[i]
    elseif isodd(length(SWc)) && i == lastindex(SWc.A) + 1
        return zero(eltype(SWc))
    else
        return -SWc.A[SWc.N + 1 - i]
    end
end

function ShapiroWilkCoefs(N::Integer)
    if N < 3
        throw(ArgumentError("N must be greater than or equal to 3, got $N"))
    elseif N == 3 # exact
        return ShapiroWilkCoefs(N, [sqrt(2.0) / 2.0])
    else
        # Weisberg&Bingham 1975 statistic; store only positive half of m:
        # it is (anti-)symmetric; hence '2' factor below
        m = [-quantile(Normal(), (i - 3 / 8) / (N + 1 / 4)) for i in 1:div(N, 2)]
        mᵀm = 2sum(abs2, m)
        x = 1 / sqrt(N)
        a₁ = m[1] / sqrt(mᵀm) + __RS92_C1(x) # aₙ = cₙ + (...)
        if N ≤ 5
            ϕ = (mᵀm - 2m[1]^2) / (1 - 2a₁^2)
            m .= m ./ sqrt(ϕ) # A, but reusing m to save allocs
            m[1] = a₁
        else
            a₂ = m[2] / sqrt(mᵀm) + __RS92_C2(x) # aₙ₋₁ = cₙ₋₁ + (...)
            ϕ = (mᵀm - 2m[1]^2 - 2m[2]^2) / (1 - 2a₁^2 - 2a₂^2)
            m .= m ./ sqrt(ϕ) # A, but reusing m to save allocs
            m[1], m[2] = a₁, a₂
        end
        return ShapiroWilkCoefs(N, m)
    end
end

function swstat(X::AbstractVector{<:Real}, A::ShapiroWilkCoefs)
    if last(X) - first(X) < length(X) * eps()
        throw(ArgumentError("sample is constant (up to numerical accuracy)"))
    end
    AX = dot(view(A, 1:length(X)), X)
    m = mean(X)
    S² = sum(x -> abs2(x - m), X)
    W = AX^2 / S²
    return min(W, one(W)) # to guard against numeric errors
end

struct ShapiroWilkTest <: HypothesisTest
    coefs::ShapiroWilkCoefs # expectation of order statistics for Shapiro-Wilk test
    W::Float64              # test statistic
    censored::Int           # only the smallest N₁ samples were used
end

testname(::ShapiroWilkTest) = "Shapiro-Wilk normality test"
function population_param_of_interest(t::ShapiroWilkTest)
    return ("Squared correlation of data and expected order statistics of N(0,1) (W)",
            1.0, t.W)
end
default_tail(::ShapiroWilkTest) = :left
censored_ratio(t::ShapiroWilkTest) = t.censored / length(t.coefs)

function show_params(io::IO, t::ShapiroWilkTest, indent)
    l = 24
    println(io, indent, rpad("number of observations:", l), length(t.coefs))
    println(io, indent, rpad("censored ratio:", l), censored_ratio(t))
    return println(io, indent, rpad("W-statistic:", l), t.W)
end

function StatsAPI.pvalue(t::ShapiroWilkTest)
    n = length(t.coefs)
    W = t.W

    if iszero(censored_ratio(t))
        if n == 3 # exact by Shapiro&Wilk 1965
            # equivalent to 6/π * (asin(sqrt(W)) - asin(sqrt(3/4)))
            return 1 - 6acos(sqrt(W)) / π
        elseif n ≤ 11 # Royston 1992
            γ = __RS92_G(n)
            if γ ≤ log1p(-W)
                return zero(W)
            end
            w = -log(γ - log1p(-W))
            μ = __RS92_C3(n)
            σ = exp(__RS92_C4(n))
        elseif 12 ≤ n # Royston 1992
            w = log1p(-W)
            μ = __RS92_C5(log(n))
            σ = exp(__RS92_C6(log(n)))
        end
        return ccdf(Normal(μ, σ), w)
    else
        throw("censored samples not implemented yet")
        # to implement censored samples follow Royston 1993 Section 3.3
    end
end

"""
    ShapiroWilkTest(X::AbstractVector{<:Real},
                    swc::ShapiroWilkCoefs=ShapiroWilkCoefs(length(X));
                    sorted::Bool=issorted(X),
                    censored::Integer=0)

Perform a Shapiro-Wilk test of the null hypothesis that sample `X` come from a
normal distribution.

This implementation is based on the method by Royston (1992).
The calculation of the p-value is exact for sample size `N = 3`, and for ranges
`4 ≤ N ≤ 11` and `12 ≤ N ≤ 5000` (Royston 1992) two separate approximations
for p-values are used.

# Keyword arguments
The following keyword arguments may be passed.
* `sorted::Bool=issorted(X)`: to indicate that sample `X` is already sorted.
* `censored::Integer=0`: to censor the largest samples from `X`
  (so called upper-tail censoring)

Implements: [`pvalue`](@ref)

# Warning
As noted by Royston (1993), (approximated) W-statistic will be accurate
but returned p-values may not be reliable if either of these apply:
* Sample size is large  (`N > 2000`) or small (`N < 20`)
* Too much data is censored (`censored / N > 0.8`)

# Implementation notes
* The current implementation DOES NOT implement p-values for censored data.
* If multiple Shapiro-Wilk tests are to be performed on samples of same
  cardinality it is beneficial to pass `swc` for re-use.
* For maximal performance sorted `X` should be passed and indicated with
  `sorted=true` keyword argument.

# References
Shapiro, S. S., & Wilk, M. B. (1965). An Analysis of Variance Test for Normality
(Complete Samples). *Biometrika*, 52, 591–611.
[doi:10.1093/BIOMET/52.3-4.591](https://doi.org/10.1093/BIOMET/52.3-4.591).

Royston, P. (1992). Approximating the Shapiro-Wilk W-test for non-normality.
*Statistics and Computing*, 2(3), 117–119.
[doi:10.1007/BF01891203](https://doi.org/10.1007/BF01891203)

Royston, P. (1993). A Toolkit for Testing for Non-Normality in Complete and
Censored Samples. Journal of the Royal Statistical Society Series D
(The Statistician), 42(1), 37–43.
[doi:10.2307/2348109](https://doi.org/10.2307/2348109)

Royston, P. (1995). Remark AS R94: A Remark on Algorithm AS 181: The W-test for
Normality. *Journal of the Royal Statistical Society Series C
(Applied Statistics)*, 44(4), 547–551.
[doi:10.2307/2986146](https://doi.org/10.2307/2986146).
"""
function ShapiroWilkTest(sample::AbstractVector{<:Real},
                         swcoefs::ShapiroWilkCoefs=ShapiroWilkCoefs(length(sample));
                         sorted::Bool=issorted(sample),
                         censored::Integer=0)
    N = length(sample)
    if N < 3
        throw(ArgumentError("at least 3 samples are required, got $N"))
    elseif censored ≥ N
        throw(ArgumentError("`censored` must be less than `length(sample)`, " *
                            "got `censored = $censored > $N = length(sample)`"))
    elseif length(swcoefs) ≠ length(sample)
        throw(DimensionMismatch("length of sample and Shapiro-Wilk coefficients " *
                                "differ, got $N and $(length(swcoefs))"))
    end

    W = if !sorted
        swstat(sort!(sample[1:(end - censored)]), swcoefs)
    else
        swstat(@view(sample[1:(end - censored)]), swcoefs)
    end

    return ShapiroWilkTest(swcoefs, W, censored)
end

#=
# To compare with the standard ALGORITHM AS R94 fortran subroutine
#  * grab scipys (swilk.f)[https://raw.githubusercontent.com/scipy/scipy/main/scipy/stats/statlib/swilk.f];
#  * compile
#  ```
#  gfortran -shared -fPIC -o swilk.so swilk.f
#  gfortran -fdefault-integer-8 -fdefault-real-8 -shared -fPIC swilk.f -o swilk64.so
#  ```

for (lib, I, F) in (("./swilk64.so", Int64, Float64),
                    ("./swilk.so"  , Int32, Float32))
    @eval begin
        function swilkfort!(X::AbstractVector{$F}, A::AbstractVector{$F}, computeA=true)
            X = issorted(X) ? X : sort(X)
            w, pval = Ref{$F}(0.0), Ref{$F}(0.0)
            ifault = Ref{$I}(0)

            ccall((:swilk_, $lib),
                Cvoid,
                (
                Ref{Bool},  # INIT if false compute SWCoeffs in A, else use A
                Ref{$F},    # X    sample
                Ref{$I},    # N    samples length
                Ref{$I},    # N1   (upper) uncensored data length
                Ref{$I},    # N2   length of A
                Ref{$F},    # A    A
                Ref{$F},    # W    W-statistic
                Ref{$F},    # PW   p-value
                Ref{$I},    # IFAULT error code (see swilk.f for meaning)
                ),
                !computeA, X, length(X), length(X), div(N,2), A, w, pval, ifault)
            return (w[], pval[], ifault[], A)
        end
        swilkfort(X::Vector{$F}) = swilkfort!(X, zeros($F, div(length(X),2)))
    end
end
=#
