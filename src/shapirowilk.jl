export ShapiroWilkTest, SWCoeffs

#=
From:
PATRICK ROYSTON
Approximating the Shapiro-Wilk W-test for non-normality
*Statistics and Computing* (1992) **2**, 117-119
DOI: [10.1007/BF01891203](https://doi.org/10.1007/BF01891203)
=#

# TODO: Rerun simulation and polynomial fitting

const ROYSTON_COEFFS = Dict{String, Vector{Float64}}(
"C1" => [0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1],
"C2" => [0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1],
"C3" => [0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3],
"C4" => [0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2],
"C5" => [-0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2],
"C6" => [-0.4803E0, -0.82676E-1, 0.30302E-2],
"C7" => [0.164E0, 0.533E0],
"C8" => [0.1736E0, 0.315E0],
"C9" => [0.256E0, -0.635E-2],
"G"  => [-0.2273E1, 0.459E0]
)

for (s,c) in ROYSTON_COEFFS
    @eval $(Symbol("_"*s))(x) = Base.Math.@horner(x, $(c...))
end

#=
The following hardcoded constants has been replaced by more precise values:

SQRTH = sqrt(2.0)/2.0 # 0.70711E0
TH = 3/8 # 0.375E0
SMALL = eps(1.0) # 1E-19
PI6 = π/6 # 0.1909859E1
STQR = asin(sqrt(0.75)) # 0.1047198E1
=#

struct SWCoeffs
    N::Int
    A::Vector{Float64}
end

Base.length(SWc::SWCoeffs) = SWc.N
Base.endof(SWc::SWCoeffs) = length(SWc)

function Base.getindex(SWc::SWCoeffs, i::Int)
    if i <= endof(SWc.A)
        return SWc.A[i]
    elseif i <= endof(SWc)
        if isodd(SWc.N) && i == div(SWc.N, 2) + 1
            return 0.0
        else
            return -SWc.A[SWc.N+1-i]
        end
    else
        throw(BoundsError(SWc, i))
    end
end

function SWCoeffs(N::Int)
    if N == 3 # exact
        A = [sqrt(2.0)/2.0]
    elseif N > 3 # Weisberg&Bingham 1975 statistic
        m = [norminvcdf((i - 3/8)/(N + 1/4)) for i in 1:div(N,2)]

        mm = 2sum(abs2, m)
        sqrt_mm = sqrt(mm)

        x = 1/sqrt(N)
        a₁ = _C1(x) - m[1]/sqrt_mm

        if N ≤ 5
            ϕ = (mm - 2m[1]^2)/(1 - 2a₁^2)
            A = -m/sqrt(ϕ)
            A[1] = a₁
        else
            a₂ = -m[2]/sqrt_mm + _C2(x)
            ϕ = (mm - 2m[1]^2 - 2m[2]^2)/(1 - 2a₁^2 - 2a₂^2)
            A = -m/sqrt(ϕ)
            A[1], A[2] = a₁, a₂
        end
    else
        throw(ArgumentError("N must be greater than or equal to 3: $N"))
    end
    return SWCoeffs(N, A)
end

function swstat(X::AbstractArray{T}, A::SWCoeffs) where T<:Real

    if X[end] - X[1] < endof(X)*eps(eltype(X))
        throw("Data seems to be constant!")
    end

    AX = sum([A[i]*X[i] for i in 1:endof(X)])
    S² = sum(abs2, X-mean(X))

    return AX^2/S²
end

function pvalue(W::Float64, A::SWCoeffs, N1=A.N)
    N = A.N
    logN = log(A.N)
    if N == 3 # exact by Shapiro&Wilk 1965
        return π/6 * (asin(sqrt(W)) - asin(sqrt(0.75)))
    elseif N ≤ 11

        γ = _G(N)
        if log(1 - W) > γ
            return eps(Float64)
        end

        w = -log(γ - log(1 - W))
        m = _C3(N)
        sd = exp(_C4(N))
    else
        w = log(1-W)
        m = _C5(logN)
        sd = exp(_C6(logN))
    end

    if (N - N1) > 0
        throw("Not implemented yet!")
    end

    return normccdf((w - m)/sd)
end

struct ShapiroWilkTest <: HypothesisTest
    SWc::SWCoeffs            # Expectation of order statistics for Shapiro-Wilk test
    W::Float64            # test statistic
    N1::Int               #
end

testname(::ShapiroWilkTest) = "Shapiro-Wilk normality test"
population_param_of_interest(t::ShapiroWilkTest) =
("Squared correlation of data and SWCoeffes (W)", 1.0, t.W)

default_tail(::ShapiroWilkTest) = :left

function show_params(io::IO, t::ShapiroWilkTest, ident)
    println(io, ident, "number of observations:         ", t.SWc.N)
    println(io, ident, "censored ratio:                 ", (t.SWc.N - t.N1)/t.SWc.N)
    println(io, ident, "W-statistic:                    ", t.W)
end

pvalue(t::ShapiroWilkTest) = pvalue(t.W, t.SWc, t.N1)

function ShapiroWilkTest(
    X::AbstractArray{T}, SWc::SWCoeffs=SWCoeffs(length(X)), N1=length(X)) where {T<:Real}

    N = length(X)
    #fatal errors
    if N < 3
        throw("Need at least 3 samples.")
    elseif N1 > N
        throw("N1 must be less than or equal to length(X)")
    elseif length(SWc) ≠ length(X)
        throw("Length of the sample differs from Shapiro-Wilk coefficients!")
    end

    #non-fatal errors
    if N > 5000
        warn("p-value may be unreliable for samples larger than 5000 points")
    elseif (N1 < N) && (N < 20)
        warn("Number of samples is < 20. Censoring may produce unreliable p-value.")
    elseif (N - N1)/N > 0.8
        warn("(N - N1)/N > 0.8. Censoring too much data may produce unreliable p-value.")
    end

    if !issorted(view(X, 1:N1))
        warn("Shapiro-Wilk requires sorted data")
        W = swstat(sort(X[1:N1]), SWc)
    else
        W = swstat(view(X, 1:N1), SWc)
    end

    return ShapiroWilkTest(SWc, W, N1)
end

#=
# To compare with the standard ALGORITHM AS R94 fortran subroutine
#  * grab scipys (swilk.f)[https://github.com/scipy/scipy/blob/master/scipy/stats/statlib/swilk.f];
#  * compile (to use Int64/Float64 You need a few trivial changes to swilk.f)
#  ```
#  gfortran -shared -fPIC -o swilk.so swilk.f
#  gfortran -shared -fPIC -o swilk64.so swilk64.f
#  ```
#

for (lib, I, F) in (("./swilk64.so", Int64, Float64),
                    ("./swilk.so"  , Int32, Float32))
    @eval begin
        function swilkfort!(X::AbstractVector{$F}, A::AbstractVector{$F}, computeA=true)

            w, pval = Ref{$F}(0.0), Ref{$F}(0.0)
            ifault = Ref{$I}(0)

            ccall((:swilk_, $lib),
                Void,
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
                !computeA, X, length(X), length(X), length(A), A, w, pval, ifault)
            return (w[], pval[], ifault[], A)
        end
        swilkfort(X::Vector{$F}) = swilkfort!(X, zeros($F, div(length(X),2)))
    end
end
=#
