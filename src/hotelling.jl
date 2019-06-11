# Hotelling T² test

export OneSampleHotellingT2Test, EqualCovHotellingT2Test, UnequalCovHotellingT2Test

abstract type HotellingT2TestTest <: HypothesisTest end

default_tail(::HotellingT2TestTest) = :right
pvalue(T::HotellingT2TestTest; tail=:right) = pvalue(FDist(dof(T)...), T.F, tail=tail)

function show_params(io::IO, T::HotellingT2TestTest, indent="")
    println(io, indent, "number of observations: ", nobs(T))
    println(io, indent, "number of variables:    ", T.p)
    println(io, indent, "T² statistic:           ", T.T²)
    println(io, indent, "transformed statistic:  ", T.F)
    println(io, indent, "degrees of freedom:     ", dof(T))
    println(io, indent, "covariance matrix:")
    Base.print_matrix(io, T.S, indent^2)
end

## Utility functions

function checkdims(X::AbstractMatrix, Y::AbstractMatrix)
    nx, p = size(X)
    ny, q = size(Y)
    p == q || throw(DimensionMismatch("Inconsistent number of variables: $p, $q"))
    (nx > 0 && ny > 0) || throw(ArgumentError("Inputs must be non-empty"))
    return (p, nx, ny)
end

# Helper function for computing A'B⁻¹A when B is a covariance matrix
At_Binv_A(A::AbstractArray, B::AbstractArray) = A'*(B \ A)

## One sample

struct OneSampleHotellingT2Test <: HotellingT2TestTest
    T²::Real
    F::Real
    n::Int
    p::Int
    μ₀::Vector
    x̄::Vector
    S::Matrix
end

StatsBase.nobs(T::OneSampleHotellingT2Test) = T.n
StatsBase.dof(T::OneSampleHotellingT2Test) = (T.p, T.n - T.p)

"""
    OneSampleHotellingT2Test(X::AbstractMatrix, μ₀=<zero vector>)

Perform a one sample Hotelling's ``T^2`` test of the hypothesis that the vector of
column means of `X` is equal to `μ₀`.
"""
function OneSampleHotellingT2Test(X::AbstractMatrix{T},
                              μ₀::AbstractVector=fill(middle(zero(T)), size(X, 2))) where T
    n, p = size(X)
    p == length(μ₀) ||
        throw(DimensionMismatch("Number of variables does not match number of means"))
    n > 0 || throw(ArgumentError("The input must be non-empty"))
    x̄ = vec(mean(X, dims=1))
    S = cov(X)
    T² = n * At_Binv_A(x̄ .- μ₀, S)
    F = (n - p) * T² / (p * (n - 1))
    return OneSampleHotellingT2Test(T², F, n, p, μ₀, x̄, S)
end

"""
    OneSampleHotellingT2Test(X::AbstractMatrix, Y::AbstractMatrix, μ₀=<zero vector>)

Perform a paired Hotelling's ``T^2`` test of the hypothesis that the vector of mean
column differences between `X` and `Y` is equal to `μ₀`.
"""
function OneSampleHotellingT2Test(X::AbstractMatrix{T}, Y::AbstractMatrix{S},
                              μ₀::AbstractVector=fill(middle(zero(T), zero(S)), size(X, 2))) where {T,S}
    p, nx, ny = checkdims(X, Y)
    nx == ny || throw(DimensionMismatch("Inconsistent number of observations: $nx, $ny"))
    p == length(μ₀) ||
        throw(DimensionMismatch("Number of variables does not match number of means"))
    return OneSampleHotellingT2Test(X .- Y, μ₀)
end

testname(::OneSampleHotellingT2Test) = "One sample Hotelling's T² test"
population_param_of_interest(T::OneSampleHotellingT2Test) = ("Mean vector", T.μ₀, T.x̄)

## Two sample

## Two sample: equal covariance

struct EqualCovHotellingT2Test <: HotellingT2TestTest
    T²::Real
    F::Real
    nx::Int
    ny::Int
    p::Int
    Δ::Vector
    S::Matrix
end

"""
    EqualCovHotellingT2Test(X::AbstractMatrix, Y::AbstractMatrix)

Perform a two sample Hotelling's ``T^2`` test of the hypothesis that the difference in
the mean vectors of `X` and `Y` is zero, assuming that `X` and `Y` have equal covariance
matrices.
"""
function EqualCovHotellingT2Test(X::AbstractMatrix, Y::AbstractMatrix)
    p, nx, ny = checkdims(X, Y)
    Δ = vec(mean(X, dims=1) .- mean(Y, dims=1))
    S = poolcov!(cov(X), nx - 1, cov(Y), ny - 1) .* (inv(nx) .+ inv(ny))
    T² = At_Binv_A(Δ, S)
    F = T² * (nx + ny - p - 1) / (p * (nx + ny - 2))
    return EqualCovHotellingT2Test(T², F, nx, ny, p, Δ, S)
end

StatsBase.nobs(T::EqualCovHotellingT2Test) = (T.nx, T.ny)
StatsBase.dof(T::EqualCovHotellingT2Test) = (T.p, T.nx + T.ny - T.p - 1)

testname(::EqualCovHotellingT2Test) =
    "Two sample Hotelling's T² test (equal covariance matrices)"
population_param_of_interest(T::EqualCovHotellingT2Test) =
    ("Difference in mean vectors", zeros(eltype(T.Δ), T.p), T.Δ)

## Two sample: unequal covariance

# Store the denominator degrees of freedom in the type, since the computation
# is expensive and we don't want to redo it every time the user calls dof
struct UnequalCovHotellingT2Test <: HotellingT2TestTest
    T²::Real
    F::Real
    nx::Int
    ny::Int
    p::Int
    ν::Int
    Δ::Vector
    S::Matrix
end

"""
    UnequalCovHotellingT2Test(X::AbstractMatrix, Y::AbstractMatrix)

Perform a two sample Hotelling's ``T^2`` test of the hypothesis that the difference in
the mean vectors of `X` and `Y` is zero, without assuming that `X` and `Y` have equal
covariance matrices.
"""
function UnequalCovHotellingT2Test(X::AbstractMatrix, Y::AbstractMatrix)
    p, nx, ny = checkdims(X, Y)
    Δ = vec(mean(X, dims=1) .- mean(Y, dims=1))
    Sx = cov(X) ./ nx
    Sy = cov(Y) ./ ny
    ST = Sx .+ Sy
    T² = At_Binv_A(Δ, ST)
    F = (nx + ny - p - 1) * T² / (p * (nx + ny - 2))
    tmp = Symmetric(ST) \ Δ
    iν = (dot(tmp, Sx * tmp) / T²)^2 / (nx - 1) + (dot(tmp, Sy * tmp) / T²)^2 / (ny - 1)
    ν = trunc(Int, inv(iν))
    return UnequalCovHotellingT2Test(T², F, nx, ny, p, ν, Δ, ST)
end

StatsBase.nobs(T::UnequalCovHotellingT2Test) = (T.nx, T.ny)
StatsBase.dof(T::UnequalCovHotellingT2Test) = (T.p, T.ν)

testname(::UnequalCovHotellingT2Test) =
    "Two sample Hotelling's T² test (unequal covariance matrices)"
population_param_of_interest(T::UnequalCovHotellingT2Test) =
    ("Difference in mean vectors", zeros(eltype(T.Δ), T.p), T.Δ)
