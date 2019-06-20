# Tests for equality of covariance matrices

export BartlettTest

abstract type CovarianceEqualityTest <: HypothesisTest end

population_param_of_interest(::CovarianceEqualityTest) =
    ("Equality of covariance matrices", NaN, NaN)

## Utility functions

# Finite population correction factor
@inline _correction(p::Int, nxm1::Int, nym1::Int) =
    1 - ((2p^2 + 3p - 1) * (inv(nxm1) + inv(nym1) - inv(nxm1 + nym1)) / 6(p + 1))

## Bartlett's test

struct BartlettTest <: CovarianceEqualityTest
    L′::Real
    p::Int
    nx::Int
    ny::Int
end

"""
    BartlettTest(X::AbstractMatrix, Y::AbstractMatrix)

Perform Bartlett's test of the hypothesis that the covariance matrices of `X` and `Y`
are equal.

!!! note
    Bartlett's test is sensitive to departures from multivariate normality.
"""
function BartlettTest(X::AbstractMatrix, Y::AbstractMatrix)
    nx, p = size(X)
    ny, q = size(Y)
    p == q || throw(DimensionMismatch("Inconsistent number of variables"))
    a = nx - 1
    b = ny - 1
    Sx = cov(X)
    Sy = cov(Y)
    L′ = -a * logdet(Sx) - b * logdet(Sy)
    L′ += (a + b) * logdet(poolcov!(Sx, a, Sy, b))
    L′ *= _correction(p, a, b)
    return BartlettTest(L′, p, nx, ny)
end

StatsBase.nobs(B::BartlettTest) = (B.nx, B.ny)
StatsBase.dof(B::BartlettTest) = div(B.p * (B.p + 1), 2)

testname(::BartlettTest) = "Bartlett's Test for Equality of Covariance Matrices"
default_tail(::BartlettTest) = :right
pvalue(B::BartlettTest; tail=:right) = pvalue(Chisq(dof(B)), B.L′, tail=tail)

function show_params(io::IO, B::BartlettTest, indent="")
    println(io, indent, "number of observations: ", nobs(B))
    println(io, indent, "number of variables:    ", B.p)
    println(io, indent, "χ² statistic:           ", B.L′)
    println(io, indent, "degrees of freedom:     ", dof(B))
end
