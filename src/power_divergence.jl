export PowerDivergenceTest, ChisqTest, MultinomialLRTest

const Levels{T} = Tuple{UnitRange{T},UnitRange{T}}

function boundproportion(x::T) where T<:Real
    max(min(convert(Float64,x),1.0),0.0)
end

struct PowerDivergenceTest <: HypothesisTest
    lambda::Float64
    theta0::Vector{Float64}
    stat::Float64
    df::Int64
    observed::Matrix{Int64}
    n::Int64
    thetahat::Vector{Float64}

    expected::Matrix{Float64}
    residuals::Matrix{Float64}
    stdresiduals::Matrix{Float64}
end

#test name
function testname(x::PowerDivergenceTest)
    if x.lambda == 1
        return "Pearson's Chi-square Test"
    elseif x.lambda == 0
        return "Multinomial Likelihood Ratio Test"
    elseif x.lambda == -1
        return "Minimum Discrimination Information Test"
    elseif x.lambda == -2
        return "Neyman Modified Chi-square Test"
    elseif x.lambda == -0.5
        return "Freeman-Tukey Test"
    else
        return "Power Divergence Test"
    end
end

# parameter of interest: name, value under h0, point estimate
population_param_of_interest(x::PowerDivergenceTest) = ("Multinomial Probabilities", x.theta0, x.thetahat)
default_tail(test::PowerDivergenceTest) = :right

pvalue(x::PowerDivergenceTest; tail=:right) = pvalue(Chisq(x.df),x.stat; tail=tail)

"""
    confint(test::PowerDivergenceTest; level = 0.95, tail = :both, method = :auto)

Compute a confidence interval with coverage `level` for multinomial proportions using
one of the following methods. Possible values for `method` are:

  - `:auto` (default): If the minimum of the expected cell counts exceeds 100, Quesenberry-Hurst
    intervals are used, otherwise Sison-Glaz.
  - `:sison_glaz`: Sison-Glaz intervals
  - `:bootstrap`: Bootstrap intervals
  - `:quesenberry_hurst`: Quesenberry-Hurst intervals
  - `:gold`: Gold intervals (asymptotic simultaneous intervals)

# References

  * Agresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013.
  * Sison, C.P and Glaz, J. Simultaneous confidence intervals and sample size determination
    for multinomial proportions. Journal of the American Statistical Association,
    90:366-369, 1995.
  * Quesensberry, C.P. and Hurst, D.C. Large Sample Simultaneous Confidence Intervals for
    Multinational Proportions. Technometrics, 6:191-195, 1964.
  * Gold, R. Z. Tests Auxiliary to ``χ^2`` Tests in a Markov Chain. Annals of
    Mathematical Statistics, 30:56-74, 1963.
"""
function StatsBase.confint(x::PowerDivergenceTest; level::Float64=0.95,
                           tail::Symbol=:both, method::Symbol=:auto, correct::Bool=true,
                           bootstrap_iters::Int64=10000, GC::Bool=true)
    check_level(level)

    m  = length(x.thetahat)

    if tail == :left
        i = StatsBase.confint(x, level=1-(1-level)*2, method=method, GC=GC)
        Tuple{Float64,Float64}[(0.0, i[j][2]) for j in 1:m]
    elseif tail == :right
        i = StatsBase.confint(x, level=1-(1-level)*2, method=method, GC=GC)
        Tuple{Float64,Float64}[(i[j][1], 1.0) for j in 1:m]
    elseif tail == :both
        if method == :auto
            method = minimum(x.expected) > 100 ? :quesenberry_hurst : :sison_glaz
        end
        if method == :gold
            ci_gold(x, 1-level, correct=correct, GC=GC)
        elseif method == :sison_glaz
            ci_sison_glaz(x, 1-level, skew_correct=correct)
        elseif method == :quesenberry_hurst
            ci_quesenberry_hurst(x, 1-level, GC=GC)
        elseif method == :bootstrap
            ci_bootstrap(x, 1-level, bootstrap_iters)
        else
            throw(ArgumentError("method=$(method) is invalid or not implemented yet"))
        end
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# Bootstrap
function ci_bootstrap(x::PowerDivergenceTest,alpha::Float64, iters::Int64)
    m = mapslices(x -> quantile(x, [alpha / 2, 1 - alpha / 2]), rand(Multinomial(x.n, convert(Vector{Float64}, x.thetahat)),iters) / x.n, dims=2)
    Tuple{Float64,Float64}[(boundproportion(m[i,1]), boundproportion(m[i,2])) for i in 1:length(x.thetahat)]
end

# Quesenberry, Hurst (1964)
function ci_quesenberry_hurst(x::PowerDivergenceTest, alpha::Float64; GC::Bool=true)
    m  = length(x.thetahat)
    cv = GC ? quantile(Chisq(1), 1 - alpha / m) : quantile(Chisq(m - 1), 1 - alpha)

    u = (cv .+ 2 .* x.observed .+ sqrt.(cv .* (cv .+ 4 .* x.observed .* (x.n .- x.observed) ./ x.n))) ./ (2 .* (x.n .+ cv))
    l = (cv .+ 2 .* x.observed .- sqrt.(cv .* (cv .+ 4 .* x.observed .* (x.n .- x.observed) ./ x.n))) ./ (2 .* (x.n .+ cv))
    Tuple{Float64,Float64}[(boundproportion(l[j]), boundproportion(u[j])) for j in 1:m]
end

# asymptotic simultaneous confidence interval
# Gold (1963)
function ci_gold(x::PowerDivergenceTest, alpha::Float64; correct::Bool=true, GC::Bool=true)
    m  = length(x.thetahat)
    cv = GC ? quantile(Chisq(1), 1 - alpha / 2m) : quantile(Chisq(m - 1), 1 - alpha / 2)

    u = [x.thetahat[j] + sqrt.(cv * x.thetahat[j] * (1 - x.thetahat[j]) / x.n) + ifelse(correct, inv(2x.n), 0) for j in 1:m]
    l = [x.thetahat[j] - sqrt.(cv * x.thetahat[j] * (1 - x.thetahat[j]) / x.n) - ifelse(correct, inv(2x.n), 0) for j in 1:m]
    Tuple{Float64,Float64}[ (boundproportion(l[j]), boundproportion(u[j])) for j in 1:m]
end

# Simultaneous confidence interval
# Sison, Glaz (1995)
# adapted from SAS macro by May, Johnson (2000)
function ci_sison_glaz(x::PowerDivergenceTest, alpha::Float64; skew_correct::Bool=true)
    k = length(x.thetahat)
    probn = inv(pdf(Poisson(x.n), x.n))

    c = 0
    p = 0.0
    p_old = 0.0
    m1, m2, m3, m4, m5 = zeros(k), zeros(k), zeros(k), zeros(k), zeros(k)
    mu = zeros(4)

    for _c in 1:x.n
        #run truncpoi
        for i in 1:k
            lambda = x.observed[i]
            #run moments
            a = lambda + _c
            b = max(lambda - _c, 0)
            poislama = cdf(Poisson(lambda), a)
            poislamb = cdf(Poisson(lambda), b - 1)
            den = b > 0.0 ? poislama-poislamb : poislama

            for r in 1:4
                plar = cdf(Poisson(lambda), a - r)
                plbr = cdf(Poisson(lambda), b - r - 1)

                poisA = ifelse( (a - r) >= 0, poislama - plar, poislama)
                poisB = 0.0
                if  (b - r - 1) >= 0
                    poisB = poislamb - plbr
                end
                if (b - r - 1) < 0 && b - 1 >= 0
                    poisB = poislamb
                end
                if (b - r - 1) < 0 && (b - 1) < 0
                    poisB = 0.0
                end
                mu[r] = lambda^r * (1 - (poisA - poisB) / den)
            end
            # end of moments
            m1[i] = mu[1]
            m2[i] = mu[2] + mu[1] - mu[1]^2
            m3[i] = mu[3] + mu[2] * (3 - 3mu[1]) + (mu[1] - 3mu[1]^2 + 2mu[1]^3)
            m4[i] = mu[4] + mu[3] * (6 - 4mu[1]) + mu[2] * (7 - 12mu[1] + 6mu[1]^2) + mu[1] - 4mu[1]^2 + 6mu[1]^3 - 3mu[1]^4
            m5[i] = den
        end
        for i in 1:k
            m4[i] -= 3m2[i]^2
        end
        s1, s2, s3, s4 = sum(m1), sum(m2), sum(m3), sum(m4)
        z  = (x.n - s1) / sqrt(s2)
        g1 = s3 / s2^(3/2)
        g2 = s4 / s2^2

        poly = 1 + g1 * (z^3 - 3z) / 6 + g2 * (z^4 - 6z^2 + 3) / 24 + g1^2 * (z^6 - 15z^4 + 45z^2 - 15) / 72
        f = poly * exp(-z^2 / 2) / sqrt(2π)
        probx = prod(m5)

        p = probn * probx * f / sqrt(s2)
        # end of truncpoi

        if p > 1 - alpha && p_old < 1 - alpha
            c = _c
            break
        else
            c = x.n
        end
        p_old = p
    end

    delta = (1 - alpha - p_old) / (p - p_old)
    out = zeros(k, 5)
    num = zeros(k, 1)

    c -= 1

    vol1 = 1
    vol2 = 1
    for i in 1:k
        num[i,1] = i
        cn = c / x.n
        onen = 1 / x.n
        out[i,1] = x.thetahat[i]

        out[i,2] = max(x.thetahat[i] - cn, 0)
        out[i,3] = min(x.thetahat[i] + cn + 2delta / x.n, 1)

        out[i,4] = max(x.thetahat[i] - cn - onen, 0)
        out[i,5] = min(x.thetahat[i] + cn + onen, 1)
    end
    if skew_correct
        return Tuple{Float64,Float64}[(boundproportion(out[i,2]), boundproportion(out[i,3])) for i in 1:k]
    else
        return Tuple{Float64,Float64}[(boundproportion(out[i,4]), boundproportion(out[i,5])) for i in 1:k]
    end
end

# power divergence statistic for testing goodness of fit
# Cressie and Read 1984; Read and Cressie 1988
# lambda    =  1: Pearson's Chi-square statistic
# lambda ->  0: Converges to Likelihood Ratio test stat
# lambda -> -1: Converges to minimum discrimination information statistic (Gokhale, Kullback 1978)
# lambda =  -2: Neyman Modified chi-squared statistic (Neyman 1949)
# lambda = -.5: Freeman-Tukey statistic (Freeman, Tukey 1950)

# Under regularity conditions, their asymptotic distributions are all the same (Drost 1989)
# Chi-squared null approximation works best for lambda near 2/3
"""
    PowerDivergenceTest(x[, y]; lambda = 1.0, theta0 = ones(length(x))/length(x))

Perform a Power Divergence test.

If `y` is not given and `x` is a matrix with one row or column, or `x` is a vector, then
a goodness-of-fit test is performed (`x` is treated as a one-dimensional contingency
table). In this case, the hypothesis tested is whether the population probabilities equal
those in `theta0`, or are all equal if `theta0` is not given.

If `x` is a matrix with at least two rows and columns, it is taken as a two-dimensional
contingency table. Otherwise, `x` and `y` must be vectors of the same length. The contingency
table is calculated using the `counts` function from the `StatsBase` package. Then the power
divergence test is conducted under the null hypothesis that the joint distribution of the
cell counts in a 2-dimensional contingency table is the product of the row and column
marginals.

Note that the entries of `x` (and `y` if provided) must be non-negative integers.

Computed confidence intervals by default are Quesenberry-Hurst intervals
if the minimum of the expected cell counts exceeds 100, and Sison-Glaz intervals otherwise.
See the [`confint(::PowerDivergenceTest)`](@ref) documentation for a list of
supported methods to compute confidence intervals.

The power divergence test is given by
```math
    \\dfrac{2}{λ(λ+1)}\\sum_{i=1}^I \\sum_{j=1}^J n_{ij} \\left[(n_{ij}
    /\\hat{n}_{ij})^λ -1\\right]
```
where ``n_{ij}`` is the cell count in the ``i`` th row and ``j`` th column and ``λ`` is a
real number determining the nature of the test to be performed:

  * ``λ = 1``: equal to Pearson's chi-squared statistic
  * ``λ \\to 0``: converges to the likelihood ratio test statistic
  * ``λ \\to -1``: converges to the minimum discrimination information statistic
    (Gokhale and Kullback, 1978)
  * ``λ = -2``: equals Neyman modified chi-squared (Neyman, 1949)
  * ``λ = -1/2``: equals the Freeman-Tukey statistic (Freeman and Tukey, 1950).

Under regularity conditions, the asymptotic distributions are identical (see Drost et. al.
1989). The ``χ^2`` null approximation works best for ``λ`` near ``2/3``.

Implements: [`pvalue`](@ref), [`confint(::PowerDivergenceTest)`](@ref)

# References

  * Agresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013.
"""
function PowerDivergenceTest(x::AbstractMatrix{T}; lambda::U=1.0, theta0::Vector{U} = ones(length(x))/length(x)) where {T<:Integer,U<:AbstractFloat}

    nrows, ncols = size(x)
    n = sum(x)

    #validate date
    (any(x .< 0) || any(!isfinite, x)) && throw(ArgumentError("all entries must be nonnegative and finite"))
    n == 0 && throw(ArgumentError("at least one entry must be positive"))
    (!isfinite(nrows) || !isfinite(ncols) || !isfinite(nrows*ncols)) && throw(ArgumentError("invalid number of rows or columns"))

    if nrows > 1 && ncols > 1
        rowsums = sum(x, dims=2)
        colsums = sum(x, dims=1)
        df = (nrows - 1) * (ncols - 1)
        thetahat = x ./ n
        xhat = rowsums * colsums / n
        theta0 = xhat / n
        V = Float64[(colsums[j]/n) * (rowsums[i]/n) * (1 - rowsums[i]/n) * (n - colsums[j]) for i in 1:nrows, j in 1:ncols]
    elseif nrows == 1 || ncols == 1
        df = length(x) - 1
        xhat = reshape(n * theta0, size(x))
        thetahat = x / n
        V = reshape(n .* theta0 .* (1 .- theta0), size(x))

        abs(1 - sum(theta0)) <= sqrt(eps()) || throw(ArgumentError("Probabilities must sum to one"))
    else
        throw(ArgumentError("Number of rows or columns must be at least 1"))
    end

    stat = 0
    if lambda == 0
        for i in 1:length(x)
            stat += x[i] * (log(x[i]) - log(xhat[i]))
        end
        stat *= 2
    elseif lambda == -1
        for i in 1:length(x)
            stat += xhat[i] * (log(xhat[i]) - log(x[i]))
        end
        stat *= 2
    else
        for i in 1:length(x)
            stat += x[i] * ((x[i] / xhat[i])^lambda - 1)
        end
        stat *= 2 / (lambda * (lambda + 1))
    end

    PowerDivergenceTest(lambda, vec(theta0), stat, df, x, n, vec(thetahat), xhat, (x - xhat) ./ sqrt.(xhat), (x - xhat) ./ sqrt.(V))
end

#convenience functions

#PDT
function PowerDivergenceTest(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T}; lambda::U=1.0) where {T<:Integer,U<:AbstractFloat}
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=lambda)
end

function PowerDivergenceTest(x::AbstractVector{T}, y::AbstractVector{T}, k::T; lambda::U=1.0) where {T<:Integer,U<:AbstractFloat}
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=lambda)
end

PowerDivergenceTest(x::AbstractVector{T}; lambda::U=1.0, theta0::Vector{U} = ones(length(x))/length(x)) where {T<:Integer,U<:AbstractFloat} =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=lambda, theta0=theta0)

#ChisqTest
"""
    ChisqTest(x[, y][, theta0 = ones(length(x))/length(x)])

Perform a Pearson chi-squared test (equivalent to a [`PowerDivergenceTest`](@ref)
with ``λ = 1``).

If `y` is not given and `x` is a matrix with one row or column, or `x` is a vector, then
a goodness-of-fit test is performed (`x` is treated as a one-dimensional contingency
table). In this case, the hypothesis tested is whether the population probabilities equal
those in `theta0`, or are all equal if `theta0` is not given.

If only `y` and `x` are given and both are vectors of integer type, then once again a
goodness-of-fit test is performed. In this case, `theta0` is calculated by the proportion
of each individual values in `y`. Here, the hypothesis tested is whether the two samples
`x` and `y` come from the same population or not.

If `x` is a matrix with at least two rows and columns, it is taken as a two-dimensional
contingency table. Otherwise, `x` and `y` must be vectors of the same length. The contingency
table is calculated using `counts` function from the `StatsBase` package. Then the power
divergence test is conducted under the null hypothesis that the joint distribution of the
cell counts in a 2-dimensional contingency table is the product of the row and column
marginals.

Note that the entries of `x` (and `y` if provided) must be non-negative integers.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function ChisqTest(x::AbstractMatrix{T}) where T<:Integer
    PowerDivergenceTest(x, lambda=1.0)
end

function ChisqTest(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T}) where T<:Integer
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=1.0)
end

function ChisqTest(x::AbstractVector{T}, y::AbstractVector{T}, k::T) where T<:Integer
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=1.0)
end

function ChisqTest(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Integer}
    theta0 = y ./ sum(y)
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=1.0, theta0=theta0)
end

ChisqTest(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) where {T<:Integer,U<:AbstractFloat} =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=1.0, theta0=theta0)

#MultinomialLRTest
"""
    MultinomialLRTest(x[, y][, theta0 = ones(length(x))/length(x)])

Perform a multinomial likelihood ratio test (equivalent to a [`PowerDivergenceTest`](@ref)
with ``λ = 0``).

If `y` is not given and `x` is a matrix with one row or column, or `x` is a vector, then
a goodness-of-fit test is performed (`x` is treated as a one-dimensional contingency
table). In this case, the hypothesis tested is whether the population probabilities equal
those in `theta0`, or are all equal if `theta0` is not given.

If `x` is a matrix with at least two rows and columns, it is taken as a two-dimensional
contingency table. Otherwise, `x` and `y` must be vectors of the same length. The contingency
table is calculated using `counts` function from the `StatsBase` package. Then the power
divergence test is conducted under the null hypothesis that the joint distribution of the
cell counts in a 2-dimensional contingency table is the product of the row and column
marginals.

Note that the entries of `x` (and `y` if provided) must be non-negative integers.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function MultinomialLRTest(x::AbstractMatrix{T}) where T<:Integer
    PowerDivergenceTest(x, lambda=0.0)
end

function MultinomialLRTest(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T}) where T<:Integer
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=0.0)
end

function MultinomialLRTest(x::AbstractVector{T}, y::AbstractVector{T}, k::T) where T<:Integer
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=0.0)
end

MultinomialLRTest(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) where {T<:Integer,U<:AbstractFloat} =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=0.0, theta0=theta0)

function show_params(io::IO, x::PowerDivergenceTest, ident="")
    println(io, ident, "Sample size:        $(x.n)")
    println(io, ident, "statistic:          $(x.stat)")
    println(io, ident, "degrees of freedom: $(x.df)")
    print(io, ident, "residuals:          ")
    show(io, vec(x.residuals))
    println(io)
    print(io, ident, "std. residuals:     ")
    show(io, vec(x.stdresiduals))
    println(io)
end
