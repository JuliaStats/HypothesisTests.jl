export PowerDivergenceTest, ChisqTest, MultinomialLRT

@compat const Levels{T} = Tuple{UnitRange{T},UnitRange{T}}

function boundproportion{T<:Real}(x::T)
    max(min(convert(Float64,x),1.0),0.0)
end


immutable PowerDivergenceTest <: HypothesisTest
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

function StatsBase.confint(x::PowerDivergenceTest, alpha::Float64=0.05;
                           tail::Symbol=:both, method::Symbol=:sison_glaz, correct::Bool=true,
                           bootstrap_iters::Int64=10000, GC::Bool=true)
    check_alpha(alpha)

    m  = length(x.thetahat)

    if tail == :left
        i = StatsBase.confint(x, alpha*2,method=method, GC=GC)
        Tuple{Float64,Float64}[(0.0, i[j][2]) for j in 1:m]
    elseif tail == :right
        i = StatsBase.confint(x, alpha*2,method=method, GC=GC)
        Tuple{Float64,Float64}[(i[j][1], 1.0) for j in 1:m]
    elseif tail == :both
        if method == :gold
            ci_gold(x,alpha,correct=correct,GC=GC)
        elseif method == :sison_glaz
            ci_sison_glaz(x,alpha, skew_correct=correct)
        elseif method == :quesenberry_hurst
            ci_quesenberry_hurst(x,alpha,GC=GC)
        elseif method == :bootstrap
            ci_bootstrap(x,alpha,bootstrap_iters)
        else
            throw(ArgumentError("method=$(method) is invalid or not implemented yet"))
        end
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# Bootstrap
function ci_bootstrap(x::PowerDivergenceTest,alpha::Float64, iters::Int64)
    m = mapslices(x -> quantile(x, [alpha / 2, 1 - alpha / 2]), rand(Multinomial(x.n, convert(Vector{Float64}, x.thetahat)),iters) / x.n, 2)
    Tuple{Float64,Float64}[(boundproportion(m[i,1]), boundproportion(m[i,2])) for i in 1:length(x.thetahat)]
end

# Quesenberry, Hurst (1964)
function ci_quesenberry_hurst(x::PowerDivergenceTest,alpha::Float64; GC::Bool=true)
    m  = length(x.thetahat)
    cv = GC ? quantile(Chisq(1), 1 - alpha / m) : quantile(Chisq(m - 1), 1 - alpha)

    u = (cv + 2*x.observed + sqrt.(cv * (cv + 4 * x.observed .* (x.n - x.observed) / x.n))) / (2*(x.n + cv))
    l = (cv + 2*x.observed - sqrt.(cv * (cv + 4 * x.observed .* (x.n - x.observed) / x.n))) / (2*(x.n + cv))
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
    p_old = 0.0
    m1, m2, m3, m4, m5 = zeros(k), zeros(k), zeros(k), zeros(k), zeros(k)
    mu = zeros(4)

    for c in 1:x.n
        #run truncpoi
        for i in 1:k
            lambda = x.observed[i]
            #run moments
            a = lambda + c
            b = max(lambda - c, 0)
            if lambda > 0
                poislama = cdf(Poisson(lambda), a)
                poislamb = cdf(Poisson(lambda), b - 1)
            else
                poislama = poislamb = 1.0
            end
            den = b > 0 ? poislama-poislamb : poislama

            for r in 1:4
                if lambda > 0
                    plar = cdf(Poisson(lambda), a - r)
                    plbr = cdf(Poisson(lambda), b - r - 1)
                else
                    plar = plbr = 1.0
                end
                poisA = ifelse(a - r >= 0, poislama - plar, poislama)
                poisB = 0.0
                if  b - r - 1 >= 0
                    poisB = poislamb - plbr
                end
                if b - r - 1 < 0 && b - 1 >= 0
                    poisB = poislamb
                end
                if b - r - 1 < 0 && b - 1 < 0
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
        g2 = s4 * s2^2

        poly = 1 + g1 * (z^3 - 3z) / 6 + g2 * (z^4 - 6z^2 + 3) / 24 + g1^2 * (z^6 - 15z^4 + 45z^2 - 15) / 72
        f = poly * exp(-z^2 / 2) / sqrt(2π)
        probx = 1.0
        for i in 1:k
            probx *= probx * m5[i]
        end

        p = probn * probx * f / sqrt(s2)
        # end of truncpoi

        p > 1 - alpha && p_old < 1 - alpha && break
        p_old = p
    end

    delta = (1 - alpha - p_old) / p_old
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
function PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractMatrix{T}; lambda::U=1.0, theta0::Vector{U} = ones(length(x))/length(x))

    nrows, ncols = size(x)
    n = sum(x)

    #validate date
    (any(x .< 0) || any(!isfinite, x)) && throw(ArgumentError("all entries must be nonnegative and finite"))
    n == 0 && throw(ArgumentError("at least one entry must be positive"))
    (!isfinite(nrows) || !isfinite(ncols) || !isfinite(nrows*ncols)) && throw(ArgumentError("invalid number of rows or columns"))

    if nrows > 1 && ncols > 1
        rowsums = sum(x, 2)
        colsums = sum(x, 1)
        df = (nrows - 1) * (ncols - 1)
        thetahat = x ./ n
        xhat = rowsums * colsums / n
        theta0 = xhat / n
        V = Float64[(colsums[j]/n) * (rowsums[i]/n) * (1 - rowsums[i]/n) * (n - colsums[j]) for i in 1:nrows, j in 1:ncols]
    elseif nrows == 1 || ncols == 1
        df = length(x) - 1
        xhat = reshape(n * theta0, size(x))
        thetahat = x / n
        V = reshape(n * theta0 .* (1 - theta0), size(x))

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
function PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T}; lambda::U=1.0)
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=lambda)
end

function PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, y::AbstractVector{T}, k::T; lambda::U=1.0)
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=lambda)
end

PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}; lambda::U=1.0, theta0::Vector{U} = ones(length(x))/length(x)) =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=lambda, theta0=theta0)

#ChisqTest
function ChisqTest{T<:Integer}(x::AbstractMatrix{T})
    PowerDivergenceTest(x, lambda=1.0)
end

function ChisqTest{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T})
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=1.0)
end

function ChisqTest{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, k::T)
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=1.0)
end

ChisqTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=1.0, theta0=theta0)

#MultinomialLRT
function MultinomialLRT{T<:Integer}(x::AbstractMatrix{T})
    PowerDivergenceTest(x, lambda=0.0)
end

function MultinomialLRT{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T})
    d = counts(x, y, levels)
    PowerDivergenceTest(d, lambda=0.0)
end

function MultinomialLRT{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, k::T)
    d = counts(x, y, k)
    PowerDivergenceTest(d, lambda=0.0)
end

MultinomialLRT{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) =
    PowerDivergenceTest(reshape(x, length(x), 1), lambda=0.0, theta0=theta0)

function show_params(io::IO, x::PowerDivergenceTest, ident="")
    println(io, ident, "Sample size:        $(x.n)")
    println(io, ident, "statistic:          $(x.stat)")
    println(io, ident, "degrees of freedom: $(x.df)")
    println(io, ident, "residuals:          $(vec(x.residuals))")
    println(io, ident, "std. residuals:     $(vec(x.stdresiduals))")
end
