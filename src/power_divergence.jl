export PowerDivergenceTest, ChisqTest, MultinomialLRT

typealias Levels{T} @compat(Tuple{UnitRange{T},UnitRange{T}})

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
    return("Pearson's Chi-square Test")
  elseif x.lambda == 0
    return("Multinomial Likelihood Ratio Test")
  elseif x.lambda == -1
    return("Minimum Discrimination Information Test")
  elseif x.lambda == -2
    return("Neyman Modified Chi-square Test")
  elseif x.lambda == -.5
    return("Freeman-Tukey Test")
  else
    return("Power Divergence Test")
  end
end

# parameter of interest: name, value under h0, point estimate
population_param_of_interest(x::PowerDivergenceTest) = ("Multinomial Probabilities", x.theta0, x.thetahat)

pvalue(x::PowerDivergenceTest; tail=:both) = pvalue(Chisq(x.df),x.stat; tail=tail)

function ci(x::PowerDivergenceTest, alpha::Float64=0.05; tail::Symbol=:both, method::Symbol=:sison_glaz, correct::Bool=true, bootstrap_iters::Int64=10000, GC::Bool=true)
  check_alpha(alpha)

  m  = length(x.thetahat)
  
  if tail == :left
    i = ci(x, alpha*2,method=method, GC=GC)
    @compat Tuple{Float64,Float64}[ (0.0, i[j][2]) for j in 1:m]
  elseif tail == :right
    i = ci(x, alpha*2,method=method, GC=GC)
    @compat Tuple{Float64,Float64}[ (i[j][1], 1.0) for j in 1:m]
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
  m = mapslices(x -> quantile(x,[alpha/2,1-alpha/2]), rand(Multinomial(x.n,convert(Vector{Float64},x.thetahat)),iters)/x.n, 2)
  @compat Tuple{Float64,Float64}[ ( boundproportion(m[i,1]), boundproportion(m[i,2])) for i in 1:length(x.thetahat)]
end

# Quesenberry, Hurst (1964)
function ci_quesenberry_hurst(x::PowerDivergenceTest,alpha::Float64; GC::Bool=true)
  m  = length(x.thetahat)
  cv = GC ? quantile(Chisq(1),1-alpha/m) : quantile(Chisq(m-1), 1-alpha)

  u = (cv + 2*x.observed + sqrt(  cv * (cv + 4 * x.observed .* (x.n - x.observed)/x.n)))/( 2*(x.n + cv))
  l = (cv + 2*x.observed - sqrt(  cv * (cv + 4 * x.observed .* (x.n - x.observed)/x.n)))/( 2*(x.n + cv))
  @compat Tuple{Float64,Float64}[ (boundproportion(l[j]), boundproportion(u[j])) for j in 1:m] 
end

# asymptotic simultaneous confidence interval 
# Gold (1963)
function ci_gold(x::PowerDivergenceTest, alpha::Float64; correct::Bool=true, GC::Bool=true)
  m  = length(x.thetahat)
  cv = GC ? quantile(Chisq(1), 1-alpha/(2*m)) : quantile(Chisq(m-1), 1-alpha/2)

  u = [ x.thetahat[j] + sqrt(cv * x.thetahat[j]*(1-x.thetahat[j])/x.n) + (correct ? 1/(2*x.n) : 0 ) for j in 1:m]
  l = [ x.thetahat[j] - sqrt(cv * x.thetahat[j]*(1-x.thetahat[j])/x.n) - (correct ? 1/(2*x.n) : 0 ) for j in 1:m]
  @compat Tuple{Float64,Float64}[ (boundproportion(l[j]), boundproportion(u[j])) for j in 1:m] 
end

# Simultaneous confidence interval
# Sison, Glaz (1995)
# adapted from SAS macro by May, Johnson (2000)
function ci_sison_glaz(x::PowerDivergenceTest, alpha::Float64; skew_correct::Bool=true)
  k = length(x.thetahat)
  probn = 1/pdf( Poisson(x.n),x.n)
  
  c = 0
  p_old = 0

  for c in 1:x.n
    #run truncpoi
    m = zeros(k,5)
    for i in 1:k
      lambda = x.observed[i]
      #run moments
      a = lambda + c
      b = max(lambda - c, 0)
      poislama = lambda > 0 ? cdf(Poisson(lambda),a) : 1
      poislamb = lambda > 0 ? cdf(Poisson(lambda),b-1) : 1
      den = b > 0 ? poislama-poislamb : poislama 

      mu  = zeros(4,1)
      for r in 1:4
        poisA = 0
        poisB = 0
        plar = lambda > 0 ? cdf(Poisson(lambda),a-r) : 1
        plbr = lambda > 0 ? cdf(Poisson(lambda),b-r-1): 1
        poisA = a - r >= 0 ? poislama-plar : poislama
        if  b - r - 1 >= 0
          poisB = poislamb - plbr
        end
        if b - r - 1 < 0 && b-1 >= 0
          poisB = poislamb
        end
        if b - r - 1 < 0 && b-1 < 0
          poisB = 0
        end
        mu[r] = (lambda^r) * (1-(poisA-poisB)/den)
      end
      # end of moments
      m[i,1] = mu[1]
      m[i,2] = mu[2] + mu[1] - mu[1]^2
      m[i,3] = mu[3]+mu[2]*(3-3*mu[1])+(mu[1]-3*mu[1]^2+2*mu[1]^3)
      m[i,4] = mu[4] + mu[3] * (6-4*mu[1]) + mu[2] * (7-12*mu[1]+6*mu[1]^2) + mu[1]-4*mu[1]^2+6*mu[1]^3-3*mu[1]^4
      m[i,5] = den
    end
    [ m[i,4] = m[i,4] - 3*m[i,2]^2 for i in 1:k]
    s1,s2,s3,s4 = mapslices(sum,m,1)
    z  = (x.n-s1)/sqrt(s2)
    g1 = s3/(s2^(3/2))
    g2 = s4*(s2^2)
    
    poly = 1 + g1 * (z^3 - 3*z)/6 + g2*(z^4-6*z^2+3)/24 + g1^2 * (z^6 - 15*z^4 + 45*z^2 -15)/72
    f = poly*exp(-z^2/2)/(sqrt(2*pi))
    probx = 1

    [ probx *= probx*m[i,5] for i in 1:k]

    p = probn*probx*f/sqrt(s2)
    # end of truncpoi

    p > 1-alpha && p_old < 1-alpha ? break : nothing
    p_old = p
  end

  delta = (1-alpha - p_old)/(p_old)
  out = zeros(k,5)
  num = zeros(k,1)

  c = c-1

  vol1 = 1
  vol2 = 1
  for i in 1:k
    num[i,1] = i
    cn = c/x.n
    onen = 1/x.n
    out[i,1] = x.thetahat[i]

    out[i,2] = max(x.thetahat[i]-cn,0)
    out[i,3] = min(x.thetahat[i]+cn+2*delta/x.n,1)

    out[i,4] = max(x.thetahat[i] - cn - onen,0)
    out[i,5] = min(x.thetahat[i] + cn + onen,1)
  end
  if skew_correct
    return( @compat Tuple{Float64,Float64}[ (boundproportion(out[i,2]),boundproportion(out[i,3])) for i in 1:k] )
  else
    return( @compat Tuple{Float64,Float64}[ (boundproportion(out[i,4]),boundproportion(out[i,5])) for i in 1:k] )
  end
end

# power divergence statistic for testing goodness of fit
# Cressie and Read 1984; Read and Cressie 1988
# lambda  =  1: Pearson's Chi-square statistic
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
  any( x .< 0) || any( !isfinite(x)) ? error("all entries must be nonnegative and finite") : nothing
  n == 0 ? error("at least one entry must be positive") : nothing
  isfinite(nrows) && isfinite(ncols) && isfinite(nrows*ncols) ? nothing : error("Invalid number of rows or columns")

  df = thetahat = xhat = V = 0

  if nrows > 1 && ncols > 1
    rowsums = mapslices(sum, x, 2)
    colsums = mapslices(sum, x, 1)
    df = (nrows-1)*(ncols-1)
    thetahat = [ x[i]/n for i in 1:length(x)]
    xhat = rowsums * colsums / n
    theta0 = (xhat/n)[:]
    V = Float64[ colsums[j]*rowsums[i]*(n-rowsums[i])*(n-colsums[j])/n^3 for i in 1:nrows, j in 1:ncols]
  elseif nrows == 1 || ncols == 1
    df = length(x) - 1
    xhat = reshape(n * theta0, size(x))
    thetahat = x/n
    V = reshape(n * theta0 .* (1-theta0), size(x))

    abs( 1 - sum(theta0)) <= sqrt(eps()) ? nothing : error("Probabilities must sum to one") 
  else
    error("Number of rows or columns must be at least 1")
  end

  stat = 0
  if lambda == 0
    for i in 1:length(x)
      stat += x[i]*(log(x[i]) - log(xhat[i])) 
    end
    stat *= 2
  elseif lambda == -1
    for i in 1:length(x)
      stat += xhat[i]*(log(xhat[i]) - log(x[i]))
    end
    stat *= 2
  else
    for i in 1:length(x)
      stat += x[i]*( (x[i]/(xhat[i]))^lambda - 1)
    end
    stat *= 2/(lambda*(lambda+1))
  end
  
  PowerDivergenceTest(lambda, theta0, stat, df, x, n, thetahat[:], xhat, (x - xhat)./sqrt(xhat), (x-xhat)./sqrt(V))
end

#convenience functions

#PDT
function PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T}; lambda::U=1.0)
  d = counts(x,y,levels)
  PowerDivergenceTest(d,lambda=lambda)
end

function PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, y::AbstractVector{T}, k::T; lambda::U=1.0)
  d = counts(x,y,k)
  PowerDivergenceTest(d,lambda=lambda)
end

PowerDivergenceTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}; lambda::U=1.0, theta0::Vector{U} = ones(length(x))/length(x)) = 
  PowerDivergenceTest(reshape(x,length(x),1),lambda=lambda,theta0=theta0)

#ChisqTest
function ChisqTest{T<:Integer}(x::AbstractMatrix{T})
  PowerDivergenceTest(x,lambda=1.0)
end

function ChisqTest{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T})
  d = counts(x,y,levels)
  PowerDivergenceTest(d,lambda=1.0)
end

function ChisqTest{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, k::T)
  d = counts(x,y,k)
  PowerDivergenceTest(d,lambda=1.0)
end

ChisqTest{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) = 
  PowerDivergenceTest(reshape(x,length(x),1), lambda=1.0, theta0 = theta0)

#MultinomialLRT
function MultinomialLRT{T<:Integer}(x::AbstractMatrix{T})
  PowerDivergenceTest(x,lambda=0.0)
end

function MultinomialLRT{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, levels::Levels{T})
  d = counts(x,y,levels)
  PowerDivergenceTest(d,lambda=0.0)
end

function MultinomialLRT{T<:Integer}(x::AbstractVector{T}, y::AbstractVector{T}, k::T)
  d = counts(x,y,k)
  PowerDivergenceTest(d,lambda=0.0)
end

MultinomialLRT{T<:Integer,U<:AbstractFloat}(x::AbstractVector{T}, theta0::Vector{U} = ones(length(x))/length(x)) = 
  PowerDivergenceTest(reshape(x,length(x),1), lambda=0.0, theta0 = theta0)

function show_params(io::IO, x::PowerDivergenceTest, ident="")
  println(io,ident, "Sample size:        $(x.n)")
  println(io,ident, "statistic:          $(x.stat)")
  println(io,ident, "degrees of freedom: $(x.df)")
  println(io,ident, "residuals:          $(x.residuals[:])")
  println(io,ident, "std. residuals:     $(x.stdresiduals[:])")
end

