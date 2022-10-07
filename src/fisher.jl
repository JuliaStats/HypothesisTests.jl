# fisher.jl
# Fisher's exact test
#
# Copyright (C) 2013   Simon Kornblith
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

export FisherExactTest

"""
    FisherExactTest(a::Integer, b::Integer, c::Integer, d::Integer)

Perform Fisher's exact test of the null hypothesis that the success probabilities ``a/c``
and ``b/d`` are equal, that is the odds ratio ``(a/c) / (b/d)`` is one, against the
alternative hypothesis that they are not equal.

See [`pvalue(::FisherExactTest)`](@ref) and [`confint(::FisherExactTest)`](@ref) for details
about the computation of the default p-value and confidence interval, respectively.

The contingency table is structured as:

| -  | X1 | X2 |
|:--:|:--:|:--:|
|*Y1*| a  | b  |
|*Y2*| c  | d  |

!!! note
    The `show` function output contains the conditional maximum likelihood estimate of the
    odds ratio rather than the sample odds ratio; it maximizes the likelihood given by
    Fisher's non-central hypergeometric distribution.

Implements: [`pvalue(::FisherExactTest)`](@ref), [`confint(::FisherExactTest)`](@ref)

# References

  * Fay, M.P., Supplementary material to "Confidence intervals that match Fisher’s exact or
    Blaker’s exact tests". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374,
    [link](https://doi.org/10.1093/biostatistics/kxp050)
"""
struct FisherExactTest <: HypothesisTest
    # Format:
    # X1  X2
    # Y1  a  b
    # Y2  c  d
    a::Int
    b::Int
    c::Int
    d::Int

    # conditional maximum likehood estimate of odd ratio
    ω::Float64

    function FisherExactTest(a::Int, b::Int, c::Int, d::Int)
        ω = cond_mle_odds_ratio(a, b, c, d)
        new(a, b, c, d, ω)
    end
end

testname(::FisherExactTest) = "Fisher's exact test"
population_param_of_interest(x::FisherExactTest) = ("Odds ratio", 1.0, x.ω) # parameter of interest: name, value under h0, point estimate
default_tail(test::FisherExactTest) = :both

# The sizing argument to print_matrix was removed during the 0.5 dev period
function _print_matrix(io::IO, X::AbstractVecOrMat, pre::AbstractString)
    Base.print_matrix(io, X, pre)
end

function show_params(io::IO, x::FisherExactTest, ident="")
    println(io, ident, "contingency table:")
    _print_matrix(io, [x.a x.b; x.c x.d], repeat(ident, 2))
    println(io)
end

# DOC: for tail=:both there exist multiple ``method``s for computing a pvalue and the corresponding ci.
"""
    pvalue(x::FisherExactTest; tail = :both, method = :central)

Compute the p-value for a given Fisher exact test.

The one-sided p-values are based on Fisher's non-central hypergeometric distribution
``f_ω(i)`` with odds ratio ``ω``:
```math
    \\begin{align*}
        p_ω^{(\\text{left})} &=\\sum_{i ≤ a} f_ω(i)\\\\
        p_ω^{(\\text{right})} &=\\sum_{i ≥ a} f_ω(i)
    \\end{align*}
```
For `tail = :both`, possible values for `method` are:

  - `:central` (default): Central interval, i.e. the p-value is two times the minimum of the
    one-sided p-values.
  - `:minlike`: Minimum likelihood interval, i.e. the p-value is computed by summing all
    tables with the same marginals that are equally or less probable:
    ```math
        p_ω = \\sum_{f_ω(i)≤ f_ω(a)} f_ω(i)
    ```

# References

  * Gibbons, J.D., Pratt, J.W., P-values: Interpretation and Methodology, American
    Statistican, 29(1):20-25, 1975.
  * Fay, M.P., Supplementary material to "Confidence intervals that match Fisher’s exact or
    Blaker’s exact tests". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374,
    [link](https://doi.org/10.1093/biostatistics/kxp050)
"""
function pvalue(x::FisherExactTest; tail=:both, method=:central)
    if tail == :both && method != :central
        if method == :minlike
            p = pvalue_both_minlike(x)
        else
            throw(ArgumentError("method=$(method) is not implemented yet"))
        end
    else
        p = pvalue(Hypergeometric(x.a + x.b, x.c + x.d, x.a + x.c), x.a, tail=tail)
    end
    p = max(min(p, 1.0), 0.0)

    return p
end

function pvalue_both_minlike(x::FisherExactTest, ω::Float64=1.0)
    a, b, c, d = reorder(x.a, x.b, x.c, x.d)
    if a == c == 0 || b == d == 0
        return 1.0
    end
    dist = FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)

    p = pdf(dist, a)
    v = nextfloat(p)
    if a != 0
        p += cdf(dist, a-1)
    end

    # Add p-values of all tables in other tail equally or less probable
    for i = a+c:-1:a+1
        curp = pdf(dist, i)
        if curp > v
            break
        end
        p += curp
    end
    p
end

# confidence interval by inversion of p-value
"""
    confint(x::FisherExactTest; level::Float64=0.95, tail=:both, method=:central)

Compute a confidence interval with coverage `level`. One-sided intervals are based on
Fisher's non-central hypergeometric distribution. For `tail = :both`, the only
`method` implemented yet is the central interval (`:central`).

!!! note
    Since the p-value is not necessarily unimodal, the corresponding confidence region might
    not be an interval.

# References

  * Gibbons, J.D, Pratt, J.W. P-values: Interpretation and Methodology, American
    Statistican, 29(1):20-25, 1975.
  * Fay, M.P., Supplementary material to "Confidence intervals that match Fisher’s exact or
    Blaker’s exact tests". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374,
    [link](https://doi.org/10.1093/biostatistics/kxp050)
"""
function StatsBase.confint(x::FisherExactTest; level::Float64=0.95, tail=:both, method=:central)
    check_level(level)
    if x.a == x.c == 0 || x.b == x.d == 0
        return (0.0, Inf)
    end
    dist(ω) = FisherNoncentralHypergeometric(x.a+x.b, x.c+x.d, x.a+x.c, ω)
    obj(ω) = pvalue(dist(ω), x.a, tail=tail) - (1-level)

    if tail == :left # upper bound
        if (x.a == maximum(dist(1.0)))
            (0.0, Inf)
        else
            lower, upper = find_brackets(obj)
            (0.0, lower == upper ? lower : find_zero(obj, (lower, upper)))
        end
    elseif tail == :right # lower bound
        if (x.a == minimum(dist(1.0)))
            (0.0, Inf)
        else
            lower, upper = find_brackets(obj)
            (lower == upper ? lower : find_zero(obj, (lower, upper)), Inf)
        end
    elseif tail == :both
        if method == :central
            (StatsBase.confint(x, level=1-(1-level)/2, tail=:right)[1],
             StatsBase.confint(x, level=1-(1-level)/2, tail=:left)[2])
        else
            throw(ArgumentError("method=$(method) is not implemented yet"))
        end
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

## helpers

function reorder(a,b,c,d)
    if a + c > b + d
        a, b, c, d = b, a, d, c
    end
    if a/c > b/d
        a, b, c, d = c, d, a, b
    end
    (a, b, c, d)
end

# find values x_lower and x_upper, s.t. f(x_lower) < 0 and f(x_upper) > 0 or vice versa
function find_brackets(f::Function, x_init::Float64=1.0)
    f_init = f(x_init)

    if f_init > f(x_init + 1.0)
        find_brackets(x -> -f(x), x_init)
    else
        x_upper = x_lower = x_init
        if f_init > 0.0
            while x_lower > eps(0.0) && f(x_lower) > 0.0
                x_lower /= 2
            end
        else
            while f(x_upper) < 0.0
                x_upper *= 2
            end
        end
        (x_lower, x_upper)
    end
end

# find odds ratio ω that maximizes the (conditional) likelihood;
# since the mode and mean of Fisher's Noncentral Hypergeometric distribution
# coincide, this is equivalent to find ω, s.t., mean(dist(ω)) = a
function cond_mle_odds_ratio(a::Int, b::Int, c::Int, d::Int)
    if a == c == 0 || b == d == 0
        return 0.0
    end
    dist(ω) = FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)

    if (a == minimum(dist(1.0)))
        0.0
    elseif a == maximum(dist(1.0))
        Inf
    else
        obj(ω) = mean(dist(ω))-a
        lower, upper = find_brackets(obj)
        lower == upper ? lower : find_zero(obj, (lower, upper))
    end
end
