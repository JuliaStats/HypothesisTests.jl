# twobinomial.jl
# Two sample hypothesis tests for proportions
#
# Copyright (C) 2020 Vladimir Arnautov
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

export TwoSampleBinomialTest

struct TwoSampleBinomialTask
    x1::Int
    n1::Int
    x2::Int
    n2::Int
    δ::Float64
    ptype::Symbol
    htype::Symbol
    alpha::Float64
    method::Symbol
end

struct TwoSampleBinomialTest <: HypothesisTest
    test::TwoSampleBinomialTask
    est
    se
    ci
    result
end

function cihalpha(htype, alpha)
    if  htype == :equivalence || htype == :superiority
        cialpha = alpha * 2
    elseif htype == :equality
        cialpha = alpha
    else
        throw(ArgumentError("unknown hypothesis=$(htype)"))
    end
    cialpha
end

"""
    TwoSampleBinomialTest(x1::Real, n1::Real, x2::Real, n2::Real, δ::Real; ptype::Symbol, htype::Symbol, alpha::Real = 0.05, method::Symbol = :default)

Perform a two sample binomial hypothesis check via computing confidence intrval
with coverage `1 - alpha` for `equality` test or `1 - 2 * alpha` for `equivalence`
or `superiority` test.

Hypothesis type `htype` are:
- `:equality`
- `:equivalence`
- `:superiority`

Parameters type `ptype` are:
    - `:difference` : risk difference;
    - `:riskratio` : risk ratio;
    - `:oddratio` : odd ratio;

For each parameters can be used  following methods.

Possible values for `method` for `difference` are:

  - `:mn`
  - `:wald`
  - `:waldcc`
  - `:nhs`
  - `:nhscc`
  - `:ac`
  - `:ha`
  - `:mnmee`
  - `:mover`

Possible values for `method` for `riskratio` are:

  - `:mn`
  - `:cli`
  - `:li`
  - `:mover`

Possible values for `method` for `oddratio` are:

  - `:mn`
  - `:woolf`
  - `:awoolf`
  - `:mover`

"""
function TwoSampleBinomialTest(x1::Real, n1::Real, x2::Real, n2::Real, δ::Real;
                               ptype::Symbol, htype::Symbol, alpha::Real = 0.05, method::Symbol = :default)
    test = TwoSampleBinomialTask(x1, n1, x2, n2, δ, ptype, htype, alpha, method)
    cialpha = cihalpha(htype, alpha)
    if ptype == :difference
        est, se, lci, uci = ci_diff(test; level = 1 - cialpha, method = method)
    elseif ptype == :oddratio
        est, se, lci, uci = ci_or(test; level = 1 - cialpha, method = method)
    elseif ptype == :riskratio
        est, se, lci, uci = ci_rr(test; level = 1 - cialpha, method = method)
    else
        throw(ArgumentError("unknown parameter=$(ptype)"))
    end
    result = false
    if htype == :equivalence
        if lci > -abs(δ) && uci < abs(δ) result = true end
    elseif htype == :equality
        if lci > δ || uci < δ result = true end
    elseif htype == :superiority
        if lci > δ result = true end
    end
    TwoSampleBinomialTest(test, est, se, (lci, uci), result)
end
################################################################################
function ci_prop_wilson(x::Int, n::Int, alpha::Real)
    z   = abs(quantile(Normal(), 1 - alpha / 2))
    p   = x / n
    d   = 1 + (z^2) / n
    se  = sqrt((p * (1 - p) + (z^2) / (4 * n)) / n) / d
    est = (p + (z^2) / (2 * n)) / d
    return est, se, est - z * se, est + z * se
end
function ci_prop_wilson_cc(x::Int, n::Int, alpha::Real)
    z = abs(quantile(Normal(), 1 - alpha / 2))
    p = x / n
    l = (2*n*p+z*z-1-z*sqrt(z*z-2-1/n+4*p*(n*(1-p)+1)))/2/(n+z*z)
    u = (2*n*p+z*z+1+z*sqrt(z*z+2-1/n+4*p*(n*(1-p)-1)))/2/(n+z*z)
    return p, NaN, min(p, l), max(p, u)
end
################################################################################
function ci_diff(test; level, method)
    if method == :mn || method == :default
        ci_diff_mn(test, level)
    elseif method == :wald
        ci_diff_wald(test, level)
    elseif method == :waldcc
        ci_diff_wald_cc(test, level)
    elseif method == :nhs
        ci_diff_nhs(test, level)
    elseif method == :nhscc
        ci_diff_nhs_cc(test, level)
    elseif method == :ac
        ci_diff_ac(test, level)
    elseif method == :ha
        ci_diff_ha(test, level)
    elseif method == :mnmee
        ci_diff_mnmee(test, level)
    elseif method == :mover
        ci_diff_mover(test, level)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
function ci_rr(test; level, method)
    if method == :mn || method == :default
        ci_rr_mn(test, level)
    elseif method == :cli
        ci_rr_cli(test, level)
    elseif method == :li
        ci_rr_li(test, level)
    elseif method == :mover
        ci_rr_mover(test, level)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
function ci_or(test; level, method)
    if method == :mn || method == :default
        ci_or_mn(test, level)
    elseif method == :woolf
        ci_or_woolf(test, level)
    elseif method == :awoolf
        ci_or_awoolf(test, level)
    elseif method == :mover
        ci_or_mover(test, level)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
################################################################################
#Not implementsd yet
function pvalue(tr::TwoSampleBinomialTest)
    if tr.test.method ∈ [:wald, :waldcc, :ac, :ha]
        #return ccdf(Normal(), z)*2.
        throw(ArgumentError("not implemented"))
    elseif tr.test.method ∈ [:mn, :nhs, :nhscc, :mnmee]
        #return NaN
        throw(ArgumentError("not implemented"))
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
function confint(tr::TwoSampleBinomialTest)
    throw(ArgumentError("not implemented"))
end
################################################################################
#Method of Mee 1984 with Miettinen and Nurminen modification n / (n - 1) Newcombe 1998
#Score intervals for the difference of two binomial proportions
@inline function mle_diff(p1::Float64, n1::Int, p2::Float64, n2::Int, δ::Float64)
    if p1 - p2 - δ == 0 return 0.0 end
    θ = n2 / n1
    a = 1 + θ
    b = -(1 + θ + p1 + θ * p2 + δ * (θ + 2))
    c = δ^2 + δ * (2 * p1 + θ + 1) + p1 + θ * p2
    d = -p1 * δ * (1 + δ)
    v = (b / a / 3)^3 - b * c/(6 * a * a) + d / 2 / a
    u = sign(v) * sqrt((b / 3 / a)^2 - c / 3 / a)
    w = (pi + acos(v / u^3)) / 3
    p1n = 2 * u * cos(w) - b / 3 /a
    p2n = p1n - δ
    return p1n, p2n
end
@inline function mn_diff_z_val(p1::Real, n1::Int, p2::Real, n2::Int, est::Real, δ::Real)
    p1n, p2n = mle_diff(p1, n1, p2, n2, δ)
    return (est - δ)^2 / ((n1 + n2) / (n1 + n2 - 1) * (p1n * (1 - p1n) / n1 + p2n * (1 - p2n) / n2))
end
function ci_diff_mn(test, level; atol::Float64 = 1E-8)
    ests, ses, lcis, ucis = ci_diff_nhs_cc(test, level)
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Chisq(1), level)
    fmnd(x)  = mn_diff_z_val(p1, test.n1, p2, test.n2, est, x) - z

    if fmnd(lcis) * fmnd(est - eps()) < 0.0
        ll = lcis
        lu = est - eps()
    else
        ll = -1.0 + eps()
        lu = lcis
    end
    if fmnd(ucis) * fmnd(est + eps()) < 0.0
        ul = est + eps()
        uu = ucis
    else
        ul = ucis
        uu = 1.0 - eps()
    end
    lci = find_zero(fmnd, (ll, lu), atol = atol)
    uci = find_zero(fmnd, (ul, uu), atol = atol)
    return est, NaN, lci, uci
end
#Wald
function ci_diff_wald(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / test.n1 + p2 * (1 - p2) / test.n2)
    return est, se, est - z * se, est + z * se
end
#Wald continuity correction
function ci_diff_wald_cc(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / test.n1 + p2 * (1 - p2) / test.n2)
    cc       = 0.5*(1 / test.n1 + 1 / test.n2)
    return est, se, est - z * se - cc, est + z * se + cc
end
#Newcombes Hybrid (wilson) Score interval for the difference of proportions
#Newcombe 1998
function ci_diff_nhs(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    est1, se1, lci1, uci1 = ci_prop_wilson(test.x1, test.n1, alpha)
    est2, se2, lci2, uci2 = ci_prop_wilson(test.x2, test.n2, alpha)
    return est, NaN, est - z * sqrt(lci1 * (1 - lci1)/test.n1 + uci2 * (1 - uci2) / test.n2), est + z * sqrt(uci1 * (1 - uci1) / test.n1 + lci2 * (1 - lci2) / test.n2)
end
#Newcombes Hybrid Score continuity correction
function ci_diff_nhs_cc(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    est1, se1, lci1, uci1 = ci_prop_wilson_cc(test.x1, test.n1, alpha)
    est2, se2, lci2, uci2 = ci_prop_wilson_cc(test.x2, test.n2, alpha)
    return est, NaN, est - sqrt((p1 - lci1)^2 + (uci2 - p2)^2), est + sqrt((uci1 - p1)^2 + (p2 - lci2)^2)
end
#Agresti-Caffo interval for the difference of proportions
#Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000)
function ci_diff_ac(test, level)
    alpha    = 1 - level
    z        = quantile(Normal(), 1 - alpha / 2)
    p1I      = (test.x1 + 1) / (test.n1 + 2)
    p2I      = (test.x2 + 1) / (test.n2 + 2)
    n1I      = test.n1 + 2
    n2I      = test.n2 + 2
    est      = p1I - p2I
    se       = sqrt(p1I * (1 - p1I) / n1I + p2I * (1 - p2I) / n2I)
    return est, se, est - z * se, est + z * se
end
function ci_diff_ha(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / (test.n1 - 1) + p2 * (1 - p2) / (test.n2 - 1))
    cc       = 1 / min(test.n1, test.n2)
    return est, se, est - z * se - cc, est + z * se + cc
end
@inline function mnmee_z_val(p1::Real, n1::Int, p2::Real, n2::Int, est::Real, δ::Real)
    p1n, p2n = mle_diff(p1, n1, p2, n2, δ)
    return (est - δ)^2 / (p1n * (1 - p1n) / n1 + p2n * (1 - p2n) / n2)
end
#Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176
#MN - no correction
function ci_diff_mnmee(test, level; atol::Float64 = 1E-8)
    ests, ses, lcis, ucis = ci_diff_nhs_cc(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    z        = quantile(Chisq(1), 1 - alpha)
    fmnd(x)  = mnmee_z_val(p1, test.n1, p2, test.n2, est, x) - z
    if fmnd(lcis) * fmnd(est - eps()) < 0.0
        ll = lcis
        lu = est - eps()
    else
        ll = -1.0 + eps()
        lu = lcis
    end
    if fmnd(ucis) * fmnd(est + eps()) < 0.0
        ul = est + eps()
        uu = ucis
    else
        ul = ucis
        uu = 1.0 - eps()
    end
    lci = find_zero(fmnd, (ll, lu), atol = atol)
    uci = find_zero(fmnd, (ul, uu), atol = atol)
    return est, NaN, lci, uci
end
#Brown, Li's Jeffreys
function ci_diff_jeffreys(test, level)
    p1   = (test.x1 + 0.5) / (test.n1 + 1)
    p2   = (test.x2 + 0.5) / (test.n2+1)
    se   = sqrt(p1*(1 - p1) / test.n1 + p2 * (1 - p2) / test.n2)
    z    = quantile(Normal(), 1 - alpha / 2)
    est  = p1 - p2
    return test.x1 / test.n1 - test.x2 / test.n2, se, max(-1.0, est - z * se), min(1.0, est + z * se)
end
#Method of variance estimates recovery
function ci_diff_mover(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 - p2
    Z        = quantile(Normal(), 1 - alpha / 2)
    est1, se1, lci1, uci1   = ci_prop_wilson(test.x1, test.n1, alpha)
    est2, se2, lci2, uci2   = ci_prop_wilson(test.x2, test.n2, alpha)
    lci      = est - sqrt((p1 - lci1)^2 + (uci2 - p2)^2)
    uci      = est + sqrt((uci1 - p1)^2 + (p2 - lci2)^2)
    return est, NaN, lci, uci
end
################################################################################
@inline function mle_or(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int)
    a  = n2 * (φ-1)
    b  = φ * n1 + n2 - (x1 + x2) * (φ - 1)
    c  = -(x1 + x2)
    p2 = (-b + sqrt(b * b - 4 * a * c)) / a / 2
    p1 = p2 * φ/(1 + p2 * (φ - 1))
    return p1, p2
end
@inline function mle_or_z_val(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int)
    p1 = x1 / n1
    pmle1, pmle2 = mle_or(φ, x1, n1, x2, n2)
    return (n1 * (p1 - pmle1))^2 * (1 / (n1 * pmle1 * (1 - pmle1)) + 1/(n2 * pmle2 * (1 - pmle2))) / ((n1 + n2)/(n1 + n2 - 1))
end
################################################################################
#Miettinen O. S., Nurminen M. (1985) Comparative analysis of two rates.Statistics in Medicine4,213–226
#MN Score
function ci_or_mn(test, level; atol::Float64 = 1E-8)
    z        = quantile(Chisq(1), level)
    fmnor(x) = mle_or_z_val(x, test.x1, test.n1, test.x2, test.n2) - z
    if (test.x1==0 && test.x2==0) || (test.x1==test.n1 && test.x2==test.n2)
        return NaN, NaN, 0.0, Inf
    elseif test.x1==0 || test.x2==test.n2
        return 0.0, NaN, 0.0, find_zero(fmnor, 1e-6, atol = atol)
    elseif test.x1==test.n1 || test.x2 == 0
        return Inf, NaN, find_zero(fmnor, 1e-6, atol = atol), Inf
    else
        est  = (test.x1 / (test.n1 - test.x1)) / (test.x2 / (test.n2 - test.x2))
        return est, NaN, find_zero(fmnor, 1e-6, atol = atol), find_zero(fmnor, est+1e-6, atol = atol)
    end
end
#Woolf logit
#Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253.
function ci_or_woolf(test, level)
    alpha     = 1 - level
    xa        = test.x1
    xb        = test.n1 - test.x1
    xc        = test.x2
    xd        = test.n2 - test.x2
    est       = xa*xd/xc/xb
    estI      = log(est)
    se        = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
    z         = quantile(Normal(), 1 - alpha / 2)
    return est, se, exp(estI - z * se), exp(estI + z * se)
end
#Adjusted Woolf interval (Gart adjusted logit) Lawson, R (2005):Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation andComputation, 33, 1095-1113.
#Gart, J. J. (1966). Alternative analyses of contingency tables. Journal of the Royal Statistical Society. Series B (Methodological), 28:164-179.
function ci_or_awoolf(test, level)
    alpha     = 1 - level
    xa        = test.x1 + 0.5
    xb        = test.n1 - test.x1 + 0.5
    xc        = test.x2 + 0.5
    xd        = test.n2 - test.x2 + 0.5
    est       = xa*xd/xc/xb
    estI      = log(est)
    se        = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
    z         = quantile(Normal(), 1 - alpha/2)
    return est, se, exp(estI - z * se), exp(estI + z * se)
end
#Method of variance estimates recovery
#Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.
function ci_or_mover(test, level)
    alpha    = 1 - level
    p1       = (test.x1/(test.n1-test.x1))
    p2       = (test.x2/(test.n2-test.x2))
    est      = p1/p2
    z        = quantile(Normal(), 1-alpha/2)
    est1, se1, lci1, uci1   = ci_prop_wilson(test.x1, test.n1, alpha)
    est2, se2, lci2, uci2   = ci_prop_wilson(test.x2, test.n2, alpha)
    vl1      = lci1/(1 - lci1)
    vu1      = uci1/(1 - uci1)
    vl2      = lci2/(1 - lci2)
    vu2      = uci2/(1 - uci2)
    lci      = (p1*p2-sqrt((p1*p2)^2 - vl1*vu2*(2*p1-vl1)*(2*p2-vu2)))/(vu2*(2*p2 - vu2))
    uci      = (p1*p2+sqrt((p1*p2)^2 - vu1*vl2*(2*p1-vu1)*(2*p2-vl2)))/(vl2*(2*p2 - vl2))
    return est, NaN, lci, uci
end
################################################################################
@inline function mle_rr(φ, x1::Int, n1::Int, x2::Int, n2::Int)
    a = (n1 + n2) * φ
    b = -(φ * (x1 + n2) + x2 + n1)
    c = x1 + x2
    p2 = (-b - sqrt(b * b - 4 * a * c)) / 2 / a
    p1 = p2 * φ
    return p1, p2
end
@inline function mle_rr_z_val(φ, x1::Int, n1::Int, x2::Int, n2::Int)
    p1 = x1 / n1
    p2 = x2 / n2
    pmle1, pmle2 = mle_rr(φ, x1, n1, x2, n2)
    return (p1 - φ * p2)^2 / ((pmle1 * (1 - pmle1) / n1 + φ * φ * pmle2 * (1 - pmle2) / n2) * ((n1 + n2 - 1) / (n1 + n2)))
end
################################################################################
#Miettinen-Nurminen Score interval
#Miettinen, O. and Nurminen, M. (1985), Comparative analysis of two rates. Statist. Med., 4: 213-226. doi:10.1002/sim.4780040211
function ci_rr_mn(test, level; atol::Float64 = 1E-8)
    z        = quantile(Chisq(1), level)
    fmnrr(x) = mle_rr_z_val(x, test.x1, test.n1, test.x2, test.n2) - z
    if (test.x1==0 && test.x2==0) || (test.x1==test.n1 && test.x2==test.n2)
        return NaN, NaN, 0.0, Inf
    elseif test.x1==0 || test.x2==test.n2
        return 0.0, NaN, 0.0, find_zero(fmnrr, 1e-8, atol=atol)
    elseif test.x1==test.n1 || test.x2 == 0
        return Inf, NaN, find_zero(fmnrr, 1e-8, atol=atol), Inf
    else
        est = (test.x1 / test.n1) / (test.x2 / test.n2)
        return est, NaN, find_zero(fmnrr, 1e-8, atol=atol), find_zero(fmnrr, est+1e-6, atol=atol)
    end
end
#Crude log interval
#Gart, JJand Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338.
function ci_rr_cli(test, level)
    alpha     = 1 - level
    x1I       = test.x1 + 0.5
    x2I       = test.x2 + 0.5
    n1I       = test.n1 + 0.5
    n2I       = test.n2 + 0.5
    estI      = log((x1I / n1I) / (x2I / n2I))
    se        = sqrt(1 / x2I + 1 / x1I - 1 / n2I - 1 / n1I)
    est       = (test.x1 / test.n1) / (test.x2 / test.n2)
    z         =  quantile(Normal(), 1 - alpha / 2)
    return est, se, exp(estI - z * se), exp(estI + z * se)
end
#Katz D, Baptista J, Azen SP and Pike MC. Obtaining confidence intervals for the risk ratio in cohort studies. Biometrics 1978; 34: 469–474
function ci_rr_li(test, level)
    alpha     = 1 - level
    est       = (test.x1 / test.n1) / ( test.x2 / test.n2)
    estI      = log(est)
    se        = sqrt(1 / test.x2 + 1 / test.x1 - 1 / test.n2 - 1 / test.n1)
    z         = quantile(Normal(), 1 - alpha / 2)
    return est, se, exp(estI - z * se), exp(estI + z * se)
end
#Method of variance estimates recovery (Donner, Zou, 2012)
function ci_rr_mover(test, level)
    alpha    = 1 - level
    p1       = test.x1 / test.n1
    p2       = test.x2 / test.n2
    est      = p1 / p2
    Z        = quantile(Normal(), 1 - alpha / 2)
    est1, se1, lci1, uci1   = ci_prop_wilson(test.x1, test.n1, alpha)
    est2, se2, lci2, uci2   = ci_prop_wilson(test.x2, test.n2, alpha)
    lci      = (p1 * p2 - sqrt((p1 * p2)^2 - lci1 * uci2 * (2 * p1 - lci1) * (2 * p2 - uci2))) / (uci2 * (2 * p2 - uci2))
    uci      = (p1 * p2 + sqrt((p1 * p2)^2 - uci1 * lci2 * (2 * p1 - uci1) * (2 * p2 - lci2))) / (lci2 * (2 * p2 - lci2))
    return est, NaN, lci, uci
end
################################################################################

function Base.show(_io::IO, res::TwoSampleBinomialTest)
    println(_io, "Two sample binomial hypothesis test")
    println(_io, "Parameter type: ", res.test.ptype)
    println(_io, "Hypothesis type: ", res.test.htype)
    println(_io, "A: ", res.test.x1, " / ", res.test.n1)
    println(_io, "B: ", res.test.x2, " / ", res.test.n2)
    println(_io, "Estimate: ", res.est)
    println(_io, "SE: ", res.se)

    cialpha = cihalpha(res.test.htype, res.test.alpha)

    println(_io, "$((1 - cialpha)*100)% CI: ", res.ci)
    println(_io, "H₀ rejected: ", res.result)
end

################################################################################

#function StatsBase.confint(x::TwoSampleBinomialTest; level::Float64=0.95, tail=:both, method=:default)
#end
