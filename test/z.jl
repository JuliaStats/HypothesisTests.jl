using HypothesisTests, Base.Test
using Distributions

# This is always the null in our tests.
null = Normal(0.0, 1.0)

## ONE SAMPLE T-TEST

# One sample
x = -5:5
@test pvalue(OneSampleZTest(x)) == 1

x = -5:10
m, s, n = mean(x), std(x), length(x)
se = s / sqrt(n)
z = (m - 0) / se

tst = OneSampleZTest(x)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(pvalue(tst; tail=:left), cdf(null, z))
@test_approx_eq(pvalue(tst; tail=:right), ccdf(null, z))
show(IOBuffer(), tst)

tst = OneSampleZTest(m, s, n)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(ci(tst)[1], m + quantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst)[2], m + cquantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst, 0.10)[1], m + quantile(null, 0.10 / 2) * se)
@test_approx_eq(ci(tst, 0.10)[2], m + cquantile(null, 0.10 / 2) * se)
@test_approx_eq(ci(tst; tail=:left)[1], -Inf)
@test_approx_eq(ci(tst; tail=:left)[2], m + cquantile(null, 0.05) * se)
@test_approx_eq(ci(tst; tail=:right)[1], m + quantile(null, 0.05) * se)
@test_approx_eq(ci(tst; tail=:right)[2], Inf)
show(IOBuffer(), tst)

x = -10:5
m, s, n = mean(x), std(x), length(x)
se = s / sqrt(n)
z = (m - 0) / se

tst = OneSampleZTest(x)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(pvalue(tst; tail=:left), cdf(null, z))
@test_approx_eq(pvalue(tst; tail=:right), ccdf(null, z))
show(IOBuffer(), tst)

tst = OneSampleZTest(m, s, n)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(ci(tst)[1], m + quantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst)[2], m + cquantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst, 0.10)[1], m + quantile(null, 0.10 / 2) * se)
@test_approx_eq(ci(tst, 0.10)[2], m + cquantile(null, 0.10 / 2) * se)
@test_approx_eq(ci(tst; tail=:left)[1], -Inf)
@test_approx_eq(ci(tst; tail=:left)[2], m + cquantile(null, 0.05) * se)
@test_approx_eq(ci(tst; tail=:right)[1], m + quantile(null, 0.05) * se)
@test_approx_eq(ci(tst; tail=:right)[2], Inf)
show(IOBuffer(), tst)

# Paired samples
x, y = [1, 1, 2, 1, 0], [0, 1, 1, 1, 0]
m, s, n = mean(x - y), std(x - y), length(x - y)
se = s / sqrt(n)
z = (m - 0) / se
tst = OneSampleZTest(x, y)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))

## TWO SAMPLE Z-TESTS

a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

tst = EqualVarianceZTest(a1, a2)
m1, s1sq, n1 = mean(a1), var(a1), length(a1)
m2, s2sq, n2 = mean(a2), var(a2), length(a2)
xbar = (m1 - m2)
avg_var = (n1 - 1) / (n1 + n2 - 2) * s1sq + (n2 - 1) / (n1 + n2 - 2) * s2sq
se = sqrt(avg_var / n1 + avg_var / n2)
z = xbar / se

@test_approx_eq(tst.z, z)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(pvalue(tst; tail=:left), cdf(null, z))
@test_approx_eq(pvalue(tst; tail=:right), ccdf(null, z))
@test_approx_eq(ci(tst)[1], xbar + quantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst)[2], xbar + cquantile(null, 0.05 / 2) * se)
show(IOBuffer(), tst)

tst = UnequalVarianceZTest(a1, a2)
se = sqrt(s1sq / n1 + s2sq / n2)
z = xbar / se
@test_approx_eq(tst.z, z)
@test_approx_eq(pvalue(tst), 2 * min(cdf(null, z), ccdf(null, z)))
@test_approx_eq(pvalue(tst; tail=:left), cdf(null, z))
@test_approx_eq(pvalue(tst; tail=:right), ccdf(null, z))
@test_approx_eq(ci(tst)[1], xbar + quantile(null, 0.05 / 2) * se)
@test_approx_eq(ci(tst)[2], xbar + cquantile(null, 0.05 / 2) * se)
show(IOBuffer(), tst)
