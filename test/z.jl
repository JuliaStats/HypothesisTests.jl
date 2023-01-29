using HypothesisTests, Test
using Distributions
using HypothesisTests: default_tail

@testset "Z test" begin
# This is always the null in our tests.
null = Normal(0.0, 1.0)

@testset "One sample" begin

	x = -5:5
	@test pvalue(OneSampleZTest(x)) == 1

	x = -5:10
	m, s, n = mean(x), std(x), length(x)
	se = s / sqrt(n)
	z = (m - 0) / se

	tst = OneSampleZTest(x)
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test pvalue(tst; tail=:left) ≈ cdf(null, z)
	@test pvalue(tst; tail=:right) ≈ ccdf(null, z)
	@test default_tail(tst) == :both
	show(IOBuffer(), tst)

	tst = OneSampleZTest(m, s, n)
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test confint(tst)[1] ≈ m + quantile(null, 0.05 / 2) * se
	@test confint(tst)[2] ≈ m + cquantile(null, 0.05 / 2) * se
	@test confint(tst; level=0.9)[1] ≈ m + quantile(null, 0.10 / 2) * se
	@test confint(tst; level=0.9)[2] ≈ m + cquantile(null, 0.10 / 2) * se
	@test confint(tst; tail=:left)[1] ≈ -Inf
	@test confint(tst; tail=:left)[2] ≈ m + cquantile(null, 0.05) * se
	@test confint(tst; tail=:right)[1] ≈ m + quantile(null, 0.05) * se
	@test confint(tst; tail=:right)[2] ≈ Inf
	@test_throws ArgumentError confint(tst; tail=2)
	show(IOBuffer(), tst)

	x = -10:5
	m, s, n = mean(x), std(x), length(x)
	se = s / sqrt(n)
	z = (m - 0) / se

	tst = OneSampleZTest(x)
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test pvalue(tst; tail=:left) ≈ cdf(null, z)
	@test pvalue(tst; tail=:right) ≈ ccdf(null, z)
	show(IOBuffer(), tst)

	tst = OneSampleZTest(m, s, n)
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test confint(tst)[1] ≈ m + quantile(null, 0.05 / 2) * se
	@test confint(tst)[2] ≈ m + cquantile(null, 0.05 / 2) * se
	@test confint(tst; level=0.9)[1] ≈ m + quantile(null, 0.10 / 2) * se
	@test confint(tst; level=0.9)[2] ≈ m + cquantile(null, 0.10 / 2) * se
	@test confint(tst; tail=:left)[1] ≈ -Inf
	@test confint(tst; tail=:left)[2] ≈ m + cquantile(null, 0.05) * se
	@test confint(tst; tail=:right)[1] ≈ m + quantile(null, 0.05) * se
	@test confint(tst; tail=:right)[2] ≈ Inf
	show(IOBuffer(), tst)
end

@testset "Paired samples" begin
	x, y = [1, 1, 2, 1, 0], [0, 1, 1, 1, 0]
	m, s, n = mean(x - y), std(x - y), length(x - y)
	se = s / sqrt(n)
	z = (m - 0) / se
	tst = OneSampleZTest(x, y)
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
end

@testset "Two sample" begin
	a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
	a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

	tst = EqualVarianceZTest(a1, a2)
	m1, s1sq, n1 = mean(a1), var(a1), length(a1)
	m2, s2sq, n2 = mean(a2), var(a2), length(a2)
	xbar = (m1 - m2)
	avg_var = (n1 - 1) / (n1 + n2 - 2) * s1sq + (n2 - 1) / (n1 + n2 - 2) * s2sq
	se = sqrt(avg_var / n1 + avg_var / n2)
	z = xbar / se

	@test tst.z ≈ z
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test pvalue(tst; tail=:left) ≈ cdf(null, z)
	@test pvalue(tst; tail=:right) ≈ ccdf(null, z)
	@test default_tail(tst) == :both
	@test confint(tst)[1] ≈ xbar + quantile(null, 0.05 / 2) * se
	@test confint(tst)[2] ≈ xbar + cquantile(null, 0.05 / 2) * se
	show(IOBuffer(), tst)

	tst = UnequalVarianceZTest(a1, a2)
	se = sqrt(s1sq / n1 + s2sq / n2)
	z = xbar / se
	@test tst.z ≈ z
	@test pvalue(tst) ≈ 2 * min(cdf(null, z), ccdf(null, z))
	@test pvalue(tst; tail=:left) ≈ cdf(null, z)
	@test pvalue(tst; tail=:right) ≈ ccdf(null, z)
	@test confint(tst)[1] ≈ xbar + quantile(null, 0.05 / 2) * se
	@test confint(tst)[2] ≈ xbar + cquantile(null, 0.05 / 2) * se
	show(IOBuffer(), tst)
end
end
