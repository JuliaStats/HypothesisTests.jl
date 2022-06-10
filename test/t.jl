using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "T-test" begin

@testset "One sample" begin
	@test pvalue(OneSampleTTest(-5:5)) == 1

	tst = OneSampleTTest(-5:10)
	@test abs(pvalue(tst) - 0.0530) <= 1e-4
	@test abs(pvalue(tst; tail=:left) - 0.9735) <= 1e-4
	@test abs(pvalue(tst; tail=:right) - 0.0265) <= 1e-4
	@test default_tail(tst) == :both
	show(IOBuffer(), tst)

	tst = OneSampleTTest(mean(-5:10), std(-5:10), 16)
	@test abs(pvalue(tst) - 0.0530) <= 1e-4

	@test all(abs.([confint(tst)...;] - [-0.0369, 5.0369]) .<= 1e-4)
	@test all(abs.([confint(tst, level=0.9)...;] - [0.4135, 4.5865]) .<= 1e-4)
	c = confint(tst; tail=:left)
	@test c[1] == -Inf
	@test abs(c[2] - 4.5865) .<= 1e-4
	c = confint(tst; tail=:right)
	@test abs(c[1] - 0.4135) .<= 1e-4
	@test c[2] == Inf
	show(IOBuffer(), tst)

	tst = OneSampleTTest(-10:5)
	@test abs(pvalue(tst) - 0.0530) <= 1e-4
	@test abs(pvalue(tst; tail=:left) - 0.0265) <= 1e-4
	@test abs(pvalue(tst; tail=:right) - 0.9735) <= 1e-4
	@test all(abs.([confint(tst)...] - [-5.0369, 0.0369]) .<= 1e-4)
	@test abs.(confint(tst; tail=:left)[2] - (-0.4135)) .<= 1e-4
	@test abs.(confint(tst; tail=:right)[1] - (-4.5865)) .<= 1e-4
	show(IOBuffer(), tst)
end

@testset "Paired" begin
	@test abs(pvalue(OneSampleTTest([1, 1, 2, 1, 0], [0, 1, 1, 1, 0])) - 0.1778) <= 1e-4
end

@testset "Two sample" begin
	# From http://en.wikipedia.org/w/index.php?title=Student%27s_t-test&oldid=526762741

	a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
	a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

	tst = EqualVarianceTTest(a1, a2)
	@test tst.df == 10
	@test abs(tst.t - 1.959) <= 1e-3
	@test abs(pvalue(tst) - 0.078) <= 1e-3
	@test all(abs.([confint(tst)...] - [-0.0131, 0.2031]) .<= 1e-4)
	@test default_tail(tst) == :both
	show(IOBuffer(), tst)

	n1 = length(a1)
	n2 = length(a2)
	m1 = mean(a1)
	m2 = mean(a2)
	v1 = var(a1)
	v2 = var(a2)

	tst2 = EqualVarianceTTest(n1, n2, m1, m2, v1, v2)
	@test tst.df == tst2.df
	@test tst.t == tst2.t
	@test pvalue(tst) == pvalue(tst2)

	tst = UnequalVarianceTTest(a1, a2)
	@test abs(tst.df - 7.03) <= 0.01
	@test abs(tst.t - 1.959) <= 1e-3
	@test abs(pvalue(tst) - 0.091) <= 1e-3
	@test all(abs.([confint(tst)...] - [-0.0196, 0.2096]) .<= 1e-4)
	@test default_tail(tst) == :both
	show(IOBuffer(), tst)

end
end
