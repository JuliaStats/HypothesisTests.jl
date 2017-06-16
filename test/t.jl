using HypothesisTests, Base.Test
using HypothesisTests: tail

iobuffer = IOBuffer() # define once to ease temporary redirection to STDOUT

## ONE SAMPLE T-TEST

# One sample
@test pvalue(OneSampleTTest(-5:5)) == 1

tst = OneSampleTTest(-5:10)
@test abs(pvalue(tst) - 0.0530) <= 1e-4
@test abs(pvalue(tst; tail=:left) - 0.9735) <= 1e-4
@test abs(pvalue(tst; tail=:right) - 0.0265) <= 1e-4
@test tail(tst) == :both
show(iobuffer, tst)

tst = OneSampleTTest(-5:10, tail=:left) # specify the tail for the test
@test abs(pvalue(tst; tail=:both) - 0.0530) <= 1e-4
@test abs(pvalue(tst) - 0.9735) <= 1e-4 # now tail=:left is the default in pvalue(tst)
@test abs(pvalue(tst; tail=:right) - 0.0265) <= 1e-4
@test tail(tst) == :left
show(iobuffer, tst)

tst = OneSampleTTest(-5:10, tail=:right) # specify the tail for the test
@test abs(pvalue(tst; tail=:both) - 0.0530) <= 1e-4
@test abs(pvalue(tst; tail=:left) - 0.9735) <= 1e-4 
@test abs(pvalue(tst) - 0.0265) <= 1e-4 # now tail=:right is the default in pvalue(tst)
@test tail(tst) == :right
show(iobuffer, tst)

tst = OneSampleTTest(mean(-5:10), std(-5:10), 16)
@test abs(pvalue(tst) - 0.0530) <= 1e-4

@test all(abs.([confint(tst)...;] - [-0.0369, 5.0369]) .<= 1e-4)
@test all(abs.([confint(tst, 0.1)...;] - [0.4135, 4.5865]) .<= 1e-4)
c = confint(tst; tail=:left)
@test c[1] == -Inf
@test abs(c[2] - 4.5865) .<= 1e-4
c = confint(tst; tail=:right)
@test abs(c[1] - 0.4135) .<= 1e-4
@test c[2] == Inf
show(iobuffer, tst)

tst = OneSampleTTest(-10:5)
@test abs(pvalue(tst) - 0.0530) <= 1e-4
@test abs(pvalue(tst; tail=:left) - 0.0265) <= 1e-4
@test abs(pvalue(tst; tail=:right) - 0.9735) <= 1e-4
@test all(abs.([confint(tst)...] - [-5.0369, 0.0369]) .<= 1e-4)
@test abs.(confint(tst; tail=:left)[2] - (-0.4135)) .<= 1e-4
@test abs.(confint(tst; tail=:right)[1] - (-4.5865)) .<= 1e-4
show(iobuffer, tst)

tst = OneSampleTTest(-10:5, tail=:left)
@test abs.(confint(tst)[2] - (-0.4135)) .<= 1e-4
@test abs.(confint(tst; tail=:right)[1] - (-4.5865)) .<= 1e-4
show(iobuffer, tst)

tst = OneSampleTTest(-10:5, tail=:right)
@test abs.(confint(tst; tail=:left)[2] - (-0.4135)) .<= 1e-4
@test abs.(confint(tst)[1] - (-4.5865)) .<= 1e-4
show(iobuffer, tst)

# Paired samples
@test abs(pvalue(OneSampleTTest([1, 1, 2, 1, 0], [0, 1, 1, 1, 0])) - 0.1778) <= 1e-4

## TWO SAMPLE T-TESTS

# From http://en.wikipedia.org/w/index.php?title=Student%27s_t-test&oldid=526762741

a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

tst = EqualVarianceTTest(a1, a2)
@test tst.df == 10
@test abs(tst.t - 1.959) <= 1e-3
@test abs(pvalue(tst) - 0.078) <= 1e-3
@test all(abs.([confint(tst)...] - [-0.0131, 0.2031]) .<= 1e-4)
@test tail(tst) == :both
show(iobuffer, tst)

tst = EqualVarianceTTest(a1, a2, tail=:left)
@test abs(pvalue(tst) - 0.960) <= 1e-3
@test abs(pvalue(tst, tail=:both) - 0.078) <= 1e-3
@test tail(tst) == :left
show(iobuffer, tst)

tst = UnequalVarianceTTest(a1, a2)
@test abs(tst.df - 7.03) <= 0.01
@test abs(tst.t - 1.959) <= 1e-3
@test abs(pvalue(tst) - 0.091) <= 1e-3
@test all(abs.([confint(tst)...] - [-0.0196, 0.2096]) .<= 1e-4)
@test tail(tst) == :both
show(iobuffer, tst)

tst = UnequalVarianceTTest(a1, a2, tail=:right)
@test abs(pvalue(tst) - 0.045) <= 1e-3
@test abs(pvalue(tst, tail=:both) - 0.091) <= 1e-3
@test tail(tst) == :right
show(iobuffer, tst)
