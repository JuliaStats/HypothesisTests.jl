using HypothesisTests, Base.Test
using StatsBase
using HypothesisTests: tail

#Example 1 in R
#Agresti (2007) p. 39

d = [[762,484] [327,239] [468,477]]

m = PowerDivergenceTest(d)

@test m.theta0 ≈ [0.25523082406125785,0.19670969099133556,0.11593952361049113,0.08935608756107216,0.1935739395970214,0.1491899341788219]
@test m.thetahat ≈ [0.2763873775843308,0.1755531374682626,0.11860718171926006,0.08668842945230323,0.16974972796517954,0.17301414581066377]

c = confint(m)
c0 = [(0.23322451940515054, 0.31882480957222203),
      (0.1323902792890823, 0.21799056945615383),
      (0.0754443235400798, 0.16104461370715129),
      (0.04352557127312297, 0.12912586144019447),
      (0.12658686978599926, 0.21218715995307078),
      (0.12985128763148351, 0.21545157779855498)]

for i = 1:length(c)
    @test c[i][1] ≈ c0[i][1]
    @test c[i][2] ≈ c0[i][2]
end

@test pvalue(m) ≈ 2.9535891832117357e-7
@test tail(m) == :right
# @test HypothesisTests.default_tail(m) == :right
@test m.stat ≈ 30.070149095754687
@test m.df ≈ 2
@test m.n ≈ 2757
@test m.residuals ≈ reshape([2.198855766015898,-2.504669492560728,0.4113701700566286,-0.46858294710127296,-2.8432397155451494,3.2386734435365825], size(d))
@test m.stdresiduals ≈ reshape([4.502053521086705,-4.502053521086705,0.6994517329844298,-0.6994517329844298,-5.315945542704929,5.315945542704929], size(d))

m = PowerDivergenceTest(d,lambda=0.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-1.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-2.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-0.5)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=2/3)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = ChisqTest(d)
m = MultinomialLRT(d)

confint(m, method = :bootstrap)
confint(m, method = :bootstrap, tail=:left)
confint(m, method = :bootstrap, tail=:right)

confint(m, method = :gold)
confint(m, method = :gold, tail=:left)
confint(m, method = :gold, tail=:right)

confint(m, method = :quesenberry_hurst)
confint(m, method = :quesenberry_hurst, tail=:left)
confint(m, method = :quesenberry_hurst, tail=:right)

confint(m, method = :sison_glaz)
confint(m, method = :sison_glaz, correct=false)
confint(m, method = :sison_glaz, tail=:left)
confint(m, method = :sison_glaz, tail=:right)

@test_throws ArgumentError confint(m, method=:FOO)
@test_throws ArgumentError confint(m, tail=:fox)




#Example 3 in R

d = [ 20, 15, 25 ]
m = PowerDivergenceTest(d)

@test m.theta0 ≈ [0.3333333333333333,0.3333333333333333,0.3333333333333333]
@test m.thetahat ≈ [0.3333333333333333,0.25,0.4166666666666667]

c = confint(m)
c0 = [(0.04999999999999999,0.5833301356192295),(0.0,0.49999680228589616),(0.13333333333333336,0.6666634689525628)]

for i in 1:length(c)
    @test c[i][1] ≈ c0[i][1]
    @test c[i][2] ≈ c0[i][2]
end

@test pvalue(m) ≈ 0.2865047968601901
@test m.stat ≈ 2.5
@test m.df ≈ 2
@test m.n ≈ 60
@test m.residuals ≈ [0.0,-1.118033988749895,1.118033988749895]
@test m.stdresiduals ≈ [0.0,-1.3693063937629153,1.3693063937629153]

m = PowerDivergenceTest(d,lambda=0.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-1.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-2.0)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=-0.5)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = PowerDivergenceTest(d,lambda=2/3)
testname(m)
pvalue(m)
show(IOBuffer(), m)

m = ChisqTest(d)
m = MultinomialLRT(d)

confint(m, method = :bootstrap)
confint(m, method = :bootstrap, tail=:left)
confint(m, method = :bootstrap, tail=:right)

confint(m, method = :gold)
confint(m, method = :gold, tail=:left)
confint(m, method = :gold, tail=:right)

confint(m, method = :quesenberry_hurst)
confint(m, method = :quesenberry_hurst, tail=:left)
confint(m, method = :quesenberry_hurst, tail=:right)

confint(m, method = :sison_glaz)
confint(m, method = :sison_glaz, correct=false)
confint(m, method = :sison_glaz, tail=:left)
confint(m, method = :sison_glaz, tail=:right)

@test_throws ArgumentError confint(m, method=:FOO)
@test_throws ArgumentError confint(m, tail=:fox)

#
x=[1,2,3,1,2,3]
y=[1,1,1,2,2,3]

d = counts(x,y,3)

ChisqTest(d)
MultinomialLRT(d)

PowerDivergenceTest(x,y,3)
PowerDivergenceTest(x,y,(1:3,1:3))

ChisqTest(x,y,3)
ChisqTest(x,y,(1:3,1:3))

MultinomialLRT(x,y,3)
MultinomialLRT(x,y,(1:3,1:3))

# Test that large counts don't cause overflow (issue #43)
d = [113997 1031298
     334453 37471]
PowerDivergenceTest(d)
