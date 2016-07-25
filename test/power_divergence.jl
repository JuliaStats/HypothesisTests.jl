using HypothesisTests
using StatsBase
using Base.Test

#Example 1 in R
#Agresti (2007) p. 39

d = [[762,484] [327,239] [468,477]]

m = PowerDivergenceTest(d)

@test_approx_eq m.theta0 [0.25523082406125785,0.19670969099133556,0.11593952361049113,0.08935608756107216,0.1935739395970214,0.1491899341788219]
@test_approx_eq m.thetahat [0.2763873775843308,0.1755531374682626,0.11860718171926006,0.08668842945230323,0.16974972796517954,0.17301414581066377]

c = ci(m)
c0 =
[(0.23322451940515054,0.31882480957222203),(0.1323902792890823,0.21799056945615383),(0.0754443235400798,0.16104461370715129),(0.04352557127312297,0.12912586144019447),(0.12658686978599926,0.21218715995307078),(0.12985128763148351,0.21545157779855498)]

[@test_approx_eq c[i][1] c0[i][1] for i in 1:length(c)]
[@test_approx_eq c[i][2] c0[i][2] for i in 1:length(c)]

@test_approx_eq pvalue(m) 2.9535891832117357e-7
@test_approx_eq m.stat 30.070149095754687
@test_approx_eq m.df 2
@test_approx_eq m.n 2757
@test_approx_eq m.residuals [2.198855766015898,-2.504669492560728,0.4113701700566286,-0.46858294710127296,-2.8432397155451494,3.2386734435365825]
@test_approx_eq m.stdresiduals [4.502053521086705,-4.502053521086705,0.6994517329844298,-0.6994517329844298,-5.315945542704929,5.315945542704929]

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

ci(m, method = :bootstrap)
ci(m, method = :bootstrap, tail=:left)
ci(m, method = :bootstrap, tail=:right)

ci(m, method = :gold)
ci(m, method = :gold, tail=:left)
ci(m, method = :gold, tail=:right)

ci(m, method = :quesenberry_hurst)
ci(m, method = :quesenberry_hurst, tail=:left)
ci(m, method = :quesenberry_hurst, tail=:right)

ci(m, method = :sison_glaz)
ci(m, method = :sison_glaz, correct=false)
ci(m, method = :sison_glaz, tail=:left)
ci(m, method = :sison_glaz, tail=:right)

@test_throws ArgumentError ci(m, method=:FOO)
@test_throws ArgumentError ci(m, tail=:fox)




#Example 3 in R

d = [20, 15, 25]
m = PowerDivergenceTest(d)

@test_approx_eq m.theta0 [0.3333333333333333,0.3333333333333333,0.3333333333333333]
@test_approx_eq m.thetahat [0.3333333333333333,0.25,0.4166666666666667]

c = ci(m)
c0 = [(0.04999999999999999,0.5833301356192295),(0.0,0.49999680228589616),(0.13333333333333336,0.6666634689525628)]

[@test_approx_eq c[i][1] c0[i][1] for i in 1:length(c)]
[@test_approx_eq c[i][2] c0[i][2] for i in 1:length(c)]

@test_approx_eq pvalue(m) 0.2865047968601901
@test_approx_eq m.stat 2.5
@test_approx_eq m.df 2
@test_approx_eq m.n 60
@test_approx_eq m.residuals [0.0,-1.118033988749895,1.118033988749895]
@test_approx_eq m.stdresiduals [0.0,-1.3693063937629153,1.3693063937629153]

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

ci(m, method = :bootstrap)
ci(m, method = :bootstrap, tail=:left)
ci(m, method = :bootstrap, tail=:right)

ci(m, method = :gold)
ci(m, method = :gold, tail=:left)
ci(m, method = :gold, tail=:right)

ci(m, method = :quesenberry_hurst)
ci(m, method = :quesenberry_hurst, tail=:left)
ci(m, method = :quesenberry_hurst, tail=:right)

ci(m, method = :sison_glaz)
ci(m, method = :sison_glaz, correct=false)
ci(m, method = :sison_glaz, tail=:left)
ci(m, method = :sison_glaz, tail=:right)

@test_throws ArgumentError ci(m, method=:FOO)
@test_throws ArgumentError ci(m, tail=:fox)

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

# Issue #43 - large numbers caused overflow
PowerDivergenceTest(rand(1_000_000+(1:10), 3, 3))
