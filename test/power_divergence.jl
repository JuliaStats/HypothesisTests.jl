using HypothesisTests, Test
using StatsBase
using HypothesisTests: default_tail

@testset "Power Divergence" begin
#Example 1 in R
#Agresti (2007) p. 39

d = [[762,484] [327,239] [468,477]]

m = PowerDivergenceTest(d)

@test m.theta0 ≈ [0.25523082406125785,0.19670969099133556,0.11593952361049113,0.08935608756107216,0.1935739395970214,0.1491899341788219]
@test m.thetahat ≈ [0.2763873775843308,0.1755531374682626,0.11860718171926006,0.08668842945230323,0.16974972796517954,0.17301414581066377]

c = confint(m, method=:sison_glaz)
c0 = [(.25788900979, .29554669435),
      (.15705476968, .19471245423),
      (.10010881393, .13776649848),
      (.06819006166, .10584774622),
      (.15125136017, .18890904473),
      (.15451577802, .19217346257)]

for i = 1:length(c)
    @test c[i][1] ≈ c0[i][1]
    @test c[i][2] ≈ c0[i][2]
end

@test pvalue(m) ≈ 2.9535891832117357e-7
@test default_tail(m) == :right
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
m = MultinomialLRTest(d)

confint(m)
confint(m, tail=:left)
confint(m, tail=:right)

confint(m, method = :auto)
confint(m, method = :auto, tail=:left)
confint(m, method = :auto, tail=:right)

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

@test confint(m, method = :quesenberry_hurst) == confint(m, method = :auto) == confint(m)


#Example 3 in R

d = [ 20, 15, 25 ]
m = PowerDivergenceTest(d)

@test m.theta0 ≈ [0.3333333333333333,0.3333333333333333,0.3333333333333333]
@test m.thetahat ≈ [0.3333333333333333,0.25,0.4166666666666667]

c = confint(m)
c0 = [(.21666666667, .48098082062),
      (.13333333333, .39764748728),
      (.30000000000, .56431415395)]

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
m = MultinomialLRTest(d)

confint(m)
confint(m, tail=:left)
confint(m, tail=:right)

confint(m, method = :auto)
confint(m, method = :auto, tail=:left)
confint(m, method = :auto, tail=:right)

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

@test confint(m, method = :sison_glaz) == confint(m, method = :auto) == confint(m)

#
x=[1,2,3,1,2,3]
y=[1,1,1,2,2,3]

d = counts(x,y,3)

ChisqTest(d)
MultinomialLRTest(d)

PowerDivergenceTest(x,y,3)
PowerDivergenceTest(x,y,(1:3,1:3))

ChisqTest(x,y,3)
ChisqTest(x,y,(1:3,1:3))

MultinomialLRTest(x,y,3)
MultinomialLRTest(x,y,(1:3,1:3))

# Test that large counts don't cause overflow (issue #43)
d = [113997 1031298
     334453 37471]
PowerDivergenceTest(d)

# Pearson's Chi-Squared Test
# As in https://en.wikipedia.org/wiki/Pearson's_chi-squared_test#Fairness_of_dice
O = [5,8,9,8,10,20]
E = fill(10,6)
m = ChisqTest(O,E)

@test pvalue(m) ≈ 0.01990522033477436
@test m.stat ≈ 13.4
@test m.df == 5
@test m.n == 60
@test m.residuals ≈ reshape([-1.5811388300841895,-0.6324555320336759,-0.31622776601683794,-0.6324555320336759, 0.0, 3.162277660168379],size(O))

# As in https://en.wikipedia.org/wiki/Pearson's_chi-squared_test#Goodness_of_fit
O = [44,56]
E = fill(50,2)
m = ChisqTest(O,E)

@test pvalue(m) ≈ 0.23013934044341544
@test m.stat ≈ 1.44
@test m.df == 1
@test m.n == 100
@test m.residuals ≈ reshape([-0.848528137423857,0.848528137423857],size(O))

end
