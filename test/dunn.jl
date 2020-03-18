using HypothesisTests, Test
using HypothesisTests: default_tail, zscores, pvalues
using StatsBase: dof

@testset "Dunn" begin

# https://cran.r-project.org/web/packages/dunn.test/dunn.test.pdf
x = [2.9, 3.0, 2.5, 2.6, 3.2] # normal subjects
y = [3.8, 2.7, 4.0, 2.4]      # with obstructive airway disease
z = [2.8, 3.4, 3.7, 2.2, 2.0] # with asbestosis
t1 = HypothesisTests.DunnTest(x, y, z)
t2 = HypothesisTests.DunnTest(x, y, z; adj=:bonferroni)
t3 = HypothesisTests.DunnTest(x, y, z; adj=:sidak)

@test dof(t1) == 2
@test zscores(t1)[1=>2] ≈ 0.6414269805898184
@test zscores(t1)[1=>3] ≈ -0.22677868380553654
@test zscores(t1)[2=>3] ≈ -0.8552359741197582
@test t1.adjustment == :none
@test collect(values(pvalues(t1))) ≈ [ 0.4102979198777204, 0.2606226540574074, 0.19621026223829088]
@test zscores(t1) == zscores(t2)
@test t2.adjustment == :bonferroni
@test collect(values(pvalues(t2))) ≈ [ 1.0, 0.7818679621722222, 0.5886307867148727]
@test t3.adjustment == :sidak
@test collect(values(pvalues(t3))) ≈ [ 0.7949319606561769, 0.5957980356371838, 0.480689179999353]

show(IOBuffer(), t1)

# http://www.brightstat.com/index.php?option=com_content&task=view&id=41&Itemid=1&limit=1&limitstart=2
city1 = [68, 93, 123, 83, 108, 122]
city2 = [119, 116, 101, 103, 113, 84]
city3 = [70, 68, 54, 73, 81, 68]
city4 = [61, 54, 59, 67, 59, 70]
t = HypothesisTests.DunnTest(city1, city2, city3, city4)

# significant differences
pv = pvalues(t)
@testset for ij in [1=>3, 1=>4, 2=>3, 2=>4]
    @test pv[ij] < 0.05
end

end
