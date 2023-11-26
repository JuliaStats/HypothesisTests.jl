using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Friedman" begin
a1 = [4, 6, 3, 4, 3, 2, 2, 7, 6, 5]
a2 = [5, 6, 8, 7, 7, 8, 4, 6, 4, 5]
a3 = [2, 4, 4, 3, 2, 2, 1, 4, 3, 2]
t = FriedmanTest(a1, a2, a3)
@test t.n == length(a1) == length(a2) == length(a3)
@test t.df == 2
@test t.rank_sums ==  [21.5, 27, 11.5]
@test t.chi_sq ≈ 13.351351351351344
@test pvalue(t) ≈ 0.0012612201221243592
show(IOBuffer(), t)

b1 = [1.2, 1.9, 2.1]
b2 = [3.4, 1.9, 5.6]
b3 = [1.3, 4.4, 9.9]
t = FriedmanTest(b1, b2, b3)
@test t.n == length(b1) == length(b2) == length(b3)
@test t.df == 2
@test t.rank_sums ==  [3.5, 6.5, 8.0]
@test t.chi_sq ≈ 3.8181818181818183
@test pvalue(t) ≈ 0.14821506633752016
show(IOBuffer(), t)

c1 = [0.37551722, -0.91571614, -0.72111682, 0.9376363, 0.38910423]
c2 = [0.54212944, 0.66137236, -0.67373387, 0.4027261, -0.32210339]
c3 = [-0.98970611, 0.05812077, 0.34603985, -0.2855241, 0.65564481]
c4 = [-0.65447731, -0.02035673, 0.89509799, 0.86444013, 0.74524022]
t = FriedmanTest(c1, c2, c3, c4)
@test t.n == length(c1) == length(c2) == length(c3) == length(c4)
@test t.df == 3
@test t.rank_sums ==  [11, 13, 11, 15]
@test t.chi_sq ≈ 1.3199999999999932
@test pvalue(t) ≈ 0.7243894461724858
show(IOBuffer(), t)
end
