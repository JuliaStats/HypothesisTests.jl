using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Friedman" begin
a1 = [4, 6, 3, 4, 3, 2, 2, 7, 6, 5, 3, 8, 6, 1, 5, 7]
a2 = [5, 6, 8, 7, 7, 8, 4, 6, 4, 5, 2, 4, 4, 3, 8, 9]
a3 = [2, 4, 4, 3, 2, 2, 1, 4, 3, 2, 5, 3, 5, 7, 6, 1]
t = FriedmanTest(a1, a2, a3)
@test t.n == length(a1) == length(a2) == length(a3)
@test t.df == 2
@test t.rank_sums == [33.5, 39, 23.5]
@test t.Q ≈ 8.098360655737705
@test pvalue(t) ≈ 0.017436661128762587
@test default_tail(t) == :right
show(IOBuffer(), t)

# https://github.com/scipy/scipy/blob/main/scipy/stats/_stats_py.py
b1 = [72, 96, 88, 92, 74, 76, 82]
b2 = [120, 120, 132, 120, 101, 96, 112]
b3 = [76, 95, 104, 96, 84, 72, 76]
t = FriedmanTest(b1, b2, b3)
@test t.n == length(b1) == length(b2) == length(b3)
@test t.df == 2
@test t.rank_sums == [10, 21, 11]
@test t.Q ≈ 10.57142857142857
@test pvalue(t) ≈ 0.005063414171757498
show(IOBuffer(), t)

c1 = [-0.263, 0.216, 0.667, 0.487, 0.376, 0.333, 0.353, 0.132, -0.864, 0.567]
c2 = [0.384, -0.432, 0.987, 0.084, 0.974, 0.415, 0.121, -0.369, 0.847, 0.076]
c3 = [0.129, 0.312, -0.453, 0.927, 0.594, 0.767, -0.009, 0.735, 0.853, 0.697]
c4 = [0.292, 0.203, 0.983, -0.690, 0.906, -0.035, 0.209, 0.309, 0.516, 0.850]
c5 = [0.256, 0.941, 0.260, 0.628, -0.189, 0.687, 0.443, 0.767, 0.636, -0.931]
t = FriedmanTest(c1, c2, c3, c4, c5)
@test t.n == length(c1) == length(c2) == length(c3) == length(c4) == length(c5)
@test t.df == 4
@test t.rank_sums == [24, 30, 34, 29, 33]
@test t.Q ≈ 2.4799999999999999898
@test pvalue(t) ≈ 0.6482206481834754
show(IOBuffer(), t)
end