using HypothesisTests, Base.Test

@test SWCoeffs(3).A == [sqrt(2.0)/2.0]
@test SWCoeffs(3).N == 3

S = SWCoeffs(10)
@test S.N == 10
@test length(S) == 10
@test endof(S) == 10
@test length(S.A) == 5

@test S.A[4] == S[4] == -S[7]
@test S.A[2] == S[2] == -S[9]

S = SWCoeffs(11)
@test S.N == 11
@test length(S) == 11
@test endof(S) == 11
@test length(S.A) == 5

@test S[6] == 0.0
@test S.A[5] == S[5] == -S[7]
@test S.A[3] == S[3] == -S[9]
@test S.A[1] == S[1] == -S[11]

#anti-symmetry
S = SWCoeffs(20)
@test all([S[i] == -S[end-i+1] for i in 1:length(S)])


# Values obtained by calling `_swilk` fortran subroutine directly.

SWc10 = SWCoeffs(10)
a = SWc10.A .- [0.573715, 0.32897, 0.214349, 0.122791, 0.0400887]
@test norm(a,1) ≈ 0.0 atol=S.N*eps(Float32)

SWc20 = SWCoeffs(20)
b = SWc20.A .- [0.473371, 0.32174, 0.255663, 0.208297, 0.16864, 0.133584, 0.101474, 0.0712893, 0.0423232, 0.0140351]
@test norm(b,1) ≈ 0.0 atol=S.N*eps(Float32)


srand(1)
X = sort(rand(Float64, 10))
@test HypothesisTests.swstat(X, SWc10) ≈ 0.8453033192812377 atol=1.2e-8
@test pvalue(HypothesisTests.swstat(X, SWc10), SWc10) ≈ 0.05106360837817192 atol=4.9e-8

# **Worked Example** (Section 4) from
# PATRICK ROYSTON Approximating the Shapiro-Wilk W-test for non-normality
# *Statistics and Computing* (1992) **2**, 117-119

X = [48.4, 49.0, 59.5, 59.6, 60.7, 88.8, 98.2, 109.4, 169.1, 227.1]
SWc = SWCoeffs(length(X))
@test all(abs.(SWc.A .- [0.5737, 0.3290, 0.2143, 0.1228, 0.0401]) .< 5.0e-5)
W = HypothesisTests.swstat(X, SWc)
@test W ≈ 0.8078 atol=2.9e-5
pval = pvalue(W, SWc)
@test pval ≈ 0.018 atol=4.7e-5

t = ShapiroWilkTest(X)
@test t.W == W
@test pvalue(t) == pval
@test t.N1 == length(X)
@test t.SWc.N == length(X)
@test length(t.SWc) == length(X)

@test pvalue(t) == pval
