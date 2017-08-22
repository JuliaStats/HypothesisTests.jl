using HypothesisTests, Base.Test
using HypothesisTests: tail

# Basic exact test
@test abs(pvalue(ExactSignedRankTest([1:10;], [2:2:20;])) - 0.0020) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([2:2:20;], [1:10;])) - 0.0020) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([1:10;], [2:2:16; -1; 1])) - 0.4316) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([2:2:16; -1; 1], [1:10;])) - 0.4316) <= 1e-4
@test tail(ExactSignedRankTest([1:10;], [2:2:20;])) == :both
show(IOBuffer(), ExactSignedRankTest([1:10;], [2:2:20;]))

# Exact with ties
@test abs(pvalue(ExactSignedRankTest([1:10;], [1:10;])) - 1) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([1:10;], [2:11;])) - 0.0020) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([2:11;], [1:10;])) - 0.0020) <= 1e-4
@test abs(pvalue(ExactSignedRankTest(1:10, [1:5; ones(5)])) - 0.0625) <= 1e-4
@test abs(pvalue(ExactSignedRankTest([1:5; ones(5)], [1:10;])) - 0.0625) <= 1e-4
show(IOBuffer(), ExactSignedRankTest([1:10;], [1:10;]))

# Approximate test
@test abs(pvalue(ApproximateSignedRankTest([1:10;], [2:2:20;])) - 0.005922) <= 1e-6
@test abs(pvalue(ApproximateSignedRankTest([2:2:20;], [1:10;])) - 0.005922) <= 1e-6
@test abs(pvalue(ApproximateSignedRankTest([1:10;], [2:2:16; -1; 1])) - 0.4148) <= 1e-4
@test abs(pvalue(ApproximateSignedRankTest([2:2:16; -1; 1], [1:10;])) - 0.4148) <= 1e-4
@test tail(ApproximateSignedRankTest([1:10;], [2:2:20;])) == :both
show(IOBuffer(), ApproximateSignedRankTest([1:10;], [2:2:20;]))

# Approximate with ties
@test abs(pvalue(ApproximateSignedRankTest([1:10;], [1:10;])) - 1) <= 1e-4
@test abs(pvalue(ApproximateSignedRankTest([1:10;], [2:11;])) - 0.001904) <= 1e-6
@test abs(pvalue(ApproximateSignedRankTest([2:11;], [1:10;])) - 0.001904) <= 1e-6
@test abs(pvalue(ApproximateSignedRankTest([1:10;], [1:5; ones(5)])) - 0.05906) <= 1e-5
@test abs(pvalue(ApproximateSignedRankTest([1:5; ones(5)], 1:10)) - 0.05906) <= 1e-5
show(IOBuffer(), ApproximateSignedRankTest([1:10;], [1:10;]))

# # Tests for automatic selection
@test abs(pvalue(SignedRankTest([1:10;], [2:2:20;])) - 0.0020) <= 1e-4
@test abs(pvalue(SignedRankTest([1:10;], [2:11;])) - 0.0020) <= 1e-4
@test tail(SignedRankTest([1:10;], [2:2:20;])) == :both
show(IOBuffer(), SignedRankTest([1:10;], [2:2:20;]))

# One Sample tests
# P-value computed using R wilcox.test
@test abs(pvalue(SignedRankTest([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] - 10.1)) - 0.09460449) <= 1e-4
# P-value computed using R wilcox.test
@test abs(pvalue(SignedRankTest([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] - 10.1)) - 0.1928101) <= 1e-4

# One Sample tests with ties
# P-value computed using R package exactRankTests wilcox.exact
@test abs(pvalue(SignedRankTest([1,2,3,4,5,6,7,10,10,10,10,10,13,14,15] - 10.1)) - 0.04052734) <= 1e-4
# P-value computed using R wilcox.test
@test abs(pvalue(SignedRankTest([1,2,3,4,5,6,7,10,10,10,10,10,13,14,15,16] - 10.1)) - 0.1021964) <= 1e-4

# Test confidence interval
x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
@test isapprox(confint(ExactSignedRankTest(x))[1], 3.3, atol=1e-4)
@test isapprox(confint(ExactSignedRankTest(x))[2], 15.5, atol=1e-4)
@test isapprox(confint(ApproximateSignedRankTest(x))[1], 3.3, atol=1e-4)
@test isapprox(confint(ApproximateSignedRankTest(x))[2], 15.5, atol=1e-4)
@test isapprox(confint(SignedRankTest(x); tail=:left)[1], 4.45, atol=1e-4)
@test isapprox(confint(SignedRankTest(x); tail=:right)[2], 14.45, atol=1e-4)
