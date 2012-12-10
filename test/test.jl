load("test")
load("src/Wilcoxon")
using Wilcoxon, Test

## MANN-WHITNEY U

# Basic exact test
@test abs(p_value([1:10], [2.1:2:21], ExactMannWhitneyUTest) - 0.0232) <= 1e-4
@test abs(p_value([2.1:2:21], [1:10], ExactMannWhitneyUTest) - 0.0232) <= 1e-4
@test abs(p_value([1.5:10:100], [2.1:2:21], ExactMannWhitneyUTest) - 0.0068) <= 1e-4
@test abs(p_value([2.1:2:21], [1.5:10:100], ExactMannWhitneyUTest) - 0.0068) <= 1e-4

# Exact with ties
@test abs(p_value([1:10], [1:10], ExactMannWhitneyUTest) - 1) <= 1e-4
@test abs(p_value([1:10], [2:11], ExactMannWhitneyUTest) - 0.5096) <= 1e-4
@test abs(p_value([2:11], [1:10], ExactMannWhitneyUTest) - 0.5096) <= 1e-4
@test abs(p_value([1:10], [1:5, ones(5)], ExactMannWhitneyUTest) - 0.0057) <= 1e-4
@test abs(p_value([1:5, ones(5)], [1:10], ExactMannWhitneyUTest) - 0.0057) <= 1e-4

# Approximate test
@test abs(p_value([1:10], [2.1:2:21], ApproximateMannWhitneyUTest) - 0.0257) <= 1e-4
@test abs(p_value([2.1:2:21], [1:10], ApproximateMannWhitneyUTest) - 0.0257) <= 1e-4
@test abs(p_value([1.5:10:100], [2.1:2:21], ApproximateMannWhitneyUTest) - 0.0091) <= 1e-4
@test abs(p_value([2.1:2:21], [1.5:10:100], ApproximateMannWhitneyUTest) - 0.0091) <= 1e-4

# Approximate with ties
@test abs(p_value([1:10], [1:10], ApproximateMannWhitneyUTest) - 1) <= 1e-4
@test abs(p_value([1:10], [2:11], ApproximateMannWhitneyUTest) - 0.4948) <= 1e-4
@test abs(p_value([2:11], [1:10], ApproximateMannWhitneyUTest) - 0.4948) <= 1e-4
@test abs(p_value([1:10], [1:5, ones(5)], ApproximateMannWhitneyUTest) - 0.0076) <= 1e-4
@test abs(p_value([1:5, ones(5)], [1:10], ApproximateMannWhitneyUTest) - 0.0076) <= 1e-4

## WILCOXON SIGNED RANK

# Basic exact test
@test abs(p_value([1:10], [2:2:20], ExactSignedRankTest) - 0.0020) <= 1e-4
@test abs(p_value([2:2:20], [1:10], ExactSignedRankTest) - 0.0020) <= 1e-4
@test abs(p_value([1:10], [2:2:16, -1, 1], ExactSignedRankTest) - 0.4316) <= 1e-4
@test abs(p_value([2:2:16, -1, 1], [1:10], ExactSignedRankTest) - 0.4316) <= 1e-4

# Exact with ties
@test abs(p_value([1:10], [1:10], ExactSignedRankTest) - 1) <= 1e-4
@test abs(p_value([1:10], [2:11], ExactSignedRankTest) - 0.0020) <= 1e-4
@test abs(p_value([2:11], [1:10], ExactSignedRankTest) - 0.0020) <= 1e-4
@test abs(p_value([1:10], [1:5, ones(5)], ExactSignedRankTest) - 0.0625) <= 1e-4
@test abs(p_value([1:5, ones(5)], [1:10], ExactSignedRankTest) - 0.0625) <= 1e-4

# Approximate test
@test abs(p_value([1:10], [2:2:20], ApproximateSignedRankTest) - 0.005922) <= 1e-6
@test abs(p_value([2:2:20], [1:10], ApproximateSignedRankTest) - 0.005922) <= 1e-6
@test abs(p_value([1:10], [2:2:16, -1, 1], ApproximateSignedRankTest) - 0.4148) <= 1e-4
@test abs(p_value([2:2:16, -1, 1], [1:10], ApproximateSignedRankTest) - 0.4148) <= 1e-4

# Approximate with ties
@test abs(p_value([1:10], [2:11], ApproximateSignedRankTest) - 0.001904) <= 1e-6
@test abs(p_value([1:10], [1:10], ApproximateSignedRankTest) - 1) <= 1e-4
@test abs(p_value([1:10], [1:5, ones(5)], ApproximateSignedRankTest) - 0.05906) <= 1e-5