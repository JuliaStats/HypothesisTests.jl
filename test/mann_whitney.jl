using HypothesisTests, Base.Test

# Basic exact test
@test abs(pvalue(ExactMannWhitneyUTest([1:10], [2.1:2:21])) - 0.0232) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([2.1:2:21], [1:10])) - 0.0232) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([1.5:10:100], [2.1:2:21])) - 0.0068) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([2.1:2:21], [1.5:10:100])) - 0.0068) <= 1e-4
show(IOBuffer(), ExactMannWhitneyUTest([1:10], [2.1:2:21]))

# Exact with ties
@test abs(pvalue(ExactMannWhitneyUTest([1:10], [1:10])) - 1) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([1:10], [2:11])) - 0.5096) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([2:11], [1:10])) - 0.5096) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([1:10], [1:5, ones(5)])) - 0.0057) <= 1e-4
@test abs(pvalue(ExactMannWhitneyUTest([1:5, ones(5)], [1:10])) - 0.0057) <= 1e-4
show(IOBuffer(), ExactMannWhitneyUTest([1:10], [1:10]))

# Approximate test
@test abs(pvalue(ApproximateMannWhitneyUTest([1:10], [2.1:2:21])) - 0.0257) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([2.1:2:21], [1:10])) - 0.0257) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([1.5:10:100], [2.1:2:21])) - 0.0091) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([2.1:2:21], [1.5:10:100])) - 0.0091) <= 1e-4
show(IOBuffer(), ApproximateMannWhitneyUTest([1:10], [2.1:2:21]))

# Approximate with ties
@test abs(pvalue(ApproximateMannWhitneyUTest([1:10], [1:10])) - 1) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([1:10], [2:11])) - 0.4948) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([2:11], [1:10])) - 0.4948) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([1:10], [1:5, ones(5)])) - 0.0076) <= 1e-4
@test abs(pvalue(ApproximateMannWhitneyUTest([1:5, ones(5)], [1:10])) - 0.0076) <= 1e-4
show(IOBuffer(), ApproximateMannWhitneyUTest([1:10], [1:10]))

# Tests for automatic selection
@test abs(pvalue(MannWhitneyUTest([1:10], [2.1:2:21])) - 0.0232) <= 1e-4
@test abs(pvalue(MannWhitneyUTest([1:10], [2:11])) - 0.4948) <= 1e-4
show(IOBuffer(), MannWhitneyUTest([1:10], [2.1:2:21]))