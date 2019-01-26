var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "HypothesisTests package",
    "title": "HypothesisTests package",
    "category": "page",
    "text": ""
},

{
    "location": "#HypothesisTests-package-1",
    "page": "HypothesisTests package",
    "title": "HypothesisTests package",
    "category": "section",
    "text": "This package implements several hypothesis tests in Julia.Pages = [\"methods.md\", \"parametric.md\", \"nonparametric.md\", \"time_series.md\"]\nDepth = 2"
},

{
    "location": "methods/#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "methods/#Methods-1",
    "page": "Methods",
    "title": "Methods",
    "category": "section",
    "text": ""
},

{
    "location": "methods/#StatsBase.confint",
    "page": "Methods",
    "title": "StatsBase.confint",
    "category": "function",
    "text": "confint(test::HypothesisTest, alpha = 0.05; tail = :both)\n\nCompute a confidence interval C with coverage 1-alpha.\n\nIf tail is :both (default), then a two-sided confidence interval is returned. If tail is :left or :right, then a one-sided confidence interval is returned.\n\nnote: Note\nMost of the implemented confidence intervals are strongly consistent, that is, the confidence interval with coverage 1-alpha does not contain the test statistic under h_0 if and only if the corresponding test rejects the null hypothesis h_0 θ = θ_0:    C (x 1  α) = θ  p_θ (x)  αwhere p_θ is the pvalue of the corresponding test.\n\n\n\n\n\n"
},

{
    "location": "methods/#StatsBase.confint-Tuple{BinomialTest}",
    "page": "Methods",
    "title": "StatsBase.confint",
    "category": "method",
    "text": "confint(test::BinomialTest, alpha = 0.05; tail = :both, method = :clopper_pearson)\n\nCompute a confidence interval with coverage 1-alpha for a binomial proportion using one of the following methods. Possible values for method are:\n\n:clopper_pearson (default): Clopper-Pearson interval is based on the binomial distribution. The empirical coverage is never less than the nominal coverage of 1-alpha; it is usually too conservative.\n:wald: Wald (or normal approximation) interval relies on the standard approximation of the actual binomial distribution by a normal distribution. Coverage can be erratically poor for success probabilities close to zero or one.\n:wilson: Wilson score interval relies on a normal approximation. In contrast to :wald, the standard deviation is not approximated by an empirical estimate, resulting in good empirical coverages even for small numbers of draws and extreme success probabilities.\n:jeffrey: Jeffreys interval is a Bayesian credible interval obtained by using a non-informative Jeffreys prior. The interval is very similar to the Wilson interval.\n:agresti_coull: Agresti-Coull interval is a simplified version of the Wilson interval; both are centered around the same value. The Agresti Coull interval has higher or equal coverage.\n:arcsine: Confidence interval computed using the arcsine transformation to make var(p) independent of the probability p.\n\nReferences\n\nBrown, L.D., Cai, T.T., and DasGupta, A. Interval estimation for a binomial proportion. Statistical Science, 16(2):101–117, 2001.\n\nExternal links\n\nBinomial confidence interval on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "methods/#StatsBase.confint-Tuple{PowerDivergenceTest}",
    "page": "Methods",
    "title": "StatsBase.confint",
    "category": "method",
    "text": "confint(test::PowerDivergenceTest, alpha = 0.05; tail = :both, method = :sison_glaz)\n\nCompute a confidence interval with coverage 1-alpha for multinomial proportions using one of the following methods. Possible values for method are:\n\n:sison_glaz (default): Sison-Glaz intervals\n:bootstrap: Bootstrap intervals\n:quesenberry_hurst: Quesenberry-Hurst intervals\n:gold: Gold intervals (asymptotic simultaneous intervals)\n\nReferences\n\nAgresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013.\nSison, C.P and Glaz, J. Simultaneous confidence intervals and sample size determination for multinomial proportions. Journal of the American Statistical Association, 90:366-369, 1995.\nQuesensberry, C.P. and Hurst, D.C. Large Sample Simultaneous Confidence Intervals for Multinational Proportions. Technometrics, 6:191-195, 1964.\nGold, R. Z. Tests Auxiliary to χ^2 Tests in a Markov Chain. Annals of Mathematical Statistics, 30:56-74, 1963.\n\n\n\n\n\n"
},

{
    "location": "methods/#StatsBase.confint-Tuple{FisherExactTest}",
    "page": "Methods",
    "title": "StatsBase.confint",
    "category": "method",
    "text": "confint(x::FisherExactTest, alpha::Float64=0.05; tail=:both, method=:central)\n\nCompute a confidence interval with coverage 1 - alpha. One-sided intervals are based on Fisher\'s non-central hypergeometric distribution. For tail = :both, the only method implemented yet is the central interval (:central).\n\nnote: Note\nSince the p-value is not necessarily unimodal, the corresponding confidence region might not be an interval.\n\nReferences\n\nGibbons, J.D, Pratt, J.W. P-values: Interpretation and Methodology, American Statistican, 29(1):20-25, 1975.\nFay, M.P., Supplementary material to \"Confidence intervals that match Fisher’s exact or Blaker’s exact tests\". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374, link\n\n\n\n\n\n"
},

{
    "location": "methods/#Confidence-interval-1",
    "page": "Methods",
    "title": "Confidence interval",
    "category": "section",
    "text": "confint\nconfint(::BinomialTest)\nconfint(::PowerDivergenceTest)\nconfint(::FisherExactTest)"
},

{
    "location": "methods/#HypothesisTests.pvalue",
    "page": "Methods",
    "title": "HypothesisTests.pvalue",
    "category": "function",
    "text": "pvalue(test::HypothesisTest; tail = :both)\n\nCompute the p-value for a given significance test.\n\nIf tail is :both (default), then the p-value for the two-sided test is returned. If tail is :left or :right, then a one-sided test is performed.\n\n\n\n\n\n"
},

{
    "location": "methods/#HypothesisTests.pvalue-Tuple{FisherExactTest}",
    "page": "Methods",
    "title": "HypothesisTests.pvalue",
    "category": "method",
    "text": "pvalue(x::FisherExactTest; tail = :both, method = :central)\n\nCompute the p-value for a given Fisher exact test.\n\nThe one-sided p-values are based on Fisher\'s non-central hypergeometric distribution f_ω(i) with odds ratio ω:\n\n    beginalign*\n        p_ω^(textleft) =sum_i  a f_ω(i)\n        p_ω^(textright) =sum_i  a f_ω(i)\n    endalign*\n\nFor tail = :both, possible values for method are:\n\n:central (default): Central interval, i.e. the p-value is two times the minimum of the one-sided p-values.\n:minlike: Minimum likelihood interval, i.e. the p-value is computed by summing all tables with the same marginals that are equally or less probable:\n    p_ω = sum_f_ω(i) f_ω(a) f_ω(i)\n\nReferences\n\nGibbons, J.D., Pratt, J.W., P-values: Interpretation and Methodology, American Statistican, 29(1):20-25, 1975.\nFay, M.P., Supplementary material to \"Confidence intervals that match Fisher’s exact or Blaker’s exact tests\". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374, link\n\n\n\n\n\n"
},

{
    "location": "methods/#p-value-1",
    "page": "Methods",
    "title": "p-value",
    "category": "section",
    "text": "pvalue\npvalue(::FisherExactTest)"
},

{
    "location": "parametric/#",
    "page": "Parametric tests",
    "title": "Parametric tests",
    "category": "page",
    "text": ""
},

{
    "location": "parametric/#Parametric-tests-1",
    "page": "Parametric tests",
    "title": "Parametric tests",
    "category": "section",
    "text": ""
},

{
    "location": "parametric/#HypothesisTests.PowerDivergenceTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.PowerDivergenceTest",
    "category": "type",
    "text": "PowerDivergenceTest(x[, y]; lambda = 1.0, theta0 = ones(length(x))/length(x))\n\nPerform a Power Divergence test.\n\nIf y is not given and x is a matrix with one row or column, or x is a vector, then a goodness-of-fit test is performed (x is treated as a one-dimensional contingency table). In this case, the hypothesis tested is whether the population probabilities equal those in theta0, or are all equal if theta0 is not given.\n\nIf x is a matrix with at least two rows and columns, it is taken as a two-dimensional contingency table. Otherwise, x and y must be vectors of the same length. The contingency table is calculated using the counts function from the StatsBase package. Then the power divergence test is conducted under the null hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency table is the product of the row and column marginals.\n\nNote that the entries of x (and y if provided) must be non-negative integers.\n\nThe power divergence test is given by\n\n    dfrac2λ(λ+1)sum_i=1^I sum_j=1^J n_ij left(n_ij\n    hatn_ij)^λ -1right\n\nwhere n_ij is the cell count in the i th row and j th column and λ is a real number determining the nature of the test to be performed:\n\nλ = 1: equal to Pearson\'s chi-squared statistic\nλ to 0: converges to the likelihood ratio test statistic\nλ to -1: converges to the minimum discrimination information statistic (Gokhale and Kullback, 1978)\nλ = -2: equals Neyman modified chi-squared (Neyman, 1949)\nλ = -12: equals the Freeman-Tukey statistic (Freeman and Tukey, 1950).\n\nUnder regularity conditions, the asymptotic distributions are identical (see Drost et. al. 1989). The χ^2 null approximation works best for λ near 23.\n\nImplements: pvalue, confint(::PowerDivergenceTest)\n\nReferences\n\nAgresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013.\n\n\n\n\n\n"
},

{
    "location": "parametric/#Power-divergence-test-1",
    "page": "Parametric tests",
    "title": "Power divergence test",
    "category": "section",
    "text": "PowerDivergenceTest"
},

{
    "location": "parametric/#HypothesisTests.ChisqTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.ChisqTest",
    "category": "function",
    "text": "ChisqTest(x[, y][, theta0 = ones(length(x))/length(x)])\n\nPerform a Pearson chi-squared test (equivalent to a PowerDivergenceTest with λ = 1).\n\nIf y is not given and x is a matrix with one row or column, or x is a vector, then a goodness-of-fit test is performed (x is treated as a one-dimensional contingency table). In this case, the hypothesis tested is whether the population probabilities equal those in theta0, or are all equal if theta0 is not given.\n\nIf x is a matrix with at least two rows and columns, it is taken as a two-dimensional contingency table. Otherwise, x and y must be vectors of the same length. The contingency table is calculated using counts function from the StatsBase package. Then the power divergence test is conducted under the null hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency table is the product of the row and column marginals.\n\nNote that the entries of x (and y if provided) must be non-negative integers.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#Pearson-chi-squared-test-1",
    "page": "Parametric tests",
    "title": "Pearson chi-squared test",
    "category": "section",
    "text": "ChisqTest"
},

{
    "location": "parametric/#HypothesisTests.MultinomialLRT",
    "page": "Parametric tests",
    "title": "HypothesisTests.MultinomialLRT",
    "category": "function",
    "text": "MultinomialLRT(x[, y][, theta0 = ones(length(x))/length(x)])\n\nPerform a multinomial likelihood ratio test (equivalent to a PowerDivergenceTest with λ = 0).\n\nIf y is not given and x is a matrix with one row or column, or x is a vector, then a goodness-of-fit test is performed (x is treated as a one-dimensional contingency table). In this case, the hypothesis tested is whether the population probabilities equal those in theta0, or are all equal if theta0 is not given.\n\nIf x is a matrix with at least two rows and columns, it is taken as a two-dimensional contingency table. Otherwise, x and y must be vectors of the same length. The contingency table is calculated using counts function from the StatsBase package. Then the power divergence test is conducted under the null hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency table is the product of the row and column marginals.\n\nNote that the entries of x (and y if provided) must be non-negative integers.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#Multinomial-likelihood-ratio-test-1",
    "page": "Parametric tests",
    "title": "Multinomial likelihood ratio test",
    "category": "section",
    "text": "MultinomialLRT"
},

{
    "location": "parametric/#HypothesisTests.OneSampleTTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.OneSampleTTest",
    "category": "type",
    "text": "OneSampleTTest(xbar::Real, stddev::Real, n::Int, μ0::Real = 0)\n\nPerform a one sample t-test of the null hypothesis that n values with mean xbar and sample standard deviation stddev  come from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\nOneSampleTTest(v::AbstractVector{T<:Real}, μ0::Real = 0)\n\nPerform a one sample t-test of the null hypothesis that the data in vector v comes from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\nOneSampleTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, μ0::Real = 0)\n\nPerform a paired sample t-test of the null hypothesis that the differences between pairs of values in vectors x and y come from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#HypothesisTests.EqualVarianceTTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.EqualVarianceTTest",
    "category": "type",
    "text": "EqualVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})\n\nPerform a two-sample t-test of the null hypothesis that x and y come from distributions with equal means and variances against the alternative hypothesis that the distributions have different means but equal variances.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#HypothesisTests.UnequalVarianceTTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.UnequalVarianceTTest",
    "category": "type",
    "text": "UnequalVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})\n\nPerform an unequal variance two-sample t-test of the null hypothesis that x and y come from distributions with equal means against the alternative hypothesis that the distributions have different means.\n\nThis test is sometimes known as Welch\'s t-test. It differs from the equal variance t-test in that it computes the number of degrees of freedom of the test using the Welch-Satterthwaite equation:\n\n    ν_χ  fracleft(sum_i=1^n k_i s_i^2right)^2sum_i=1^n\n        frac(k_i s_i^2)^2ν_i\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#t-test-1",
    "page": "Parametric tests",
    "title": "t-test",
    "category": "section",
    "text": "OneSampleTTest\nEqualVarianceTTest\nUnequalVarianceTTest"
},

{
    "location": "parametric/#HypothesisTests.OneSampleZTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.OneSampleZTest",
    "category": "type",
    "text": "OneSampleZTest(xbar::Real, stddev::Real, n::Int, μ0::Real = 0)\n\nPerform a one sample z-test of the null hypothesis that n values with mean xbar and population standard deviation stddev  come from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\nOneSampleZTest(v::AbstractVector{T<:Real}, μ0::Real = 0)\n\nPerform a one sample z-test of the null hypothesis that the data in vector v comes from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\nOneSampleZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, μ0::Real = 0)\n\nPerform a paired sample z-test of the null hypothesis that the differences between pairs of values in vectors x and y come from a distribution with mean μ0 against the alternative hypothesis that the distribution does not have mean μ0.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#HypothesisTests.EqualVarianceZTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.EqualVarianceZTest",
    "category": "type",
    "text": "EqualVarianceZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})\n\nPerform a two-sample z-test of the null hypothesis that x and y come from distributions with equal means and variances against the alternative hypothesis that the distributions have different means but equal variances.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#HypothesisTests.UnequalVarianceZTest",
    "page": "Parametric tests",
    "title": "HypothesisTests.UnequalVarianceZTest",
    "category": "type",
    "text": "UnequalVarianceZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})\n\nPerform an unequal variance two-sample z-test of the null hypothesis that x and y come from distributions with equal means against the alternative hypothesis that the distributions have different means.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "parametric/#z-test-1",
    "page": "Parametric tests",
    "title": "z-test",
    "category": "section",
    "text": "OneSampleZTest\nEqualVarianceZTest\nUnequalVarianceZTest"
},

{
    "location": "nonparametric/#",
    "page": "Nonparametric tests",
    "title": "Nonparametric tests",
    "category": "page",
    "text": ""
},

{
    "location": "nonparametric/#Nonparametric-tests-1",
    "page": "Nonparametric tests",
    "title": "Nonparametric tests",
    "category": "section",
    "text": ""
},

{
    "location": "nonparametric/#HypothesisTests.OneSampleADTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.OneSampleADTest",
    "category": "type",
    "text": "OneSampleADTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)\n\nPerform a one-sample Anderson–Darling test of the null hypothesis that the data in vector x come from the distribution d against the alternative hypothesis that the sample is not drawn from d.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.KSampleADTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.KSampleADTest",
    "category": "type",
    "text": "KSampleADTest(xs::AbstractVector{<:Real}...; modified = true, nsim = 0)\n\nPerform a k-sample Anderson–Darling test of the null hypothesis that the data in the k vectors xs come from the same distribution against the alternative hypothesis that the samples come from different distributions.\n\nmodified parameter enables a modified test calculation for samples whose observations do not all coincide.\n\nIf nsim is equal to 0 (the default) the asymptotic calculation of p-value is used. If it is greater than 0, an estimation of p-values is used by generating nsim random splits of the pooled data on k samples, evaluating the AD statistics for each split, and computing the proportion of simulated values which are greater or equal to observed. This proportion is reported as p-value estimate.\n\nImplements: pvalue\n\nReferences\n\nF. W. Scholz and M. A. Stephens, K-Sample Anderson-Darling Tests, Journal of the American Statistical Association, Vol. 82, No. 399. (Sep., 1987), pp. 918-924.\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Anderson-Darling-test-1",
    "page": "Nonparametric tests",
    "title": "Anderson-Darling test",
    "category": "section",
    "text": "Available are both one-sample and k-sample tests.OneSampleADTest\nKSampleADTest"
},

{
    "location": "nonparametric/#HypothesisTests.BinomialTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.BinomialTest",
    "category": "type",
    "text": "BinomialTest(x::Integer, n::Integer, p::Real = 0.5)\nBinomialTest(x::AbstractVector{Bool}, p::Real = 0.5)\n\nPerform a binomial test of the null hypothesis that the distribution from which x successes were encountered in n draws (or alternatively from which the vector x was drawn) has success probability p against the alternative hypothesis that the success probability is not equal to p.\n\nComputed confidence intervals (confint) by default are Clopper-Pearson intervals.\n\nImplements: pvalue, confint(::BinomialTest)\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Binomial-test-1",
    "page": "Nonparametric tests",
    "title": "Binomial test",
    "category": "section",
    "text": "BinomialTest"
},

{
    "location": "nonparametric/#HypothesisTests.FisherExactTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.FisherExactTest",
    "category": "type",
    "text": "FisherExactTest(a::Integer, b::Integer, c::Integer, d::Integer)\n\nPerform Fisher\'s exact test of the null hypothesis that the success probabilities ac and bd are equal, that is the odds ratio (ac)  (bd) is one, against the alternative hypothesis that they are not equal.\n\nSee pvalue(::FisherExactTest) and confint(::FisherExactTest) for details about the computation of the default p-value and confidence interval, respectively.\n\nThe contingency table is structured as:\n\n- X1 X2\nY1 a b\nY2 c d\n\nnote: Note\nThe show function output contains the conditional maximum likelihood estimate of the odds ratio rather than the sample odds ratio; it maximizes the likelihood given by Fisher\'s non-central hypergeometric distribution.\n\nImplements: pvalue(::FisherExactTest), confint(::FisherExactTest)\n\nReferences\n\nFay, M.P., Supplementary material to \"Confidence intervals that match Fisher’s exact or Blaker’s exact tests\". Biostatistics, Volume 11, Issue 2, 1 April 2010, Pages 373–374, link\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Fisher-exact-test-1",
    "page": "Nonparametric tests",
    "title": "Fisher exact test",
    "category": "section",
    "text": "FisherExactTest"
},

{
    "location": "nonparametric/#HypothesisTests.ExactOneSampleKSTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ExactOneSampleKSTest",
    "category": "type",
    "text": "ExactOneSampleKSTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)\n\nPerform a one-sample exact Kolmogorov–Smirnov test of the null hypothesis that the data in vector x comes from the distribution d against the alternative hypothesis that the sample is not drawn from d.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ApproximateOneSampleKSTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ApproximateOneSampleKSTest",
    "category": "type",
    "text": "ApproximateOneSampleKSTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)\n\nPerform an asymptotic one-sample Kolmogorov–Smirnov test of the null hypothesis that the data in vector x comes from the distribution d against the alternative hypothesis that the sample is not drawn from d.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ApproximateTwoSampleKSTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ApproximateTwoSampleKSTest",
    "category": "type",
    "text": "ApproximateTwoSampleKSTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})\n\nPerform an asymptotic two-sample Kolmogorov–Smirnov-test of the null hypothesis that x and y are drawn from the same distribution against the alternative hypothesis that they come from different distributions.\n\nImplements: pvalue\n\nExternal links\n\nApproximation of one-sided test (Encyclopedia of Mathematics) \n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Kolmogorov–Smirnov-test-1",
    "page": "Nonparametric tests",
    "title": "Kolmogorov–Smirnov test",
    "category": "section",
    "text": "Available are an exact one-sample test and approximate (i.e. asymptotic) one- and two-sample tests.ExactOneSampleKSTest\nApproximateOneSampleKSTest\nApproximateTwoSampleKSTest"
},

{
    "location": "nonparametric/#HypothesisTests.KruskalWallisTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.KruskalWallisTest",
    "category": "type",
    "text": "KruskalWallisTest(groups::AbstractVector{<:Real}...)\n\nPerform Kruskal-Wallis rank sum test of the null hypothesis that the groups mathcalG come from the same distribution against the alternative hypothesis that that at least one group stochastically dominates one other group.\n\nThe Kruskal-Wallis test is an extension of the Mann-Whitney U test to more than two groups.\n\nThe p-value is computed using a χ^2 approximation to the distribution of the test statistic H_c=fracHC:\n\n    beginalign*\n    H  = frac12n(n+1) sum_g  mathcalG fracR_g^2n_g - 3(n+1)\n    C  = 1-frac1n^3-nsum_t  mathcalT (t^3-t)\n    endalign*\n\nwhere mathcalT is the set of the counts of tied values at each tied position, n is the total number of observations across all groups, and n_g and R_g are the number of observations and the rank sum in group g, respectively. See references for further details.\n\nImplements: pvalue\n\nReferences\n\nMeyer, J.P, Seaman, M.A., Expanded tables of critical values for the Kruskal-Wallis H statistic. Paper presented at the annual meeting of the American Educational Research Association, San Francisco, April 2006.\n\nExternal links\n\nKruskal-Wallis test on Wikipedia \n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Kruskal-Wallis-rank-sum-test-1",
    "page": "Nonparametric tests",
    "title": "Kruskal-Wallis rank sum test",
    "category": "section",
    "text": "KruskalWallisTest"
},

{
    "location": "nonparametric/#HypothesisTests.MannWhitneyUTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.MannWhitneyUTest",
    "category": "function",
    "text": "MannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})\n\nPerform a Mann-Whitney U test of the null hypothesis that the probability that an observation drawn from the same population as x is greater than an observation drawn from the same population as y is equal to the probability that an observation drawn from the same population as y is greater than an observation drawn from the same population as x against the alternative hypothesis that these probabilities are not equal.\n\nThe Mann-Whitney U test is sometimes known as the Wilcoxon rank-sum test.\n\nWhen there are no tied ranks and ≤50 samples, or tied ranks and ≤10 samples, MannWhitneyUTest performs an exact Mann-Whitney U test. In all other cases, MannWhitneyUTest performs an approximate Mann-Whitney U test. Behavior may be further controlled by using ExactMannWhitneyUTest or ApproximateMannWhitneyUTest directly.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ExactMannWhitneyUTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ExactMannWhitneyUTest",
    "category": "type",
    "text": "ExactMannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})\n\nPerform an exact Mann-Whitney U test of the null hypothesis that the probability that an observation drawn from the same population as x is greater than an observation drawn from the same population as y is equal to the probability that an observation drawn from the same population as y is greater than an observation drawn from the same population as x against the alternative hypothesis that these probabilities are not equal.\n\nWhen there are no tied ranks, the exact p-value is computed using the pwilcox function from the Rmath package. In the presence of tied ranks, a p-value is computed by exhaustive enumeration of permutations, which can be very slow for even moderately sized data sets.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ApproximateMannWhitneyUTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ApproximateMannWhitneyUTest",
    "category": "type",
    "text": "ApproximateMannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})\n\nPerform an approximate Mann-Whitney U test of the null hypothesis that the probability that an observation drawn from the same population as x is greater than an observation drawn from the same population as y is equal to the probability that an observation drawn from the same population as y is greater than an observation drawn from the same population as x against the alternative hypothesis that these probabilities are not equal.\n\nThe p-value is computed using a normal approximation to the distribution of the Mann-Whitney U statistic:\n\n    beginalign*\n        μ  = fracn_x n_y2\n        σ  = fracn_x n_y12left(n_x + n_y + 1 - fraca(n_x + n_y)(n_x +\n            n_y - 1)right)\n        a  = sum_t in mathcalT t^3 - t\n    endalign*\n\nwhere mathcalT is the set of the counts of tied values at each tied position.\n\nImplements: pvalue\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Mann-Whitney-U-test-1",
    "page": "Nonparametric tests",
    "title": "Mann-Whitney U test",
    "category": "section",
    "text": "MannWhitneyUTest\nExactMannWhitneyUTest\nApproximateMannWhitneyUTest"
},

{
    "location": "nonparametric/#HypothesisTests.SignTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.SignTest",
    "category": "type",
    "text": "SignTest(x::AbstractVector{T<:Real}, median::Real = 0)\nSignTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, median::Real = 0)\n\nPerform a sign test of the null hypothesis that the distribution from which x (or x - y if y is provided) was drawn has median median against the alternative hypothesis that the median is not equal to median.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Sign-test-1",
    "page": "Nonparametric tests",
    "title": "Sign test",
    "category": "section",
    "text": "SignTest"
},

{
    "location": "nonparametric/#HypothesisTests.SignedRankTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.SignedRankTest",
    "category": "function",
    "text": "SignedRankTest(x::AbstractVector{<:Real})\nSignedRankTest(x::AbstractVector{<:Real}, y::AbstractVector{T<:Real})\n\nPerform a Wilcoxon signed rank test of the null hypothesis that the distribution of x (or the difference x - y if y is provided) has zero median against the alternative hypothesis that the median is non-zero.\n\nWhen there are no tied ranks and ≤50 samples, or tied ranks and ≤15 samples, SignedRankTest performs an exact signed rank test. In all other cases, SignedRankTest performs an approximate signed rank test. Behavior may be further controlled by using ExactSignedRankTest or ApproximateSignedRankTest directly.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ExactSignedRankTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ExactSignedRankTest",
    "category": "type",
    "text": "ExactSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])\n\nPerform a Wilcoxon exact signed rank U test of the null hypothesis that the distribution of x (or the difference x - y if y is provided) has zero median against the alternative hypothesis that the median is non-zero.\n\nWhen there are no tied ranks, the exact p-value is computed using the psignrank function from the Rmath package. In the presence of tied ranks, a p-value is computed by exhaustive enumeration of permutations, which can be very slow for even moderately sized data sets.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ApproximateSignedRankTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ApproximateSignedRankTest",
    "category": "type",
    "text": "ApproximateSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])\n\nPerform a Wilcoxon approximate signed rank U test of the null hypothesis that the distribution of x (or the difference x - y if y is provided) has zero median against the alternative hypothesis that the median is non-zero.\n\nThe p-value is computed using a normal approximation to the distribution of the signed rank statistic:\n\n    beginalign*\n        μ  = fracn(n + 1)4\n        σ  = fracn(n + 1)(2 * n + 1)24 - fraca48\n        a  = sum_t in mathcalT t^3 - t\n    endalign*\n\nwhere mathcalT is the set of the counts of tied values at each tied position.\n\nImplements: pvalue, confint\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Wilcoxon-signed-rank-test-1",
    "page": "Nonparametric tests",
    "title": "Wilcoxon signed rank test",
    "category": "section",
    "text": "SignedRankTest\nExactSignedRankTest\nApproximateSignedRankTest"
},

{
    "location": "nonparametric/#HypothesisTests.ExactPermutationTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ExactPermutationTest",
    "category": "function",
    "text": "ExactPermutationTest(x::Vector, y::Vector, f::Function)\n\nPerform a permutation test (a.k.a. randomization test) of the null hypothesis that f(x) is equal to f(y).  All possible permutations are sampled.\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#HypothesisTests.ApproximatePermutationTest",
    "page": "Nonparametric tests",
    "title": "HypothesisTests.ApproximatePermutationTest",
    "category": "function",
    "text": "ApproximatePermutationTest(x::Vector, y::Vector, f::Function, n::Int)\n\nPerform a permutation test (a.k.a. randomization test) of the null hypothesis that f(x) is equal to f(y).  n of the factorial(length(x)+length(y)) permutations are sampled at random.\n\n\n\n\n\n"
},

{
    "location": "nonparametric/#Permutation-test-1",
    "page": "Nonparametric tests",
    "title": "Permutation test",
    "category": "section",
    "text": "ExactPermutationTest\nApproximatePermutationTest"
},

{
    "location": "time_series/#",
    "page": "Time series tests",
    "title": "Time series tests",
    "category": "page",
    "text": ""
},

{
    "location": "time_series/#Time-series-tests-1",
    "page": "Time series tests",
    "title": "Time series tests",
    "category": "section",
    "text": ""
},

{
    "location": "time_series/#HypothesisTests.DurbinWatsonTest",
    "page": "Time series tests",
    "title": "HypothesisTests.DurbinWatsonTest",
    "category": "type",
    "text": "DurbinWatsonTest(X::AbstractArray, e::AbstractVector; p_compute::Symbol = :ndep)\n\nCompute the Durbin-Watson test for serial correlation in the residuals of a regression model.\n\nX is the matrix of regressors from the original regression model and e the vector of residuals. Note that the Durbin-Watson test is not valid if X includes a lagged dependent variable. The test statistic is computed as\n\nDW = fracsum_t=2^n (e_t - e_t-1)^2sum_t=1^n e_t^2\n\nwhere n is the number of observations.\n\nBy default, the choice of approach to compute p-values depends on the sample size (p_compute=:ndep). For small samples (n<100), Pan\'s algorithm (Farebrother, 1980) is employed. For larger samples, a normal approximation is used (Durbin and Watson, 1950). To always use Pan\'s algorithm, set p_compute=:exact. p_compute=:approx will always use the normal approximation.\n\nDefault is a two-sided p-value for the alternative hypothesis of positive or negative serial correlation. One-sided p-values can be requested by calling pvalue(x::DurbinWatsonTest; tail=) with the options :left (negative serial correlation) and :right (positive serial correlation).\n\nReferences\n\nJ. Durbin and G. S. Watson, 1951, \"Testing for Serial Correlation in Least Squares Regression: II\", Biometrika, Vol. 38, No. 1/2, pp. 159-177, http://www.jstor.org/stable/2332325.\nJ. Durbin and G. S. Watson, 1950, \"Testing for Serial Correlation in Least Squares Regression: I\", Biometrika, Vol. 37, No. 3/4, pp. 409-428, http://www.jstor.org/stable/2332391.\nR. W. Farebrother, 1980, \"Algorithm AS 153: Pan\'s Procedure for the Tail Probabilities of the Durbin-Watson Statistic\", Journal of the Royal Statistical Society, Series C (Applied Statistics), Vol. 29, No. 2, pp. 224-227, http://www.jstor.org/stable/2986316.\n\nExternal links\n\nDurbin-Watson test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#Durbin-Watson-test-1",
    "page": "Time series tests",
    "title": "Durbin-Watson test",
    "category": "section",
    "text": "DurbinWatsonTest"
},

{
    "location": "time_series/#HypothesisTests.BoxPierceTest",
    "page": "Time series tests",
    "title": "HypothesisTests.BoxPierceTest",
    "category": "type",
    "text": "BoxPierceTest(y, lag, dof=0)\n\nCompute the Box-Pierce Q statistic to test the null hypothesis of independence in a time series y.\n\nlag specifies the number of lags used in the construction of Q. When testing the residuals of an estimated model, dof has to be set to the number of estimated parameters. E.g., when testing the residuals of an ARIMA(p,0,q) model, set dof=p+q.\n\nExternal links\n\nBox-Pierce test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#HypothesisTests.LjungBoxTest",
    "page": "Time series tests",
    "title": "HypothesisTests.LjungBoxTest",
    "category": "type",
    "text": "LjungBoxTest(y, lag, dof=0)\n\nCompute the Ljung-Box Q statistic to test the null hypothesis of independence in a time series y.\n\nlag specifies the number of lags used in the construction of Q. When testing the residuals of an estimated model, dof has to be set to the number of estimated parameters. E.g., when testing the residuals of an ARIMA(p,0,q) model, set dof=p+q.\n\nExternal links\n\nLjung-Box test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#Box-Pierce-and-Ljung-Box-tests-1",
    "page": "Time series tests",
    "title": "Box-Pierce and Ljung-Box tests",
    "category": "section",
    "text": "BoxPierceTest\nLjungBoxTest"
},

{
    "location": "time_series/#HypothesisTests.BreuschGodfreyTest",
    "page": "Time series tests",
    "title": "HypothesisTests.BreuschGodfreyTest",
    "category": "type",
    "text": "BreuschGodfreyTest(X, e, lag, start0 = true)\n\nCompute the Breusch-Godfrey test for serial correlation in the residuals of a regression model.\n\nX is the matrix of regressors from the original model and e the vector of residuals. lag determines the number of lagged residuals included in the auxiliary regression. Set start0 to specify how the starting values for the lagged residuals are handled. start0 = true (default) sets them to zero (as in Godfrey, 1978); start0 = false uses the first lag residuals as starting values, i.e. shortening the sample by lag.\n\nExternal links\n\nBreusch-Godfrey test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#Breusch-Godfrey-test-1",
    "page": "Time series tests",
    "title": "Breusch-Godfrey test",
    "category": "section",
    "text": "BreuschGodfreyTest"
},

{
    "location": "time_series/#HypothesisTests.JarqueBeraTest",
    "page": "Time series tests",
    "title": "HypothesisTests.JarqueBeraTest",
    "category": "type",
    "text": "JarqueBeraTest(y::AbstractVector)\n\nCompute the Jarque-Bera statistic to test the null hypothesis that a real-valued vector y is normally distributed.\n\nNote that the approximation by the Chi-squared distribution does not work well and the speed of convergence is slow. In small samples, the test tends to be over-sized for nominal levels up to about 3% and under-sized for larger nominal levels (Mantalos, 2010).\n\nReferences\n\nPanagiotis Mantalos, 2011, \"The three different measures of the sample skewness and kurtosis and the effects to the Jarque-Bera test for normality\", International Journal of Computational Economics and Econometrics, Vol. 2, No. 1, link.\n\nExternal links\n\nJarque-Bera test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#Jarque-Bera-test-1",
    "page": "Time series tests",
    "title": "Jarque-Bera test",
    "category": "section",
    "text": "JarqueBeraTest"
},

{
    "location": "time_series/#HypothesisTests.ADFTest",
    "page": "Time series tests",
    "title": "HypothesisTests.ADFTest",
    "category": "type",
    "text": "ADFTest(y, deterministic, lag)\n\nCompute the augmented Dickey-Fuller unit root test.\n\ny is the time series to be tested, deterministic determines the deterministic terms (options: none, constant, trend, squared_trend) and lag the number of lagged first-differences included in the test regression, respectively.\n\nCritical values and asymptotic p-values are computed based on response surface regressions following MacKinnon (2010) and MacKinnon (1994), respectively. These may differ slightly from those reported in other regression packages as different algorithms might be used.\n\nReferences\n\nJames G. MacKinnon, 2010, \"Critical values for cointegration tests,\" QED Working Paper No. 1227, 2010, link.\nJames G. MacKinnon, 1994, \"Approximate Asymptotic Distribution Functions for Unit-Root and Cointegration Tests\", Journal of Business & Economic Statistics, Vol. 12, No. 2, pp. 167-176, link.\n\nExternal links\n\nAugmented Dickey-Fuller test on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "time_series/#Augmented-Dickey-Fuller-test-1",
    "page": "Time series tests",
    "title": "Augmented Dickey-Fuller test",
    "category": "section",
    "text": "ADFTest"
},

]}
