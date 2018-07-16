using HypothesisTests, Compat.Test

@testset "Show" begin
# based on power_divergence.jl tests
d = [[762,484] [327,239] [468,477]]
m = PowerDivergenceTest(d)

@test sprint(show, m) ==
    """
    Pearson's Chi-square Test
    -------------------------
    Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.255231, 0.19671, 0.11594, 0.0893561, 0.193574, 0.14919]
        point estimate:          [0.276387, 0.175553, 0.118607, 0.0866884, 0.16975, 0.173014]
        95% confidence interval: Tuple{Float64,Float64}[(0.2579, 0.2955), (0.1571, 0.1947), (0.1001, 0.1378), (0.0682, 0.1058), (0.1513, 0.1889), (0.1545, 0.1922)]

    Test summary:
        outcome with 95% confidence: reject h_0
        one-sided p-value:           <1e-6

    Details:
        Sample size:        2757
        statistic:          30.070149095754687
        degrees of freedom: 2
        residuals:          [2.19886, -2.50467, 0.41137, -0.468583, -2.84324, 3.23867]
        std. residuals:     [4.50205, -4.50205, 0.699452, -0.699452, -5.31595, 5.31595]
    """

d = [ 20, 15, 25 ]
m = PowerDivergenceTest(d)

@test sprint(show, m) ==
    """
    Pearson's Chi-square Test
    -------------------------
    Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.333333, 0.333333, 0.333333]
        point estimate:          [0.333333, 0.25, 0.416667]
        95% confidence interval: Tuple{Float64,Float64}[(0.2167, 0.481), (0.1333, 0.3976), (0.3, 0.5643)]

    Test summary:
        outcome with 95% confidence: fail to reject h_0
        one-sided p-value:           0.2865

    Details:
        Sample size:        60
        statistic:          2.5
        degrees of freedom: 2
        residuals:          [0.0, -1.11803, 1.11803]
        std. residuals:     [0.0, -1.36931, 1.36931]
    """

# based on t.jl tests
tst = OneSampleTTest(-5:10)

@test sprint(show, tst) ==
    """
    One sample t-test
    -----------------
    Population details:
        parameter of interest:   Mean
        value under h_0:         0
        point estimate:          2.5
        95% confidence interval: (-0.0369, 5.0369)

    Test summary:
        outcome with 95% confidence: fail to reject h_0
        two-sided p-value:           0.0530

    Details:
        number of observations:   16
        t-statistic:              2.100420126042015
        degrees of freedom:       15
        empirical standard error: 1.1902380714238083
    """

end
