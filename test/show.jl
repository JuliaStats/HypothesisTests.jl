using HypothesisTests, Test

@testset "Show" begin
# based on power_divergence.jl tests
d = [[762,484] [327,239] [468,477]]
m = PowerDivergenceTest(d)

if VERSION < v"1.4"
    @test sprint(show, m, context=:compact => true) ==
    """
    Pearson's Chi-square Test
    -------------------------
    Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.255231, 0.19671, 0.11594, 0.0893561, 0.193574, 0.14919]
        point estimate:          [0.276387, 0.175553, 0.118607, 0.0866884, 0.16975, 0.173014]
        95% confidence interval: Tuple{Float64,Float64}[(0.2545, 0.2994), (0.1573, 0.1955), (0.1033, 0.1358), (0.07357, 0.1019), (0.1517, 0.1894), (0.1548, 0.1928)]

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

    @test sprint(show, m, context=:compact => true) ==
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
elseif VERSION < v"1.6"
    @test sprint(show, m, context=:compact => true) ==
    """
    Pearson's Chi-square Test
    -------------------------
    Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.255231, 0.19671, 0.11594, 0.0893561, 0.193574, 0.14919]
        point estimate:          [0.276387, 0.175553, 0.118607, 0.0866884, 0.16975, 0.173014]
        95% confidence interval: [(0.2545, 0.2994), (0.1573, 0.1955), (0.1033, 0.1358), (0.07357, 0.1019), (0.1517, 0.1894), (0.1548, 0.1928)]

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

    @test sprint(show, m, context=:compact => true) ==
    """
    Pearson's Chi-square Test
    -------------------------
    Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.333333, 0.333333, 0.333333]
        point estimate:          [0.333333, 0.25, 0.416667]
        95% confidence interval: [(0.2167, 0.481), (0.1333, 0.3976), (0.3, 0.5643)]

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
else
    @test sprint(show, m, context=:compact => true) ==
        """
        Pearson's Chi-square Test
        -------------------------
        Population details:
            parameter of interest:   Multinomial Probabilities
            value under h_0:         [0.255231, 0.19671, 0.11594, 0.0893561, 0.193574, 0.14919]
            point estimate:          [0.276387, 0.175553, 0.118607, 0.0866884, 0.16975, 0.173014]
            95% confidence interval: [(0.2545, 0.2994), (0.1573, 0.1955), (0.1033, 0.1358), (0.07357, 0.1019), (0.1517, 0.1894), (0.1548, 0.1928)]

        Test summary:
            outcome with 95% confidence: reject h_0
            one-sided p-value:           <1e-06

        Details:
            Sample size:        2757
            statistic:          30.070149095754687
            degrees of freedom: 2
            residuals:          [2.19886, -2.50467, 0.41137, -0.468583, -2.84324, 3.23867]
            std. residuals:     [4.50205, -4.50205, 0.699452, -0.699452, -5.31595, 5.31595]
        """

    d = [ 20, 15, 25 ]
    m = PowerDivergenceTest(d)

    @test sprint(show, m, context=:compact => true) ==
        """
        Pearson's Chi-square Test
        -------------------------
        Population details:
            parameter of interest:   Multinomial Probabilities
            value under h_0:         [0.333333, 0.333333, 0.333333]
            point estimate:          [0.333333, 0.25, 0.416667]
            95% confidence interval: [(0.2167, 0.481), (0.1333, 0.3976), (0.3, 0.5643)]

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
end

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
        95% confidence interval: (-0.03693, 5.037)

    Test summary:
        outcome with 95% confidence: fail to reject h_0
        two-sided p-value:           0.0530

    Details:
        number of observations:   16
        t-statistic:              2.100420126042015
        degrees of freedom:       15
        empirical standard error: 1.1902380714238083
    """

# issue #248
x = [
    6.598170000000001e-7,
    6.9452e-7,
    2.41933e-7,
    2.4264999999999997e-7,
    4.830650000000001e-7,
    2.67262e-7,
    2.5027699999999996e-7,
    2.51241e-7,
    2.67511e-7,
    2.27148e-7,
    2.41169e-7,
    2.56646e-7,
    4.31067e-7,
    2.2686500000000001e-7,
    2.35553e-7,
    2.32062e-7,
    7.36284e-7,
]
y = [
    3.59147e-7,
    2.75594e-7,
    1.63942e-7,
    1.7980399999999999e-7,
    2.82113e-7,
    1.6574299999999998e-7,
    1.60492e-7,
    1.61266e-7,
    2.12196e-7,
    1.51524e-7,
    1.86578e-7,
    2.1346e-7,
    1.59902e-7,
    1.50073e-7,
    1.64e-7,
    1.42769e-7,
    1.70032e-7,
]
tst = UnequalVarianceTTest(x, y)

@test sprint(show, tst) ==
    """
    Two sample t-test (unequal variance)
    ------------------------------------
    Population details:
        parameter of interest:   Mean difference
        value under h_0:         0
        point estimate:          1.55673e-7
        95% confidence interval: (5.93e-8, 2.52e-7)

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           0.0031

    Details:
        number of observations:   [17,17]
        t-statistic:              3.3767280623082523
        degrees of freedom:       19.363987783845342
        empirical standard error: 4.610162387563106e-8
    """
end
