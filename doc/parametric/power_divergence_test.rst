Power Divergence Test
=============================================

.. function:: PowerDivergenceTest(x [,y] [, lambda] [,theta0] )

    If ``x`` is a matrix with one row or column, or if ``x`` is a vector and ``y`` is
    not given, then a goodness-of-fit test is performed (``x`` is treated as a one-
    dimensional contingency table. The entries of ``x`` must be non-negative integers.
    In this case, the hypothesis tested is whether the population probabilities equal
    those in ``theta0``, or are all equal if ``theta0`` is not given.

    If ``x`` is a matrix with at least two rows and columns, it is taken as a two-dimensional
    contingency table: the entries of ``x`` must be non-negative integers. Otherwise, ``x``
    and ``y`` must be vectors of the same length. The contingency table is calculated using
    ``counts`` from ``Statsbase``. Then the power divergence test is performed of the null
    hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency
    table is the product of the row and column marginals.

    The power divergence test is given by

    .. math::

      \dfrac{2}{\lambda(\lambda+1)}\sum_{i=1}^I \sum_{j=1}^J n_{ij}\left[(n_{ij}/\hat{n}_{ij})^\lambda -1\right]

    where :math:`n_{ij}` is the cell count in the :math:`{i}` the row and :math:`{j}` the column and :math:`\lambda` is a real number.
    Note that when :math:`\lambda = 1`, this is equal to Pearson's chi-squared statistic, as :math:`\lambda \to 0`, it converges
    to the likelihood ratio test statistic, as :math:`\lambda \to -1` it converges to the minimum discrimination information
    statistic (Gokhale and Kullback 1978), for :math:`\lambda=-2` it equals Neyman modified chi-squared (Neyman 1949), and for
    :math:`\lambda=-1/2` it equals the Freeman-Tukey statistic (Freeman and Tukey 1950). Under regularity conditions, their
    asymptotic distributions are identical (see Drost et. al. 1989). The chi-squared null approximation works best for
    :math:`\lambda` near :math:`{2/3}`.

    Implements: :ref:`pvalue<pvalue>`, :ref:`confint<confint>`

    References:

    - Agresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013.


.. function:: ChisqTest(x [,y] [,theta0])

    Convenience function for power divergence test with :math:`\lambda=1`.

    .. code-block:: none

        julia> using FreqTables

        julia> using HypothesisTests

        julia> using RDatasets

        julia> data = dataset("MASS", "survey");

        julia> table = freqtable(data, :Smoke, :Exer)
        4×3 Named Array{Int64,2}
        Smoke ╲ Exer │ Freq  None  Some
        ─────────────┼─────────────────
        Heavy        │    7     1     3
        Never        │   87    18    84
        Occas        │   12     3     4
        Regul        │    9     1     7

        julia> ChisqTest(table)
        Pearson's Chi-square Test
        -------------------------
        Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.0227126,0.390243,0.0392308,0.0351013,0.00454252,0.0780487,0.00784616,0.00702025,0.0193551,0.332555,0.0334315,0.0299124]
        point estimate:          [0.029661,0.368644,0.0508475,0.0381356,0.00423729,0.0762712,0.0127119,0.00423729,0.0127119,0.355932,0.0169492,0.029661]
        95% confidence interval: Tuple{Float64,Float64}[(0.0,0.177966),(0.211864,0.516949),(0.0,0.199152),(0.0,0.186441),(0.0,0.152542),(0.0,0.224576),(0.0,0.161017),(0.0,0.152542),(0.0,0.161017),(0.199153,0.504237),(0.0,0.165254),(0.0,0.177966)]

        Test summary:
        outcome with 95% confidence: fail to reject h_0
        two-sided p-value:           0.48284216946545644 (not significant)

        Details:
        Sample size:        236
        statistic:          5.48854589058423
        degrees of freedom: 6
        residuals:          [0.708288,-0.531165,0.900995,0.248804,-0.0695717,-0.0977427,0.843864,-0.510255,-0.733561,0.622746,-1.38483,-0.0223272]
        std. residuals:     [1.01307,-1.66226,1.31224,0.360707,-0.0750004,-0.230546,0.926328,-0.557555,-0.982466,1.82488,-1.8886,-0.0303099]

.. function:: MultinomialLRT(x [,y] [,theta0])

    Convenience function for power divergence test with :math:`\lambda=0`.

    .. code-block:: none

        julia> using FreqTables

        julia> using HypothesisTests

        julia> using RDatasets

        julia> data = dataset("MASS", "survey");

        julia> table = freqtable(data, :Smoke, :Exer)
        4×3 Named Array{Int64,2}
        Smoke ╲ Exer │ Freq  None  Some
        ─────────────┼─────────────────
        Heavy        │    7     1     3
        Never        │   87    18    84
        Occas        │   12     3     4
        Regul        │    9     1     7

        julia> MultinomialLRT(table)
        Multinomial Likelihood Ratio Test
        ---------------------------------
        Population details:
        parameter of interest:   Multinomial Probabilities
        value under h_0:         [0.0227126,0.390243,0.0392308,0.0351013,0.00454252,0.0780487,0.00784616,0.00702025,0.0193551,0.332555,0.0334315,0.0299124]
        point estimate:          [0.029661,0.368644,0.0508475,0.0381356,0.00423729,0.0762712,0.0127119,0.00423729,0.0127119,0.355932,0.0169492,0.029661]
        95% confidence interval: Tuple{Float64,Float64}[(0.0,0.177966),(0.211864,0.516949),(0.0,0.199152),(0.0,0.186441),(0.0,0.152542),(0.0,0.224576),(0.0,0.161017),(0.0,0.152542),(0.0,0.161017),(0.199153,0.504237),(0.0,0.165254),(0.0,0.177966)]

        Test summary:
        outcome with 95% confidence: fail to reject h_0
        two-sided p-value:           0.4457935354971526 (not significant)

        Details:
        Sample size:        236
        statistic:          5.8014667453476925
        degrees of freedom: 6
        residuals:          [0.708288,-0.531165,0.900995,0.248804,-0.0695717,-0.0977427,0.843864,-0.510255,-0.733561,0.622746,-1.38483,-0.0223272]
        std. residuals:     [1.01307,-1.66226,1.31224,0.360707,-0.0750004,-0.230546,0.926328,-0.557555,-0.982466,1.82488,-1.8886,-0.0303099]
