Power Divergence Test
=============================================

.. function:: PowerDivergenceTest(x [,y] [, lambda] [,theta0] )

    If ``x`` is a matrix with one row or column, or if ``x`` is a vector and ``y`` is
    not given, then a goodness-of-fit test is performed (``x`` is treated as a one-
    dimensional contingency table. The entries of ``x`` must be non-negative integers. 
    In this case, the hypothesis tested is whether the population probabilities equal 
    those in ``theta0``, or are all equal if ``theta0`` is not given.

    If ``x`` is a matrix with at least two rows and columns, it is taken as a two-dimensional
    contigency table: the entries of ``x`` must be non-negative integers. Otherwise, ``x``
    and ``y`` must be vectors of the same length. The contigency table is calculated using
    ``counts`` from ``Statsbase``. Then the power divergence test is performed of the null
    hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency
    table is the product of the row and column marginals. 

    The power divergence test is given by 

    .. math::

      \dfrac{2}{\lambda(\lambda+1)}\sum_{i=1}^I \sum_{j=1}^J n_{ij}\left[(n_{ij}/\hat{n}_{ij})^\lambda -1\right]

    where :math:`n_{ij}` is the cell count in the :math:`{i}` th row and :math:`{j}` th column and :math:`\lambda` is a real number.
    Note that when :math:`\lambda = 1`, this is equal to Pearson's chi-squared statistic, as :math`\lambda \to 0`, it converges
    to the likelihood ratio test statistic, as :math:`\lambda \to -1` it converges to the minimum discrimination information 
    statistic (Gokhale and Kullback 1978), for :math:`\lambda=-2` it equals Neyman modified chi-squared (Neyman 1949), and for 
    :math:`\lambda=-1/2` it equals the Freeman-Tukey statistic (Freeman and Tukey 1950). Under regulairty conditions, their
    asymptotic distributions are identical (see Drost et. al. 1989). The chis-squared null approximation works best for 
    :math:`\lambda` near :math:`{2/3}`.  

    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

    References:

    - Agresti, Alan. Categorical Data Analysis, 3rd Edition. Wiley, 2013. 


.. function:: ChisqTest(x [,y] [,theta0])

    Convenience function for power divergence test with :math:`\lambda=1`. 

.. function:: MultinomialLRT(x [,y] [,theta0])

    Convenience function for power divergence test with :math:`\lambda=0`.
