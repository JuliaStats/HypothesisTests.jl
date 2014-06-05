Kruskal-Wallis rank sum test
=============================================

.. function:: KruskalWallisTest{T<:Real}(groups::AbstractVector{T}...)

    Perform Kruskal Wallis rank sum test of the null hypothesis that the location parameters of the 
    distribution of the :math:`n` observations are the same in each of the ``groups`` :math:`\mathcal{G}`
    against the alternative hypothesis that they differ in at least one.

    The Kruskal-Wallis test is an extension of the :ref:`Mann-Whitney U test<test_mann_whitney>` to 
    more than two groups. 

    The p-value is computed using a chi-square approximation to the distribution of the test statistic :math:`H_c=\frac H C`:

    .. math::
        H &= \frac{12}{n(n+1)} \sum_{g \in \mathcal{G}} \frac{R_g^2}{n_g} - 3(n+1)\\
        C &= 1-\frac{1}{n^3-n}\sum_{t \in \mathcal{T}} (t^3-t),

    where :math:`\mathcal{T}` is the set of the counts of tied values
    at each tied position, :math:`n_g` is the number of observations and 
    :math:`R_g` is the rank sum in group g.

    Implements: :ref:`pvalue<pvalue>`
