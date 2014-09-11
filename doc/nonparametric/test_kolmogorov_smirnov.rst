Kolmogorov–Smirnov test
=============================================

The null hypothesis of the Kolmogorov–Smirnov test is that a dataset comes from
a certain distribution; the reference distribution can be specified explicitely 
(one-sample test) or by an empirical sample (two-sample test). The alternative
hypothesis is that the cumulative distributions of the sample is different
(``tail=:both``; default), smaller (``tail=:left`), or larger (``tail:=right``)
than the reference cumulative distribution. The exact test is based on the exact
distribution of the differences whereas the approximate test is derived from its
asymptotic distribution.

.. function:: ExactOneSampleKSTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)

    Perform a one sample Kolmogorov–Smirnov-test of the null hypothesis that the data
    in vector ``x`` comes from the distribution ``d`` against
    the alternative hypothesis that the sample is not drawn from ``d``.

    Implements: :ref:`pvalue<pvalue>`

.. function:: ApproximateOneSampleKSTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)

    Perform an asymptotic one sample Kolmogorov–Smirnov-test of the null hypothesis that the data
    in vector ``x`` comes from the distribution ``d`` against
    the alternative hypothesis that the sample is not drawn from ``d``.

    Implements: :ref:`pvalue<pvalue>`

.. function:: ApproximateTwoSampleKSTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})

    Perform an asymptotic two sample Kolmogorov–Smirnov-test of the null hypothesis that ``x``
    and ``y`` are drawn from the same distribution against
    the alternative hypothesis that the distribution comes from different
    distributions. 

    Implements: :ref:`pvalue<pvalue>`

    References:

    - Approximation of one-sided test: http://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test