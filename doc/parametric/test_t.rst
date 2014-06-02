T-test
=============================================

.. function:: OneSampleTTest(v::AbstractVector{T<:Real}, mu0::Real=0)

    Perform a one sample t-test of the null hypothesis that the data
    in vector ``v`` comes from a distribution with ``mu0`` against
    the alternative hypothesis that the distribution does not have mean
    ``mu0``.

    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

.. function:: OneSampleTTest(xbar::Real, stdev::Real, n::Int, mu0::Real=0)

    Perform a one sample t-test of the null hypothesis that ``n``
    values with mean ``xbar`` and sample standard deviation
    ``stdev``  come from a distribution with ``mu0`` against
    the alternative hypothesis that the distribution does not have mean
    ``mu0``. 
    
    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

.. function:: OneSampleTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, mu0::Real=0)

    Perform a paired sample t-test of the null hypothesis that
    the differences between pairs of values in vectors ``x`` and
    ``y`` come from a distribution with ``mu0`` against the
    alternative hypothesis that the distribution does not have mean
    ``mu0``.
    
    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

.. function:: EqualVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform a two-sample t-test of the null hypothesis that
    ``x`` and ``y`` come from a distributions with the same mean
    and equal variances against the alternative hypothesis that the
    distributions have different means and but equal variances.
    
    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

.. function:: UnequalVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform an unequal variance two-sample t-test of the null
    hypothesis that ``x`` and ``y`` come from a distributions with
    the same mean against the alternative hypothesis that the
    distributions have different means.

    This test is also known as sometimes known as Welch's t-test. It
    differs from the equal variance t-test in that it computes the
    number of degrees of freedom of the test using the
    Welch-Satterthwaite equation:

    .. math::
        \nu_{\chi'} \approx \frac{\left(\sum_{i=1}^n k_i s_i^2\right)^2}
                                 {\sum_{i=1}^n \frac{(k_i s_i^2)^2}{\nu_i}}
    
    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`
