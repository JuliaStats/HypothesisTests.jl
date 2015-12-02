Anderson–Darling test
=============================================

The null hypothesis of the Anderson–Darling test is that a dataset comes from
a certain distribution; the reference distribution can be specified explicitely 
(one-sample test). K-sample Anderson–Darling tests are available for testing whether
several samples are coming from a single population drawn from the distribution function
which does not have to be specified.

.. function:: OneSampleADTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)

    Perform a one sample Anderson–Darling test of the null hypothesis that the data
    in vector ``x`` comes from the distribution ``d`` against
    the alternative hypothesis that the sample is not drawn from ``d``.

    Implements: :ref:`pvalue<pvalue>`

.. function:: KSampleADTest{T<:Real}(xs::AbstractVector{T}...; modified=true)

    Perform an k-sample Anderson–Darling test of the null hypothesis that the data
    in vectors ``xs`` comes from the same distribution against the alternative
    hypothesis that the samples comes from different distributions.

    ``modified`` paramater enables a modified test calculation for samples whose
    observations do not all coincide.

    Implements: :ref:`pvalue<pvalue>`

    References:

    - k-Sample Anderson-Darling Tests, F. W. Scholz and M. A. Stephens, Journal of the American Statistical Association, Vol. 82, No. 399. (Sep., 1987), pp. 918-924.
