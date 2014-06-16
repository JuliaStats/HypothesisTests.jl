Wilcoxon signed rank test
=============================================

.. function:: SignedRankTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform a Wilcoxon signed rank test of the null hypothesis that the
    distribution of the difference ``x - y`` has zero median against the
    alternative hypothesis that the median is non-zero.

    When there are no tied ranks and ≤50 samples, or tied ranks and ≤15
    samples, ``SignedRankTest`` performs an exact signed rank test. In
    all other cases, ``SignedRankTest`` performs an approximate signed
    rank test. Behavior may be further controlled by using
    ``ExactSignedRankTest`` or ``ApproximateSignedRankTest`` directly.
    See below for further algorithmic details.

    Implements: :ref:`pvalue<pvalue>`
    
.. function:: SignedRankTest(x::AbstractVector{T<:Real})

    Perform a Wilcoxon signed rank test of the null hypothesis that the
    distribution from which ```x``` is drawn has zero median against the
    alternative hypothesis that the median is non-zero.

    Implements: :ref:`pvalue<pvalue>`
    
.. function:: ExactSignedRankTest(x::AbstractVector{T<:Real}[, y::AbstractVector{T<:Real}])

    Perform an exact signed rank U test.

    When there are no tied ranks, the exact p-value is computed using
    the ``psignrank`` function from libRmath. In the presence of tied
    ranks, a p-value is computed by exhaustive enumeration of
    permutations, which can be very slow for even moderately sized
    data sets.

    Implements: :ref:`pvalue<pvalue>`

.. function:: ApproximateSignedRank(x::AbstractVector{T<:Real}[, y::AbstractVector{T<:Real}])

    Perform an approximate signed rank U test.

    The p-value is computed using a normal approximation to the
    distribution of the signed rank statistic:

    .. math::
        \mu &= \frac{n(n + 1)}{4}\\
        \sigma &= \frac{n(n + 1)(2 * n + 1)}{24} - \frac{a}{48}\\
        a &= \sum_{t \in \mathcal{T}} t^3 - t

    where :math:`\mathcal{T}` is the set of the counts of tied values
    at each tied position.

    Implements: :ref:`pvalue<pvalue>`
