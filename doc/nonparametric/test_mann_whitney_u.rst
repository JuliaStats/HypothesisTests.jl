Mann Whitney U test
=============================================

.. function:: MannWhitneyUTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform a Mann-Whitney U test of the null hypothesis that the
    probability that an observation drawn from the same population as
    ``x`` is greater than an observation drawn from the same
    population as ``y`` is equal to the probability that an
    observation drawn from the same population as ``y`` is greater
    than an observation drawn from the same population as ``x``
    against the alternative hypothesis that these probabilities are not
    equal. 

    The Mann-Whitney U test is sometimes known as the Wilcoxon rank sum
    test.

    When there are no tied ranks and ≤50 samples, or tied ranks and
    ≤10 samples, ``MannWhitneyUTest`` performs an exact
    Mann-Whitney U test. In all other cases, ``MannWhitneyUTest``
    performs an approximate Mann-Whitney U test. Behavior may be
    further controlled by using ``ExactMannWhitneyUTest`` or
    ``ApproximateMannWhitneyUTest`` directly. See below for further
    algorithmic details.
    
    Implements: :ref:`pvalue<pvalue>`

.. function:: ExactMannWhitneyUTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform an exact Mann-Whitney U test.

    When there are no tied ranks, the exact p-value is computed using
    the ``pwilcox`` function from libRmath. In the presence of tied
    ranks, a p-value is computed by exhaustive enumeration of
    permutations, which can be very slow for even moderately sized
    data sets.

    Implements: :ref:`pvalue<pvalue>`

.. function:: ApproximateMannWhitneyUTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform an approximate Mann-Whitney U test.

    The p-value is computed using a normal approximation to the
    distribution of the Mann-Whitney *U* statistic:

    .. math::
        \mu &= \frac{n_x n_y}{2}\\
        \sigma &= \frac{n_x n_y}{12}\left(n_x + n_y + 1 - \frac{a}{(n_x + n_y)(n_x + n_y - 1)}\right)\\
        a &= \sum_{t \in \mathcal{T}} t^3 - t

    where :math:`\mathcal{T}` is the set of the counts of tied values
    at each tied position.

    Implements: :ref:`pvalue<pvalue>`
