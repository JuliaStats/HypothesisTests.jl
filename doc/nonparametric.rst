Nonparametric tests
=============================================

.. function:: BinomialTest(x::Integer, n::Integer, p::Real=0.5)

	Perform a binomial test of the null hypothesis that the distribution
	from which ``x`` successes were encountered in ``n`` draws has
	success probability ``p`` against the alternative hypothesis that
	the success probability is not equal to ``p``.

	Computed confidence intervals are Clopper-Pearson intervals.

    Implements: ``pvalue``, ``ci``

.. function:: BinomialTest(x::AbstractVector{Bool}, p::Real=0.5)

	Perform a binomial test of the null hypothesis that the distribution
	from which ``x`` was drawn has success probability ``p``` against the
	alternative hypothesis that the success probability is not equal to
	``p``.

    Implements: ``pvalue``, ``ci``

.. function:: SignTest(x::AbstractVector{T<:Real}, median::Real=0)

	Perform a sign test of the null hypothesis that the distribution
	from which ``x`` was drawn has median ``median`` against the
	alternative hypothesis that the median is not equal to ``median``.

    Implements: ``pvalue``

.. function:: SignTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, median::Real=0)

	Perform a sign test of the null hypothesis that the distribution
	of ``x - y`` has median ``median`` against the alternative
	hypothesis that the median is not equal to ``median``.

    Implements: ``pvalue``

.. function:: FisherExactTest(a::Integer, b::Integer, c::Integer, d::Integer)

	Perform Fisher's exact test of the null hypothesis that the
	success probabilities ``a/c`` and ``b/d`` are equal against the
	alternative hypothesis that they are not equal.

	The contingency table is structured as:
	+----+----+----+
	|    | X1 | X2 |
	+====+====+====+
	| Y1 |  a |  b |
    +====+----+----+
	| Y2 |  c |  d |
    +====+----+----+

    Two-sided p-values are computed by summing all tables with the same
    marginals that are equally or less probable.

    Implements: ``pvalue``

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
    
    Implements: ``pvalue``

.. function:: ExactMannWhitneyUTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

    Perform an exact Mann-Whitney U test.

    When there are no tied ranks, the exact p-value is computed using
    the ``pwilcox`` function from libRmath. In the presence of tied
    ranks, a p-value is computed by exhaustive enumeration of
    permutations, which can be very slow for even moderately sized
    data sets.

    Implements: ``pvalue``

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

    Implements: ``pvalue``

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

    Implements: ``pvalue``
    
.. function:: SignedRankTest(x::AbstractVector{T<:Real})

    Perform a Wilcoxon signed rank test of the null hypothesis that the
    distribution from which ```x``` is drawn has zero median against the
    alternative hypothesis that the median is non-zero.

    Implements: ``pvalue``
    
.. function:: ExactSignedRankTest(x::AbstractVector{T<:Real}[, y::AbstractVector{T<:Real}])

    Perform an exact signed rank U test.

    When there are no tied ranks, the exact p-value is computed using
    the ``psignrank`` function from libRmath. In the presence of tied
    ranks, a p-value is computed by exhaustive enumeration of
    permutations, which can be very slow for even moderately sized
    data sets.

    Implements: ``pvalue``

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

    Implements: ``pvalue``
