Binomial test
=============================================

.. function:: BinomialTest(x::Integer, n::Integer, p::Real=0.5)

	Perform a binomial test of the null hypothesis that the distribution
	from which ``x`` successes were encountered in ``n`` draws has
	success probability ``p`` against the alternative hypothesis that
	the success probability is not equal to ``p``.

	Computed confidence intervals are Clopper-Pearson intervals.

    Implements: :ref:`pvalue<pvalue>`, :ref:`confint<ci_binomial>`

.. function:: BinomialTest(x::AbstractVector{Bool}, p::Real=0.5)

	Perform a binomial test of the null hypothesis that the distribution
	from which ``x`` was drawn has success probability ``p`` against the
	alternative hypothesis that the success probability is not equal to
	``p``.

    Implements: :ref:`pvalue<pvalue>`, :ref:`confint<ci_binomial>`
