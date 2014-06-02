Sign test
=============================================

.. function:: SignTest(x::AbstractVector{T<:Real}, median::Real=0)

	Perform a sign test of the null hypothesis that the distribution
	from which ``x`` was drawn has median ``median`` against the
	alternative hypothesis that the median is not equal to ``median``.

    Implements: :ref:`pvalue<pvalue>`

.. function:: SignTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, median::Real=0)

	Perform a sign test of the null hypothesis that the distribution
	of ``x - y`` has median ``median`` against the alternative
	hypothesis that the median is not equal to ``median``.

    Implements: :ref:`pvalue<pvalue>`
