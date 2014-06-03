Fisher exact test
=============================================

.. function:: FisherExactTest(a::Integer, b::Integer, c::Integer, d::Integer)

	Perform Fisher's exact test of the null hypothesis that the
	success probabilities ``a/c`` and ``b/d`` are equal, that is the odds ratio
	``(a/c) / (b/d)`` is one, against the alternative hypothesis that they are
	not equal.

	The contingency table is structured as:

		.. list-table:: 
		   :widths: 10 10 10
		   :header-rows: 1
		   :stub-columns: 1

		   * - 
		     - X1
		     - X1
		   * - Y1
		     - a
		     - b
		   * - Y2
		     - c
		     - d


	Two-sided p-values are computed by summing all tables with the same
	marginals that are equally or less probable.

	Implements: :ref:`pvalue<pvalue>`

	References:

	- Fay, M.P. Supplementary material to Confidence intervals that match Fisher’s exact or Blaker’s exact tests.
	  Biostatistics, 0(0): 1-13, 2009.
