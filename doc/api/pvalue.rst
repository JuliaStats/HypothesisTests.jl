.. _pvalue:

p-value
==============================================

.. function:: pvalue(test::HypothesisTest; tail=tail(test))

    Compute the p-value for a given significance test.

    If ``tail`` is ``:both`` (default for most hypothesis tests), then the p-value for the
    two-sided test is returned. If ``tail`` is ``:left`` or
    ``:right``, then a one-sided test is performed.

.. _pvalue_fisher:

p-value for Fisher exact test
----------------------------------------------

.. function:: function pvalue(x::FisherExactTest; tail=:both, method=:central)

   Compute the p-value for a given significance test. The one-sided p-values are based on
   Fisher's non-central hypergeometric distribution :math:`f_\omega(i)` with odd-ratio :math:`\omega`:

   .. math::
        p_\omega^{(\text{left})} &=\sum_{i\leq a} f_\omega(i)\\
        p_\omega^{(\text{right})} &=\sum_{i\geq a} f_\omega(i)

   For ``tail=:both``, possible values for ``method`` are:

   - Central interval ``:central`` (default): This p-value is two times the minimum of the one-sided
     p-values.

   - Minimum likelihood interval ``:minlike``: This p-value is computed by summing all tables with the same
     marginals that are equally or less probable:

    .. math::
        p_\omega &=\sum_{f_\omega(i)\leq f_\omega(a)} f_\omega(i)

    .. note::
        Since the p-value is not necessarily unimodal, the corresponding confidence region might not be an interval.

   References:

   - Gibbons, J.D, Pratt, J.W. P-values: Interpretation and Methodology
     American Statistican, 29(1):20-25, 1975.
   - Fay, M.P. Supplementary material to Confidence intervals that match Fisher’s exact or Blaker’s exact tests.
     Biostatistics, 0(0):1-13, 2009.
