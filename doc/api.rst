Methods
===========

.. function:: pvalue(test::HypothesisTest; tail=:both)

    Compute the p-value for a given significance test.

    If ``tail`` is ``:both`` (default), then the p-value for the
    two-sided test is returned. If ``tail`` is ``:left`` or
    ``:right``, then a one-sided test is performed.

.. function:: ci(test::HypothesisTest, alpha=0.05; tail=:both)

    Compute a confidence interval with coverage 1-``alpha``.

    If ``tail`` is ``:both`` (default), then a two-sided confidence
    interval is returned. If ``tail`` is ``:left`` or
    ``:right``, then a one-sided confidence interval is returned
