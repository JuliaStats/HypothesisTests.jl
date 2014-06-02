p-value
==============================================

.. _pvalue:

.. function:: pvalue(test::HypothesisTest; tail=:both)

    Compute the p-value for a given significance test.

    If ``tail`` is ``:both`` (default), then the p-value for the
    two-sided test is returned. If ``tail`` is ``:left`` or
    ``:right``, then a one-sided test is performed.