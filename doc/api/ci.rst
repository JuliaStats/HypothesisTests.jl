.. _ci:

Confidence Interval
==============================================

.. function:: ci(test::HypothesisTest, alpha=0.05; tail=:both)

    Compute a confidence interval C with coverage 1-``alpha``.

    If ``tail`` is ``:both`` (default), then a two-sided confidence
    interval is returned. If ``tail`` is ``:left`` or
    ``:right``, then a one-sided confidence interval is returned

    .. note::
      Most of the implemented confidence intervals are *strongly consistent*, that is, the
      confidence interval with coverage 1-``alpha`` does not contain the test statistic
      under :math:`h_0` if and only if the corresponding test rejects the null hypothesis
      :math:`h_0: \theta=\theta_0`:

      .. math::
        C (x, 1 − \alpha) = \{\theta : p_\theta (x) > \alpha\},

      where :math:`p_\theta` is the :ref:`p-value<pvalue>` of the corresponding test.


.. _ci_binomial:

Confidence Interval for Binomial Proportions
----------------------------------------------

.. function:: ci(test::BinomialTest, alpha=0.05; tail=:both, method=:clopper_pearson)

   Compute a confidence interval with coverage 1-``alpha`` for a binomial proportion 
   using one of the following methods. Possible values for ``method`` are:

   - Clopper-Pearson interval ``:clopper_pearson`` (default): This interval is  based on the binomial 
     distribution. The empirical coverage is never less than the nominal coverage of 
     1-``alpha``; it is usually too conservative.
   - Wald interval ``:wald`` (normal approximation interval): This interval relies on 
     the standard approximation of the actual binomial distribution by a normal distribution. 
     Coverage can be erratically poor for success probabilities close to zero or one. 
   - Wilson score interval ``:wilson``: This interval relies on a normal approximation. 
     In contrast to ``:wald`` the standard deviation is not approximated by an empirical
     estimate resulting in good empirical coverages even for small numbers of draws and 
     extreme success probabilities.
   - Jeffreys interval ``:jeffrey``: Bayesian confidence interval obtained by using a
     non-informative Jeffreys prior. The interval is very similar to the Wilson interval. 
   - Agresti Coull interval ``:agresti_coull``: Simplified version of the Wilson interval;
     they are centered around the same value. The Agresti Coull interval has higher or 
     equal coverage.

   References:

   - Brown, L.D., Cai, T.T., and DasGupta, A. Interval estimation for a binomial proportion. 
     Statistical Science, 16(2):101–117, 2001. 
