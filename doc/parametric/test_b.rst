B-test
=============================================

.. function:: BTest{T<:Real, S<:Real}(X::AbstractMatrix{T}, Y::AbstractMatrix{S}; kernel::Symbol, blocksize::Int)

    Perform a two-sample Block-test of the null hypothesis that
    :math:`{\bf X}=({\bf x}_1,\dots,{\bf x}_n)^T\in \mathbb{R}^{n\times m}` and
    :math:`{\bf Y}=({\bf y}_1,\dots,{\bf y}_n)^T\in \mathbb{R}^{n\times m}` come
    from the same distributions against the alternative hypothesis that the
    distributions are different. The difference of the distributions is measured
    by the maximum mean discrepency (MMD), that is, the the squared distance
    between the means of the distributions in the reproducing kernel Hilbert
    space induced by the ``kernel`` function. Possible values are:

    - Radial basis function ``:rbf`` (default)
    - Laplacian kernel function ``:laplace``
    - Linear kernel function ``:linear``

    The choice of ``blocksize`` allows control over the tradeoff between test
    power and computation time (default: :math:`\lfloor \sqrt{n} \rfloor`).


    Implements: :ref:`pvalue<pvalue>`, :ref:`ci<ci>`

    References:

    - Zaremba, W., Gretton, A., Blaschko, M. B-tests: Low Variance Kernel Two-Sample Tests.
      Conference on Neural Information Processing Systems, 2013.

    - Matlab implementation: https://github.com/wojzaremba/btest.

    - Gretton, A., Borgwardt K.M., Rasch, M.J., Schölkopf, B., Smola, A.: A Kernel Two-Sample Test.
      Journal of Machine Learning Research 13:723-773,  2012.
