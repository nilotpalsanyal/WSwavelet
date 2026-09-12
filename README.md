
<img src="man/figures/logo.png" align="right" style="float:right; height:122px;" />

# WSwavelet: Bayesian Wavelet Denoising with Wendland-Semicircle Slab Mixture

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/WSwavelet)](https://CRAN.R-project.org/package=WSwavelet)
[![R-CMD-check](https://github.com/nilotpalsanyal/WSwavelet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nilotpalsanyal/WSwavelet/actions/workflows/R-CMD-check.yaml)
[![CodeFactor](https://www.codefactor.io/repository/github/nilotpalsanyal/WSwavelet/badge)](https://www.codefactor.io/repository/github/nilotpalsanyal/WSwavelet)
[![](http://cranlogs.r-pkg.org/badges/grand-total/WSwavelet)](https://cran.r-project.org/package=WSwavelet)

<!-- badges: end -->

WSwavelet implements resolution-adaptive Bayesian wavelet denoising with a
spike-and-slab prior. The continuous slab is a mixture of a compactly
supported Wendland-type density and the semicircle density. The package
supports Gaussian errors as the primary likelihood and a Laplace
working-likelihood sensitivity option.

The main function is wswavelet(). It performs an
orthogonal discrete wavelet transform, estimates the noise scale robustly,
constructs resolution-specific support scales and spike probabilities, fits
the Wendland mixture weight by empirical Bayes, shrinks detail coefficients by
posterior means, and reconstructs the signal.

## Installation

Install the dependency first:

    install.packages("wavethresh")

Then install the package archive:

    install.packages("WSwavelet_0.1.0.tar.gz", repos = NULL, type = "source")

## Minimal example

    library(WSwavelet)
    set.seed(1)
    n <- 128
    x <- seq(0, 1, length.out = n)
    y <- sin(6 * pi * x) + rnorm(n, sd = 0.5)

    fit <- wswavelet(
      y,
      likelihood = "gaussian",
      filter.number = 6L,
      quadrature_n = 24L
    )
    fit$estimate
    fit$level_summary

The input length must be a finite dyadic integer, such as 64, 128, or 256.
The scaling coefficients are retained, and the detail coefficients are
replaced by their posterior-mean estimates.

## Citation

Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
Recovery. *Axioms*, *15*(9), 678.
<DOI:10.3390/axioms15090678>
