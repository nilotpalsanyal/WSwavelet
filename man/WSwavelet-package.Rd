\name{WSwavelet-package}
\alias{WSwavelet}
\alias{WSwavelet-package}
\docType{package}
\title{WSwavelet: Bayesian Wavelet Denoising with Wendland-Semicircle Slab Mixture}
\description{
  Bayesian wavelet denoising with a resolution-adaptive spike-and-slab prior.
  The continuous slab is a mixture of a compactly supported Wendland-type
  density and the semicircle density.  The package provides Gaussian and
  Laplace working-likelihood versions, empirical-Bayes fitting of
  resolution-specific slab weights, robust noise-scale estimation,
  posterior-mean coefficient shrinkage, inverse-DWT reconstruction, and
  diagnostic tools.
}
\details{
  The primary user-facing function is
  \code{\link{wswavelet}}.  It requires a finite numeric
  vector whose length is a dyadic integer, such as 128 or 1024.  The function
  applies an orthogonal discrete wavelet transform using
  \pkg{wavethresh}, retains the scaling coefficients, shrinks detail
  coefficients, and reconstructs the denoised signal.

  The prior for a detail coefficient at level \eqn{j} is
  \deqn{
    \theta_{j,k} \sim \pi_j\delta_0 +
    (1-\pi_j)\{\omega_j g_W(\theta;\beta_j)
    +(1-\omega_j)g_S(\theta;\beta_j)\}.
  }
  Here \eqn{\pi_j} is the spike probability, \eqn{\omega_j} is the
  conditional Wendland weight among nonzero coefficients, and \eqn{\beta_j}
  is the common support scale.

  See \code{\link{wswavelet}} for the empirical-Bayes
  fitting details and \code{\link{level_posterior}} for the coefficient-level
  posterior calculation.
}
\seealso{
  \code{\link{wswavelet}},
  \code{\link{level_posterior}},
  \code{\link{evaluate_shrinkage_curve}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{
  Nilotpal Sanyal
}
\keyword{package}
