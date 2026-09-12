\name{evaluate_shrinkage_curve}
\alias{evaluate_shrinkage_curve}
\title{Evaluate a Fitted Wavelet Shrinkage Curve}
\description{
  Evaluates the fitted posterior-mean shrinkage rule and posterior component
  probabilities over a grid of hypothetical observed coefficients at one
  detail level.
}
\usage{
evaluate_shrinkage_curve(
  fit,
  level_index = 1L,
  d_grid = seq(-4, 4, length.out = 401L)
)
}
\arguments{
  \item{fit}{A fitted object returned by
  \code{\link{wswavelet}}.}
  \item{level_index}{The position of the detail level in \code{fit$detail},
  beginning at 1.  This is the position in the fitted list, not necessarily
  the numerical wavelet level.}
  \item{d_grid}{A numeric vector of hypothetical observed coefficient values
  at which the shrinkage rule is evaluated.}
}
\details{
  For each value in \code{d_grid}, the function recomputes the posterior for
  the selected level using the fitted \eqn{\pi_j}, \eqn{\omega_j},
  \eqn{\beta_j}, and likelihood.  The resulting curve is the empirical-Bayes
  version of the coefficientwise posterior-mean rule
  \deqn{
    \delta_j(d)=p_{W,j}(d)\mu_{W,j}(d)
    +p_{S,j}(d)\mu_{S,j}(d).
  }
  The output can be used to plot attenuation, compare the rule with the
  identity map \eqn{\delta(d)=d}, or inspect how posterior evidence for the
  spike changes with the observed coefficient.
}
\value{
  A data frame with columns \code{d}, \code{estimate},
  \code{posterior_spike}, \code{posterior_wendland}, and
  \code{posterior_semicircle}.
}
\examples{
\donttest{
set.seed(1)
n <- 128
x <- seq(0, 1, length.out = n)
y <- sin(6 * pi * x) + rnorm(n, sd = 0.5)
fit <- wswavelet(
  y, filter.number = 6L, quadrature_n = 24L
)
curve <- evaluate_shrinkage_curve(fit, level_index = 1L)
head(curve)
}
}
\seealso{
  \code{\link{wswavelet}},
  \code{\link{level_posterior}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{wavelet}
\keyword{diagnostics}
