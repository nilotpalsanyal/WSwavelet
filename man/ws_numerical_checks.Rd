\name{ws_numerical_checks}
\alias{ws_numerical_checks}
\title{Numerical Checks for the WS Posterior-Mean Rule}
\description{
  Runs a standard collection of numerical checks for endpoint handling,
  symmetry, support bounds, finiteness, and Gaussian monotonicity of the
  posterior-mean shrinkage rule.
}
\usage{
ws_numerical_checks(quadrature_n = 128L)
}
\arguments{
  \item{quadrature_n}{The number of Gauss-Legendre nodes used in the
  numerical integrations.  Larger values can improve accuracy but increase
  computation time.}
}
\details{
  The routine evaluates Gaussian and Laplace posterior means over a fixed
  coefficient grid for all combinations of spike probability
  \eqn{\pi\in\{0,0.5,1\}} and Wendland weight
  \eqn{\omega\in\{0,0.5,1\}}.  It checks the numerical versions of
  \deqn{\delta(-d)=-\delta(d)}
  and \deqn{|\delta(d)|\leq\beta}.
  For the Gaussian likelihood, \code{min_first_difference} checks the
  expected nondecreasing behavior on the grid.  The endpoint columns verify
  the exact no-spike and all-spike rules.
}
\value{
  A data frame with one row for each likelihood, spike probability, and slab
  weight combination.  Important columns are \code{max_oddness_error},
  \code{min_first_difference}, \code{max_support_excess},
  \code{all_finite}, \code{pi0_exact_slab}, and \code{pi1_exact_zero}.
}
\examples{
\donttest{
checks <- ws_numerical_checks(quadrature_n = 48L)
checks
}
}
\seealso{
  \code{\link{level_posterior}},
  \code{\link{wswavelet}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{diagnostics}
\keyword{wavelet}
