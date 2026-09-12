\name{level_posterior}
\alias{level_posterior}
\title{Posterior Probabilities and Means for One Wavelet Level}
\description{
  Computes the complete three-component posterior for a vector of observed
  wavelet coefficients at one resolution level.
}
\usage{
level_posterior(
  d,
  pi_j,
  omega_j,
  beta_j,
  sigma,
  likelihood,
  laplace_rate,
  quadrature_n,
  use_exact_wendland
)
}
\arguments{
  \item{d}{A numeric vector of observed wavelet coefficients.  The
  calculation is vectorized, with one posterior calculation for each
  element.}
  \item{pi_j}{The prior probability of the point-mass spike at zero for the
  resolution level.  It must lie in \eqn{[0,1]}.}
  \item{omega_j}{The conditional prior probability of the Wendland slab
  among nonzero coefficients.  It must lie in \eqn{[0,1]}.}
  \item{beta_j}{The positive common support scale for the two continuous
  slabs.}
  \item{sigma}{The positive Gaussian noise standard deviation.  It is used
  for the Gaussian likelihood and for Gaussian quadrature calculations.}
  \item{likelihood}{Either \code{"gaussian"} or \code{"laplace"}.}
  \item{laplace_rate}{The positive rate \eqn{\lambda} used by the Laplace
  working likelihood.  It is ignored when \code{likelihood =
  "gaussian"}, but should still be supplied.}
  \item{quadrature_n}{The number of Gauss-Legendre nodes used for numerical
  integration.}
  \item{use_exact_wendland}{Logical value indicating whether the exact
  finite-sum Wendland formula should be used for Laplace moments.}
}
\details{
  The continuous slab is
  \deqn{
    g_j(\theta)=\omega_j g_W(\theta;\beta_j)
    +(1-\omega_j)g_S(\theta;\beta_j).
  }
  For component \eqn{c\in\{W,S\}}, define
  \deqn{
    M_{r,c,j}(d)=\int_{-\beta_j}^{\beta_j}
    \theta^r g_c(\theta;\beta_j)L(d\mid\theta)\,d\theta,
    \qquad r\in\{0,1\}.
  }
  The posterior component probabilities are
  \deqn{
    p_{0,j}(d)=\frac{\pi_jL(d\mid0)}{D_j(d)},\quad
    p_{W,j}(d)=\frac{(1-\pi_j)\omega_jM_{0,W,j}(d)}{D_j(d)},\quad
    p_{S,j}(d)=\frac{(1-\pi_j)(1-\omega_j)M_{0,S,j}(d)}{D_j(d)},
  }
  where
  \deqn{
    D_j(d)=\pi_jL(d\mid0)+(1-\pi_j)
    \{\omega_jM_{0,W,j}(d)+(1-\omega_j)M_{0,S,j}(d)\}.
  }
  The returned posterior mean is
  \deqn{
    \delta_j(d)=p_{W,j}(d)\mu_{W,j}(d)
    +p_{S,j}(d)\mu_{S,j}(d),
    \qquad
    \mu_{c,j}(d)=M_{1,c,j}(d)/M_{0,c,j}(d).
  }
  The spike contributes zero to this posterior mean.  Computations are
  performed on the log scale where appropriate to reduce numerical
  underflow.  Endpoint cases such as \eqn{\pi_j=1} or
  \eqn{\omega_j=0} are handled explicitly.
}
\value{
  A list with one row or one entry per value of \code{d}:
  \describe{
    \item{\code{log_marginal}}{The log marginal density
    \eqn{\log D_j(d)}.}
    \item{\code{probability}}{A matrix with columns \code{spike},
    \code{wendland}, and \code{semicircle}, containing posterior component
    probabilities.}
    \item{\code{posterior_mean}}{The posterior-mean shrinkage estimate
    \eqn{\delta_j(d)}.}
    \item{\code{component_mean}}{A matrix containing the posterior means
    conditional on the Wendland and semicircle components.}
    \item{\code{component_log_m0}}{A matrix containing the log zeroth
    moments for the two continuous components.}
    \item{\code{quadrature_fallbacks}}{The number of Wendland coefficients
    for which the exact Laplace calculation used numerical fallback.}
  }
}
\examples{
d <- c(-1, 0, 1)
post <- level_posterior(
  d = d,
  pi_j = 0.75,
  omega_j = 0.5,
  beta_j = 2,
  sigma = 1,
  likelihood = "gaussian",
  laplace_rate = 0.5,
  quadrature_n = 24L,
  use_exact_wendland = TRUE
)
post$posterior_mean
post$probability
}
\seealso{
  \code{\link{wswavelet}},
  \code{\link{evaluate_shrinkage_curve}},
  \code{\link{wendland_kernel}},
  \code{\link{semicircle_kernel}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{wavelet}
\keyword{Bayesian}
