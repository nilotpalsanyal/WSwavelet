\name{wswavelet}
\alias{wswavelet}
\title{Resolution-Adaptive Bayesian Wavelet Denoising}
\description{
  Fits the resolution-adaptive Wendland-Semicircle Bayesian wavelet
  denoiser, shrinks detail coefficients by posterior means, and reconstructs
  the signal with the inverse discrete wavelet transform.
}
\usage{
wswavelet(
  y,
  filter.number = 10L,
  family = "DaubExPhase",
  bc = "periodic",
  j0 = 0L,
  likelihood = c("gaussian", "laplace"),
  supplied_sd = NULL,
  spike_offset = 2,
  spike_gamma = 2.4,
  beta_quantile = 0.99,
  beta_quantile_type = 8L,
  beta_floor = 1e-8,
  omega_model = c("adaptive", "constant", "fixed"),
  fixed_omega = NULL,
  eta_start = c(0, 0),
  eta_bounds = c(-12, 12),
  laplace_rate = NULL,
  quadrature_n = 48L,
  use_exact_wendland = TRUE,
  mad_levels = 1L,
  return_wavelet = FALSE
)
}
\arguments{
  \item{y}{A finite numeric signal vector.  Its length must be a dyadic
  integer, such as 128, 256, or 1024.}
  \item{filter.number}{The wavelet filter order passed to
  \code{\link[wavethresh]{wd}}.  The default is 10.}
  \item{family}{The wavelet family passed to \code{\link[wavethresh]{wd}}.
  Common choices include \code{"DaubExPhase"} and
  \code{"DaubLeAsymm"}.}
  \item{bc}{The boundary convention passed to \code{\link[wavethresh]{wd}}.
  The default is \code{"periodic"}.}
  \item{j0}{The coarsest detail level to shrink.  Detail levels below
  \code{j0} are not shrunk.  The default is 0.}
  \item{likelihood}{The coefficientwise working likelihood.  Use
  \code{"gaussian"} for the primary version or \code{"laplace"} for the
  robustness-oriented alternative.}
  \item{supplied_sd}{Optional known noise standard deviation.  If
  \code{NULL}, the function estimates it from high-frequency wavelet
  coefficients using a robust MAD rule.}
  \item{spike_offset}{The offset \eqn{\ell} in the level-dependent spike
  probability.  The implementation permits values at least 1; the default
  is 2.}
  \item{spike_gamma}{The positive resolution-sparsity exponent
  \eqn{\gamma}.  Larger values make the spike probability increase more
  rapidly with resolution.}
  \item{beta_quantile}{The probability \eqn{\kappa} used in the empirical
  upper-quantile rule for the support scale \eqn{\beta_j}.  It must lie
  strictly between 0 and 1; the default is 0.99.}
  \item{beta_quantile_type}{The integer quantile convention passed to
  \code{\link[stats]{quantile}}.  The default is type 8.}
  \item{beta_floor}{A positive numerical floor used to prevent a zero support
  scale when the noise estimate and coefficients are numerically zero.}
  \item{omega_model}{How the Wendland weight is modeled across resolution
  levels.  \code{"adaptive"} estimates an intercept and a resolution trend,
  \code{"constant"} estimates one common weight, and \code{"fixed"} uses
  \code{fixed_omega}.}
  \item{fixed_omega}{A fixed Wendland weight in \eqn{[0,1]}.  Supplying this
  argument forces \code{omega_model = "fixed"}; 0 gives the
  semicircle-only endpoint and 1 gives the Wendland-only endpoint.}
  \item{eta_start}{Starting values for the empirical-Bayes logistic trend
  parameters \eqn{(\eta_0,\eta_1)}.  For a constant mixture, only the first
  value is used.}
  \item{eta_bounds}{Common finite lower and upper bounds for the optimized
  logistic parameters.  The default is \code{c(-12, 12)}.}
  \item{laplace_rate}{Optional positive rate parameter for the Laplace
  working likelihood.  If \code{NULL}, the function sets
  \eqn{\widehat{\lambda}=1/(2\widehat{\sigma}^2)} once before optimization
  and holds it fixed.  This argument has no effect for the Gaussian
  likelihood.}
  \item{quadrature_n}{Number of Gauss-Legendre nodes used for numerical
  integration.  The default is 48.  Larger values may improve accuracy but
  increase computation time.}
  \item{use_exact_wendland}{Logical value.  For the Laplace likelihood,
  \code{TRUE} uses the finite-sum Wendland moment formula with a numerical
  fallback.  If \code{FALSE}, numerical quadrature is used.  For the
  Gaussian likelihood, numerical quadrature is used.}
  \item{mad_levels}{Number of finest detail levels pooled for the robust
  noise-scale estimate.  The default is 1.}
  \item{return_wavelet}{Logical value.  If \code{TRUE}, also returns the
  original and thresholded \code{wavethresh} wavelet objects.}
}
\details{
  The prior for a detail coefficient at resolution level \eqn{j} is
  \deqn{
    \theta_{j,k} \sim \pi_j\delta_0 +
    (1-\pi_j)\{\omega_j g_W(\theta;\beta_j)
    +(1-\omega_j)g_S(\theta;\beta_j)\}.
  }
  The spike probability is
  \deqn{
    \pi_j = 1-(j-j_0+\ell)^{-\gamma}.
  }
  For the adaptive model, the Wendland weight is
  \deqn{
    \omega_j(\eta)=\operatorname{logit}^{-1}(\eta_0+\eta_1 t_j),
    \qquad
    t_j=(j-j_0)/(J-1-j_0).
  }
  The support scale is selected from the observed detail coefficients using
  \deqn{
    \widehat{\beta}_j =
    \max\{\widehat{\sigma},Q_{\kappa}(|d_{j,k}|)\},
  }
  followed by the positive numerical floor
  \deqn{
    \beta_{\min}=\epsilon_\beta s_\beta,\qquad
    s_\beta=\max\{1,\widehat{\sigma},\max_{j,k}|d_{j,k}|\},\qquad
    \widehat{\beta}_j\leftarrow\max\{\widehat{\beta}_j,\beta_{\min}\},
  }
  where \eqn{\epsilon_\beta} is \code{beta_floor}.

  When \code{omega_model = "constant"}, the fitted Wendland weight is
  \eqn{\operatorname{logit}^{-1}(\eta_0)} at every level.  When
  \code{omega_model = "fixed"}, no weight optimization is performed.

  For Gaussian errors, the coefficientwise likelihood is
  \deqn{
    L_N(d\mid\theta,\sigma^2)=
    (2\pi\sigma^2)^{-1/2}
    \exp\{-(d-\theta)^2/(2\sigma^2)\}.
  }
  For the Laplace working likelihood, the implementation uses
  \deqn{
    L_L(d\mid\theta,\lambda)=(a/2)\exp\{-a|d-\theta|\},
    \qquad a=\sqrt{2\lambda}.
  }
  The empirical-Bayes parameters are fitted by bounded L-BFGS-B
  optimization of the joint log marginal likelihood.  The posterior mean
  replaces each detail coefficient, while the scaling coefficients are
  retained without shrinkage.
}
\value{
  A list with the following principal components:
  \describe{
    \item{\code{estimate}}{The reconstructed denoised signal.}
    \item{\code{sigma_hat}}{The estimated or supplied noise standard
    deviation.}
    \item{\code{laplace_rate_hat}}{The rate used for the Laplace working
    likelihood.}
    \item{\code{likelihood}}{The selected likelihood.}
    \item{\code{omega_model}}{The selected mixture-weight model.}
    \item{\code{eta_hat}}{The fitted logistic-trend parameters, padded to
    length two for the constant and fixed cases.}
    \item{\code{optimization}}{The optimization result or fixed-weight
    record.}
    \item{\code{level_summary}}{A data frame containing level-specific
    support scales, prior probabilities, fitted weights, average posterior
    component probabilities, and numerical diagnostics.}
    \item{\code{detail}}{A list containing the observed and fitted
    coefficients and posterior quantities for each detail level.}
    \item{\code{original}}{The input signal.}
    \item{\code{diagnostics}}{Quadrature, fallback, non-finite-integral,
    quantile-type, and MAD-level information.}
    \item{\code{settings}}{The principal fitting settings.}
  }
  If \code{return_wavelet = TRUE}, the list also contains
  \code{wavelet_object} and \code{thresholded_wavelet}.
}
\examples{
\donttest{
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
}
}
\seealso{
  \code{\link{level_posterior}},
  \code{\link{evaluate_shrinkage_curve}},
  \code{\link{estimate_noise_sd}},
  \code{\link[wavethresh]{wd}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{smooth}
\keyword{wavelet}
\keyword{nonparametric}
