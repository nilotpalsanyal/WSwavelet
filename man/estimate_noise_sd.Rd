\name{estimate_noise_sd}
\alias{estimate_noise_sd}
\title{Robust Wavelet-Based Noise-Scale Estimation}
\description{
  Estimates the observational noise standard deviation from high-frequency
  wavelet coefficients using a robust median absolute deviation rule.
}
\usage{
estimate_noise_sd(
  wavelet_object,
  detail_levels,
  supplied_sd = NULL,
  mad_levels = 1L
)
}
\arguments{
  \item{wavelet_object}{A wavelet object returned by
  \code{\link[wavethresh]{wd}}.}
  \item{detail_levels}{A nonempty numeric or integer vector of detail levels
  available in \code{wavelet_object}.  The finest levels are selected.}
  \item{supplied_sd}{Optional known positive noise standard deviation.  If
  supplied, it is returned directly and the MAD calculation is skipped.}
  \item{mad_levels}{The positive integer number of finest detail levels to
  pool when estimating the scale.  The default is 1.}
}
\details{
  Unless \code{supplied_sd} is provided, the function computes
  \deqn{
    \widehat{\sigma} =
    \frac{\operatorname{median}(|d-\operatorname{median}(d)|)}{0.6745}
  }
  using the finest \code{mad_levels} detail levels.  If the resulting value
  is not positive and finite, the function falls back to pooling up to the
  three finest available detail levels.  A strictly positive machine-scale
  value is returned in the remaining degenerate case.
}
\value{
  A single positive numeric value giving the supplied or estimated noise
  standard deviation.
}
\examples{
\donttest{
set.seed(1)
y <- rnorm(128)
wd <- wavethresh::wd(y, filter.number = 6L, bc = "periodic")
lev <- 0:(wavethresh::nlevelsWT(wd) - 1L)
estimate_noise_sd(wd, lev)
}
}
\seealso{
  \code{\link{wswavelet}},
  \code{\link[wavethresh]{wd}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{wavelet}
\keyword{robust}
