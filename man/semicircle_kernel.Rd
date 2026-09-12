\name{semicircle_kernel}
\alias{semicircle_kernel}
\title{Standardized Semicircle Slab Density}
\description{
  Evaluates the normalized semicircle slab density on the standardized
  interval \eqn{[-1,1]}.
}
\usage{
semicircle_kernel(u)
}
\arguments{
  \item{u}{A numeric vector of standardized coefficient values.}
}
\details{
  The density is
  \deqn{
    K_S(u)=\frac{2}{\pi}\sqrt{1-u^2}I(|u|<1).
  }
  It is symmetric, integrates to one, and is more dispersed than the
  Wendland component used in the package.  For a support scale
  \eqn{\beta>0}, the corresponding density for a coefficient \eqn{\theta} is
  \eqn{g_S(\theta;\beta)=\beta^{-1}K_S(\theta/\beta)}.
}
\value{
  A numeric vector of the same length as \code{u}, containing the density
  values.
}
\examples{
u <- seq(-1.2, 1.2, length.out = 7)
semicircle_kernel(u)
}
\seealso{
  \code{\link{wendland_kernel}},
  \code{\link{level_posterior}}
}
\references{
  Sanyal, N. (2026). Resolution-Adaptive Compact-Support Priors for Bayesian
  Wavelet Denoising: A Wendland-Semicircle Slab Mixture for Low-SNR Signal
  Recovery. \emph{Axioms}, \emph{15}(9), 678.
  <DOI:10.3390/axioms15090678>
}
\author{Nilotpal Sanyal}
\keyword{distribution}
\keyword{wavelet}
