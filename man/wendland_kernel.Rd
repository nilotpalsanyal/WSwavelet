\name{wendland_kernel}
\alias{wendland_kernel}
\title{Standardized Wendland-Type Slab Density}
\description{
  Evaluates the normalized compactly supported Wendland-type slab density on
  the standardized interval \eqn{[-1,1]}.
}
\usage{
wendland_kernel(u)
}
\arguments{
  \item{u}{A numeric vector of standardized coefficient values.}
}
\details{
  The density is
  \deqn{
    K_W(u)=\frac{3}{2}(1-|u|)^4(1+4|u|)
    I(|u|<1).
  }
  It is symmetric, integrates to one, is concentrated near zero, and is
  exactly zero outside the open unit interval.  For a support scale
  \eqn{\beta>0}, the corresponding density for a coefficient \eqn{\theta} is
  \eqn{g_W(\theta;\beta)=\beta^{-1}K_W(\theta/\beta)}.
}
\value{
  A numeric vector of the same length as \code{u}, containing the density
  values.
}
\examples{
u <- seq(-1.2, 1.2, length.out = 7)
wendland_kernel(u)
}
\seealso{
  \code{\link{semicircle_kernel}},
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
