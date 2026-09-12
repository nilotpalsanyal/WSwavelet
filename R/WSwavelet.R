
# -----------------------------------------------------------------------------
# Basic numerical helpers
# -----------------------------------------------------------------------------

.gauss_legendre_cache <- new.env(parent = emptyenv())

# Supply a fallback only when an object is NULL.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Test whether a scalar is an admissible dyadic signal length.
is_power_of_two <- function(n) {
  n >= 2L && is.finite(n) && n == 2^round(log2(n))
}

# Compute row-wise log-sum-exp values without numerical underflow.
row_log_sum_exp <- function(log_matrix) {
  if (!is.matrix(log_matrix)) {
    log_matrix <- matrix(log_matrix, nrow = 1L)
  }
  row_max <- apply(log_matrix, 1L, max)
  row_max + log(rowSums(exp(log_matrix - row_max)))
}

# This avoids an additional quadrature-package dependency and supports repeated
# evaluation of the semicircle and Gaussian-likelihood integrals.
# Generate Gauss--Legendre nodes and weights on a finite interval.
gauss_legendre <- function(n = 48L, lower = -1, upper = 1, tolerance = 1e-14) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least 2.")
  if (!(is.finite(lower) && is.finite(upper) && lower < upper)) {
    stop("The quadrature interval must satisfy lower < upper.")
  }
  cache_key <- paste(n, lower, upper, tolerance, sep = "|")
  if (exists(cache_key, .gauss_legendre_cache, inherits = FALSE))
    return(get(cache_key, .gauss_legendre_cache, inherits = FALSE))

  nodes <- numeric(n)
  weights <- numeric(n)
  half_n <- ceiling(n / 2)
  midpoint <- (lower + upper) / 2
  half_length <- (upper - lower) / 2

  for (i in seq_len(half_n)) {
    z <- cos(pi * (i - 0.25) / (n + 0.5))
    for (iteration in seq_len(100L)) {
      p0 <- 1
      p1 <- z
      if (n >= 2L) {
        for (k in 2:n) {
          p2 <- p1
          p1 <- ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
          p0 <- p2
        }
      }
      derivative <- n * (z * p1 - p0) / (z^2 - 1)
      updated_z <- z - p1 / derivative
      if (abs(updated_z - z) < tolerance) {
        z <- updated_z
        break
      }
      z <- updated_z
    }

    # Recompute the derivative at the converged root.
    p0 <- 1
    p1 <- z
    if (n >= 2L) {
      for (k in 2:n) {
        p2 <- p1
        p1 <- ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
        p0 <- p2
      }
    }
    derivative <- n * (z * p1 - p0) / (z^2 - 1)
    node_weight <- 2 / ((1 - z^2) * derivative^2)

    left <- i
    right <- n + 1L - i
    nodes[left] <- midpoint - half_length * z
    nodes[right] <- midpoint + half_length * z
    weights[left] <- half_length * node_weight
    weights[right] <- half_length * node_weight
  }

  result <- list(x = nodes, w = weights)
  assign(cache_key, result, .gauss_legendre_cache)
  result
}

# Evaluate the normalized standardized Wendland slab density.
wendland_kernel <- function(u) {
  abs_u <- abs(u)
  ifelse(
    abs_u < 1,
    1.5 * (1 - abs_u)^4 * (1 + 4 * abs_u),
    0
  )
}

# Evaluate the normalized standardized semicircle slab density.
semicircle_kernel <- function(u) {
  abs_u <- abs(u)
  ifelse(abs_u < 1, 2 / pi * sqrt(pmax(0, 1 - u^2)), 0)
}

# Evaluate a coefficient-by-node log-likelihood matrix.
likelihood_log_density <- function(d, theta, sigma, likelihood, laplace_rate) {
  residual <- outer(d, theta, "-")
  if (likelihood == "gaussian") {
    -0.5 * log(2 * pi * sigma^2) - residual^2 / (2 * sigma^2)
  } else {
    a <- sqrt(2 * laplace_rate)
    matrix(log(a / 2), nrow = length(d), ncol = length(theta)) - a * abs(residual)
  }
}

# Evaluate log likelihoods when each coefficient has distinct nodes.
likelihood_log_density_rowwise <- function(d, theta_matrix, sigma, likelihood,
                                            laplace_rate) {
  residual <- sweep(theta_matrix, 1L, d, "-")
  if (likelihood == "gaussian") {
    -0.5 * log(2 * pi * sigma^2) - residual^2 / (2 * sigma^2)
  } else {
    a <- sqrt(2 * laplace_rate)
    matrix(log(a / 2), nrow = nrow(theta_matrix), ncol = ncol(theta_matrix)) -
      a * abs(residual)
  }
}

# Convert quadrature weights into log marginals and conditional means.
weighted_moments <- function(theta_matrix, prior_weights, log_likelihood) {
  if (is.vector(theta_matrix)) {
    theta_matrix <- matrix(theta_matrix, nrow = 1L)
  }
  if (is.vector(prior_weights)) {
    prior_weights <- matrix(
      rep(prior_weights, each = nrow(theta_matrix)),
      nrow = nrow(theta_matrix)
    )
  }
  log_integrand <- log(prior_weights) + log_likelihood
  log_m0 <- row_log_sum_exp(log_integrand)
  relative_weights <- exp(log_integrand - log_m0)
  posterior_mean <- rowSums(relative_weights * theta_matrix)
  list(log_m0 = log_m0, posterior_mean = posterior_mean)
}

# -----------------------------------------------------------------------------
# Exact Wendland moments under the Laplace working likelihood
# -----------------------------------------------------------------------------

# Evaluate an elementary polynomial-exponential definite integral.
elementary_exponential_integral <- function(m, b, lower, upper) {
  if (abs(b) < 1e-10) {
    return((upper^(m + 1) - lower^(m + 1)) / (m + 1))
  }

  # Evaluate the closed-form antiderivative at one endpoint.
  antiderivative <- function(t) {
    terms <- vapply(
      0:m,
      function(r) {
        (-1)^r * factorial(m) / factorial(m - r) *
          t^(m - r) / b^(r + 1)
      },
      numeric(1)
    )
    exp(b * t) * sum(terms)
  }

  value <- antiderivative(upper) - antiderivative(lower)
  if (is.finite(value)) return(value)

  # This fallback is used only when the finite-sum expression loses floating-
  # point precision for unusually large |b| or |t|.
  stats::integrate(
    function(t) t^m * exp(b * t),
    lower = lower,
    upper = upper,
    rel.tol = 1e-10,
    abs.tol = 1e-12
  )$value
}

# Compute exact zeroth/first Wendland moments under a Laplace likelihood.
wendland_laplace_moments_exact <- function(d, beta, laplace_rate) {
  if (!(is.finite(d) && is.finite(beta) && beta > 0 && laplace_rate > 0)) {
    return(c(m0 = NA_real_, m1 = NA_real_))
  }

  a <- sqrt(2 * laplace_rate)
  polynomial_coefficients <- c(1, 0, -10, 20, -15, 4)
  breakpoints <- c(-beta, 0, beta)
  if (d > -beta && d < beta) breakpoints <- c(breakpoints, d)
  breakpoints <- sort(unique(breakpoints))

  moment_values <- c(0, 0)
  for (i in seq_len(length(breakpoints) - 1L)) {
    lower <- breakpoints[i]
    upper <- breakpoints[i + 1L]
    midpoint <- (lower + upper) / 2
    sign_theta <- sign(midpoint)
    sign_residual <- sign(d - midpoint)

    for (r in 0:1) {
      interval_sum <- 0
      for (m in 0:5) {
        interval_sum <- interval_sum +
          polynomial_coefficients[m + 1L] * beta^(-m) *
          sign_theta^m * elementary_exponential_integral(
            m = r + m,
            b = a * sign_residual,
            lower = lower,
            upper = upper
          )
      }
      moment_values[r + 1L] <- moment_values[r + 1L] +
        exp(-a * sign_residual * d) * interval_sum
    }
  }

  moment_values <- 3 * a / (4 * beta) * moment_values
  names(moment_values) <- c("m0", "m1")
  moment_values
}

# -----------------------------------------------------------------------------
# Component moment calculations
# -----------------------------------------------------------------------------

# Compute Wendland component moments by Gauss--Legendre quadrature.
wendland_moments_quadrature <- function(d, beta, sigma, likelihood,
                                        laplace_rate, quadrature_n) {
  rule <- gauss_legendre(quadrature_n, -1, 1)
  theta <- beta * rule$x
  prior_weights <- rule$w * wendland_kernel(rule$x)
  log_likelihood <- likelihood_log_density(
    d = d,
    theta = theta,
    sigma = sigma,
    likelihood = likelihood,
    laplace_rate = laplace_rate
  )
  weighted_moments(
    theta_matrix = matrix(rep(theta, each = length(d)), nrow = length(d)),
    prior_weights = prior_weights,
    log_likelihood = log_likelihood
  )
}

# Use exact Laplace Wendland moments with a quadrature fallback.
wendland_moments <- function(d, beta, sigma, likelihood, laplace_rate,
                             quadrature_n, use_exact_wendland) {
  if (likelihood != "laplace" || !use_exact_wendland) {
    result <- wendland_moments_quadrature(
      d, beta, sigma, likelihood, laplace_rate, quadrature_n
    )
    result$fallback_n <- 0L
    return(result)
  }

  exact_values <- t(vapply(
    d,
    function(value) wendland_laplace_moments_exact(value, beta, laplace_rate),
    numeric(2)
  ))
  colnames(exact_values) <- c("m0", "m1")
  valid <- is.finite(exact_values[, "m0"]) &
    is.finite(exact_values[, "m1"]) & exact_values[, "m0"] > 0

  result <- list(
    log_m0 = rep(NA_real_, length(d)),
    posterior_mean = rep(NA_real_, length(d)),
    fallback_n = sum(!valid)
  )
  result$log_m0[valid] <- log(exact_values[valid, "m0"])
  result$posterior_mean[valid] <-
    exact_values[valid, "m1"] / exact_values[valid, "m0"]

  # The quadrature fallback protects against cancellation in the exact
  # expression at extreme parameter values while preserving the exact formula
  # for the ordinary parameter range.
  if (any(!valid)) {
    fallback <- wendland_moments_quadrature(
      d[!valid], beta, sigma, likelihood, laplace_rate, quadrature_n
    )
    result$log_m0[!valid] <- fallback$log_m0
    result$posterior_mean[!valid] <- fallback$posterior_mean
  }
  result
}

# Compute semicircle component moments, splitting Laplace kinks safely.
semicircle_moments <- function(d, beta, sigma, likelihood, laplace_rate,
                               quadrature_n) {
  rule <- gauss_legendre(quadrature_n, -1, 1)
  inside <- abs(d) < beta
  log_m0 <- numeric(length(d))
  posterior_mean <- numeric(length(d))

  # Outside the support, the absolute-value kink is absent and one full
  # Gauss--Legendre rule on [0, pi] is sufficient.
  if (any(!inside)) {
    d_out <- d[!inside]
    t <- (rule$x + 1) * pi / 2
    theta <- beta * cos(t)
    prior_weights <- 2 / pi * (rule$w * pi / 2) * sin(t)^2
    log_likelihood <- likelihood_log_density(
      d = d_out,
      theta = theta,
      sigma = sigma,
      likelihood = likelihood,
      laplace_rate = laplace_rate
    )
    result <- weighted_moments(
      theta_matrix = matrix(rep(theta, each = length(d_out)), nrow = length(d_out)),
      prior_weights = prior_weights,
      log_likelihood = log_likelihood
    )
    log_m0[!inside] <- result$log_m0
    posterior_mean[!inside] <- result$posterior_mean
  }

  # For |d| < beta, split at t_d = acos(d / beta), as recommended in the
  # manuscript.  The row-specific quadrature nodes make this vectorized over
  # all coefficients at a given resolution level.
  if (any(inside)) {
    d_in <- d[inside]
    t_split <- acos(d_in / beta)
    first_t <- outer(t_split, rule$x, function(s, x) (x + 1) * s / 2)
    second_t <- outer(
      t_split,
      rule$x,
      function(s, x) s + (x + 1) * (pi - s) / 2
    )
    first_weights <- sweep(
      matrix(rule$w, nrow = length(d_in), ncol = quadrature_n, byrow = TRUE),
      1L,
      t_split / 2,
      "*"
    )
    second_weights <- sweep(
      matrix(rule$w, nrow = length(d_in), ncol = quadrature_n, byrow = TRUE),
      1L,
      (pi - t_split) / 2,
      "*"
    )
    first_weights <- 2 / pi * first_weights * sin(first_t)^2
    second_weights <- 2 / pi * second_weights * sin(second_t)^2
    theta_matrix <- cbind(beta * cos(first_t), beta * cos(second_t))
    prior_weights <- cbind(first_weights, second_weights)
    log_likelihood <- likelihood_log_density_rowwise(
      d = d_in,
      theta_matrix = theta_matrix,
      sigma = sigma,
      likelihood = likelihood,
      laplace_rate = laplace_rate
    )
    result <- weighted_moments(
      theta_matrix = theta_matrix,
      prior_weights = prior_weights,
      log_likelihood = log_likelihood
    )
    log_m0[inside] <- result$log_m0
    posterior_mean[inside] <- result$posterior_mean
  }

  list(log_m0 = log_m0, posterior_mean = posterior_mean, fallback_n = 0L)
}

# Compute both slab-component moments for one resolution level.
component_moments <- function(d, beta, sigma, likelihood, laplace_rate,
                              quadrature_n, use_exact_wendland) {
  list(
    wendland = wendland_moments(
      d = d,
      beta = beta,
      sigma = sigma,
      likelihood = likelihood,
      laplace_rate = laplace_rate,
      quadrature_n = quadrature_n,
      use_exact_wendland = use_exact_wendland
    ),
    semicircle = semicircle_moments(
      d = d,
      beta = beta,
      sigma = sigma,
      likelihood = likelihood,
      laplace_rate = laplace_rate,
      quadrature_n = quadrature_n
    )
  )
}

# Combine cached component moments into endpoint-safe posterior output.
posterior_from_moments <- function(log_spike, moments, pi_j, omega_j) {
  if (!is.finite(pi_j) || pi_j < 0 || pi_j > 1 ||
      !is.finite(omega_j) || omega_j < 0 || omega_j > 1) {
    stop("pi_j and omega_j must lie in [0, 1].")
  }
  log_terms <- cbind(
    if (pi_j == 0) -Inf else log(pi_j) + log_spike,
    if (pi_j == 1 || omega_j == 0) -Inf else
      log1p(-pi_j) + log(omega_j) + moments$wendland$log_m0,
    if (pi_j == 1 || omega_j == 1) -Inf else
      log1p(-pi_j) + log1p(-omega_j) + moments$semicircle$log_m0
  )
  log_denominator <- row_log_sum_exp(log_terms)
  probability <- exp(log_terms - log_denominator)
  colnames(probability) <- c("spike", "wendland", "semicircle")
  posterior_mean <- probability[, "wendland"] * moments$wendland$posterior_mean +
    probability[, "semicircle"] * moments$semicircle$posterior_mean
  list(log_marginal = log_denominator, probability = probability,
       posterior_mean = posterior_mean,
       component_mean = cbind(wendland = moments$wendland$posterior_mean,
                              semicircle = moments$semicircle$posterior_mean),
       component_log_m0 = cbind(wendland = moments$wendland$log_m0,
                               semicircle = moments$semicircle$log_m0),
       quadrature_fallbacks = moments$wendland$fallback_n %||% 0L)
}

# Compute a complete posterior for one level from raw coefficients.
level_posterior <- function(d, pi_j, omega_j, beta_j, sigma, likelihood,
                            laplace_rate, quadrature_n, use_exact_wendland) {
  moments <- component_moments(
    d = d,
    beta = beta_j,
    sigma = sigma,
    likelihood = likelihood,
    laplace_rate = laplace_rate,
    quadrature_n = quadrature_n,
    use_exact_wendland = use_exact_wendland
  )

  log_spike <- if (likelihood == "gaussian") {
    stats::dnorm(d, mean = 0, sd = sigma, log = TRUE)
  } else {
    a <- sqrt(2 * laplace_rate)
    log(a / 2) - a * abs(d)
  }

  posterior_from_moments(log_spike, moments, pi_j, omega_j)
}

# Estimate noise by finest-level or pooled finest-level MAD.
estimate_noise_sd <- function(wavelet_object, detail_levels, supplied_sd = NULL,
                              mad_levels = 1L) {
  if (!is.null(supplied_sd)) {
    if (!is.numeric(supplied_sd) || length(supplied_sd) != 1L ||
        !is.finite(supplied_sd) || supplied_sd <= 0) {
      stop("supplied_sd must be one positive finite number.")
    }
    return(as.numeric(supplied_sd))
  }

  if (length(mad_levels) != 1L || !is.finite(mad_levels) || mad_levels < 1L ||
      mad_levels != as.integer(mad_levels)) {
    stop("mad_levels must be one positive integer.")
  }
  mad_levels <- as.integer(mad_levels)
  used_levels <- tail(detail_levels, min(mad_levels, length(detail_levels)))
  coefficients <- unlist(lapply(used_levels, function(level)
    wavethresh::accessD(wavelet_object, level = level)), use.names = FALSE)
  sigma_hat <- stats::median(abs(coefficients - stats::median(coefficients))) / 0.6745

  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    pooled_levels <- tail(detail_levels, min(3L, length(detail_levels)))
    pooled <- unlist(lapply(
      pooled_levels,
      function(level) wavethresh::accessD(wavelet_object, level = level)
    ))
    sigma_hat <- stats::median(abs(pooled - stats::median(pooled))) / 0.6745
  }

  # A strictly positive fallback is needed for a constant or numerically
  # noiseless input; ordinary noisy data use the MAD estimate above.
  max(as.numeric(sigma_hat), sqrt(.Machine$double.eps))
}

# -----------------------------------------------------------------------------
# Main denoising function
# -----------------------------------------------------------------------------

# Fit the complete empirical-Bayes WS posterior-mean estimator.
wswavelet <- function(
    y, filter.number = 10L, family = "DaubExPhase", bc = "periodic", j0 = 0L,
    likelihood = c("gaussian", "laplace"), supplied_sd = NULL,
    spike_offset = 2, spike_gamma = 2.4, beta_quantile = .99,
    beta_quantile_type = 8L, beta_floor = 1e-8,
    omega_model = c("adaptive", "constant", "fixed"), fixed_omega = NULL,
    eta_start = c(0, 0), eta_bounds = c(-12, 12), laplace_rate = NULL,
    quadrature_n = 48L, use_exact_wendland = TRUE, mad_levels = 1L,
    return_wavelet = FALSE
) {
  likelihood <- match.arg(likelihood)
  omega_model <- if (!is.null(fixed_omega)) "fixed" else match.arg(omega_model)
  y <- as.numeric(y); n <- length(y)
  if (!is_power_of_two(n) || any(!is.finite(y)))
    stop("y must be finite and have dyadic length at least 2.")
  if (length(j0) != 1L || j0 < 0 || j0 != as.integer(j0))
    stop("j0 must be a nonnegative integer.")
  if (!is.finite(spike_offset) || spike_offset < 1 ||
      !is.finite(spike_gamma) || spike_gamma <= 0)
    stop("spike_offset must be at least 1 and spike_gamma must be positive.")
  if (!is.finite(beta_quantile) || beta_quantile <= 0 || beta_quantile >= 1)
    stop("beta_quantile must lie strictly between 0 and 1.")
  if (length(beta_quantile_type) != 1L || beta_quantile_type < 1 ||
      beta_quantile_type > 9 || beta_quantile_type != as.integer(beta_quantile_type))
    stop("beta_quantile_type must be an integer from 1 through 9.")
  if (!is.finite(beta_floor) || beta_floor <= 0)
    stop("beta_floor must be positive.")
  if (omega_model == "fixed" && (is.null(fixed_omega) || length(fixed_omega) != 1L ||
      !is.finite(fixed_omega) || fixed_omega < 0 || fixed_omega > 1))
    stop("fixed_omega must be one number in [0, 1] for omega_model='fixed'.")
  if (length(eta_start) != 2L || any(!is.finite(eta_start)) ||
      length(eta_bounds) != 2L || any(!is.finite(eta_bounds)) ||
      eta_bounds[1L] >= eta_bounds[2L])
    stop("eta_start and eta_bounds must contain valid finite pairs.")
  if (length(quadrature_n) != 1L || quadrature_n < 8L ||
      quadrature_n != as.integer(quadrature_n))
    stop("quadrature_n must be an integer of at least 8.")

  wavelet_object <- wavethresh::wd(y, filter.number = filter.number,
    family = family, type = "wavelet", bc = bc, verbose = FALSE)
  number_of_levels <- wavethresh::nlevelsWT(wavelet_object)
  if (j0 > number_of_levels - 1L) stop("j0 leaves no detail levels to shrink.")
  detail_levels <- seq.int(j0, number_of_levels - 1L)
  sigma_hat <- estimate_noise_sd(wavelet_object, detail_levels, supplied_sd,
                                 mad_levels = mad_levels)
  if (is.null(laplace_rate)) laplace_rate_hat <- 1 / (2 * sigma_hat^2) else {
    if (length(laplace_rate) != 1L || !is.finite(laplace_rate) || laplace_rate <= 0)
      stop("laplace_rate must be NULL or one positive finite number.")
    laplace_rate_hat <- as.numeric(laplace_rate)
  }

  detail_coefficients <- lapply(detail_levels, function(level)
    as.numeric(wavethresh::accessD(wavelet_object, level = level)))
  names(detail_coefficients) <- as.character(detail_levels)
  scale_reference <- max(1, sigma_hat, abs(unlist(detail_coefficients, use.names = FALSE)))
  support_scales <- vapply(detail_coefficients, function(d) max(sigma_hat,
    stats::quantile(abs(d), beta_quantile, names = FALSE, type = beta_quantile_type),
    beta_floor * scale_reference), numeric(1))
  spike_probabilities <- pmin(1, pmax(0,
    1 - (detail_levels - j0 + spike_offset)^(-spike_gamma)))
  level_coordinate <- if (length(detail_levels) == 1L) 0 else
    (detail_levels - j0) / (max(detail_levels) - j0)

  # Precompute all eta-invariant component integrals once per level.
  moment_cache <- lapply(seq_along(detail_levels), function(i) component_moments(
    detail_coefficients[[i]], support_scales[i], sigma_hat, likelihood,
    laplace_rate_hat, quadrature_n, use_exact_wendland))
  log_spike <- lapply(detail_coefficients, function(d) if (likelihood == "gaussian")
    stats::dnorm(d, sd = sigma_hat, log = TRUE) else {
      a <- sqrt(2 * laplace_rate_hat); log(a / 2) - a * abs(d)
    })

  # Map EB parameters to adaptive or constant level weights.
  omega_from_eta <- function(eta) if (omega_model == "adaptive")
    stats::plogis(eta[1L] + eta[2L] * level_coordinate) else
      rep(stats::plogis(eta[1L]), length(detail_levels))

  # Evaluate the cached joint log marginal for EB optimization.
  joint_log_marginal <- function(eta) {
    omega <- omega_from_eta(eta)
    sum(vapply(seq_along(detail_levels), function(i) sum(
      posterior_from_moments(log_spike[[i]], moment_cache[[i]],
                             spike_probabilities[i], omega[i])$log_marginal), numeric(1)))
  }

  if (omega_model == "fixed") {
    eta_hat <- c(if (fixed_omega == 0) -Inf else if (fixed_omega == 1) Inf else
                   stats::qlogis(fixed_omega), 0)
    omega_hat <- rep(fixed_omega, length(detail_levels))
    fixed_log_marginal <- sum(vapply(seq_along(detail_levels), function(i) sum(
      posterior_from_moments(log_spike[[i]], moment_cache[[i]],
        spike_probabilities[i], omega_hat[i])$log_marginal), numeric(1)))
    optimization <- list(par = eta_hat, value = -fixed_log_marginal,
                         convergence = 0L, message = "fixed mixture weight")
  } else {
    parameters <- if (omega_model == "adaptive") 2L else 1L
    start <- pmin(pmax(eta_start[seq_len(parameters)], eta_bounds[1L]), eta_bounds[2L])
    optimization <- stats::optim(start, function(eta) -joint_log_marginal(eta),
      method = "L-BFGS-B", lower = rep(eta_bounds[1L], parameters),
      upper = rep(eta_bounds[2L], parameters), control = list(maxit = 200L))
    eta_hat <- if (parameters == 1L) c(optimization$par, 0) else optimization$par
    omega_hat <- omega_from_eta(optimization$par)
  }

  fitted_details <- vector("list", length(detail_levels)); thresholded <- wavelet_object
  for (i in seq_along(detail_levels)) {
    posterior <- posterior_from_moments(log_spike[[i]], moment_cache[[i]],
      spike_probabilities[i], omega_hat[i])
    fitted_details[[i]] <- list(level = detail_levels[i], observed = detail_coefficients[[i]],
      estimate = posterior$posterior_mean, posterior_probability = posterior$probability,
      component_mean = posterior$component_mean,
      component_log_m0 = posterior$component_log_m0,
      beta = support_scales[i], spike_probability = spike_probabilities[i],
      omega = omega_hat[i], quadrature_fallbacks = posterior$quadrature_fallbacks)
    thresholded <- wavethresh::putD(thresholded, level = detail_levels[i],
                                    v = posterior$posterior_mean)
  }
  level_summary <- do.call(rbind, lapply(seq_along(fitted_details), function(i) {
    z <- fitted_details[[i]]
    data.frame(level = z$level, coordinate = level_coordinate[i],
      coefficients = length(z$observed), spike_probability = z$spike_probability,
      beta = z$beta, omega = z$omega,
      mean_p0 = mean(z$posterior_probability[, "spike"]),
      mean_pW = mean(z$posterior_probability[, "wendland"]),
      mean_pS = mean(z$posterior_probability[, "semicircle"]),
      mean_log_m0_W = mean(z$component_log_m0[, "wendland"]),
      mean_log_m0_S = mean(z$component_log_m0[, "semicircle"]),
      quadrature_fallbacks = z$quadrature_fallbacks, stringsAsFactors = FALSE)
  }))
  estimate <- as.numeric(wavethresh::wr(thresholded))
  nonfinite <- sum(!is.finite(unlist(lapply(fitted_details, function(z)
    c(z$estimate, z$component_mean, z$component_log_m0)))))
  result <- list(estimate = estimate, sigma_hat = sigma_hat,
    laplace_rate_hat = laplace_rate_hat, likelihood = likelihood,
    omega_model = omega_model, eta_hat = eta_hat, optimization = optimization,
    level_summary = level_summary, detail = fitted_details, original = y,
    diagnostics = list(quadrature_n = quadrature_n,
      quadrature_fallbacks = sum(level_summary$quadrature_fallbacks),
      nonfinite_integrals = nonfinite, beta_quantile_type = beta_quantile_type,
      mad_levels = mad_levels),
    settings = list(filter.number = filter.number, family = family, bc = bc, j0 = j0,
      spike_offset = spike_offset, spike_gamma = spike_gamma,
      beta_quantile = beta_quantile, beta_quantile_type = beta_quantile_type,
      beta_floor = beta_floor, quadrature_n = quadrature_n,
      use_exact_wendland = use_exact_wendland, mad_levels = mad_levels))
  if (return_wavelet) {
    result$wavelet_object <- wavelet_object; result$thresholded_wavelet <- thresholded
  }
  result
}

# -----------------------------------------------------------------------------
# Optional inspection helpers
# -----------------------------------------------------------------------------

# Evaluate the fitted shrinkage curve and component probabilities.
evaluate_shrinkage_curve <- function(fit, level_index = 1L,
                                     d_grid = seq(-4, 4, length.out = 401L)) {
  if (!inherits(fit, "list") || is.null(fit$detail)) {
    stop("fit must be the result returned by wswavelet().")
  }
  if (level_index < 1L || level_index > length(fit$detail)) {
    stop("level_index is outside the fitted detail levels.")
  }
  detail <- fit$detail[[level_index]]
  posterior <- level_posterior(
    d = d_grid,
    pi_j = detail$spike_probability,
    omega_j = detail$omega,
    beta_j = detail$beta,
    sigma = fit$sigma_hat,
    likelihood = fit$likelihood,
    laplace_rate = fit$laplace_rate_hat,
    quadrature_n = fit$settings$quadrature_n,
    use_exact_wendland = fit$settings$use_exact_wendland
  )
  data.frame(
    d = d_grid,
    estimate = posterior$posterior_mean,
    posterior_spike = posterior$probability[, "spike"],
    posterior_wendland = posterior$probability[, "wendland"],
    posterior_semicircle = posterior$probability[, "semicircle"],
    stringsAsFactors = FALSE
  )
}

# Check endpoint rules, symmetry, support bounds, and monotonicity.
ws_numerical_checks <- function(quadrature_n = 128L) {
  d <- seq(-8, 8, length.out = 801L); sigma <- 1; beta <- 3
  cases <- expand.grid(likelihood = c("gaussian", "laplace"),
    pi = c(0, .5, 1), omega = c(0, .5, 1), stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
    likelihood <- cases$likelihood[i]; rate <- 1 / (2 * sigma^2)
    posterior <- level_posterior(d, cases$pi[i], cases$omega[i], beta, sigma,
      likelihood, rate, quadrature_n, use_exact_wendland = TRUE)
    z <- posterior$posterior_mean
    data.frame(likelihood = likelihood, pi = cases$pi[i], omega = cases$omega[i],
      max_oddness_error = max(abs(z + rev(z))),
      min_first_difference = min(diff(z)),
      max_support_excess = max(pmax(abs(z) - beta, 0)),
      all_finite = all(is.finite(z)),
      pi0_exact_slab = if (cases$pi[i] == 0)
        max(posterior$probability[, "spike"]) == 0 else NA,
      pi1_exact_zero = if (cases$pi[i] == 1) max(abs(z)) == 0 else NA)
  }))
}

# Example usage (kept commented out so sourcing this file does not run a fit):
#
# y <- your_signal_vector
# fit <- wswavelet(
#   y = y,
#   likelihood = "laplace",
#   filter.number = 10L,
#   family = "DaubExPhase",
#   bc = "periodic",
#   beta_quantile = 0.99,
#   quadrature_n = 48L
# )
# denoised_signal <- fit$estimate
# fit$level_summary
# curve <- evaluate_shrinkage_curve(fit, level_index = 1L)
#
# See Simulation_Analysis_snowfall.R for all parallel study drivers.
