## MiSo with a variational Dirichlet distribution for the mixture weights.
##
## This file is deliberately separate from code/miso.R while the behavior of
## the fully Bayesian mixture-weight update is being investigated.

if (!exists("miso_fixed_gamma")) {
  source("code/miso.R")
}

dirichlet_expected_log <- function(parameter) {
  digamma(parameter) - digamma(sum(parameter))
}

dirichlet_kl <- function(parameter, prior_parameter) {
  if (length(parameter) != length(prior_parameter)) {
    stop("parameter and prior_parameter must have the same length.")
  }
  if (any(!is.finite(parameter)) || any(parameter <= 0) ||
      any(!is.finite(prior_parameter)) || any(prior_parameter <= 0)) {
    stop("Dirichlet parameters must be finite and strictly positive.")
  }

  lgamma(sum(parameter)) - sum(lgamma(parameter)) -
    lgamma(sum(prior_parameter)) + sum(lgamma(prior_parameter)) +
    sum(
      (parameter - prior_parameter) * dirichlet_expected_log(parameter)
    )
}

update_dirichlet_mixture <- function(component_elbo,
                                     prior_parameter,
                                     posterior_parameter_init = NULL,
                                     max_iters = 20,
                                     tol = 1e-10) {
  N = nrow(component_elbo)
  S = ncol(component_elbo)

  if (length(prior_parameter) == 1) {
    prior_parameter = rep(prior_parameter, S)
  }
  if (length(prior_parameter) != S ||
      any(!is.finite(prior_parameter)) || any(prior_parameter <= 0)) {
    stop("prior_parameter must contain one positive value per cluster.")
  }

  if (is.null(posterior_parameter_init)) {
    posterior_parameter = prior_parameter + N / S
  } else {
    if (length(posterior_parameter_init) != S ||
        any(!is.finite(posterior_parameter_init)) ||
        any(posterior_parameter_init <= 0)) {
      stop(
        "posterior_parameter_init must contain one positive value per cluster."
      )
    }
    posterior_parameter = posterior_parameter_init
  }

  omega = matrix(1 / S, nrow = N, ncol = S)
  for (iter in seq_len(max_iters)) {
    expected_log_pi = dirichlet_expected_log(posterior_parameter)
    log_responsibility = component_elbo + matrix(
      expected_log_pi, nrow = N, ncol = S, byrow = TRUE
    )
    omega_new = softmax_rows(log_responsibility)
    posterior_parameter_new = prior_parameter + colSums(omega_new)

    relative_change = max(
      abs(posterior_parameter_new - posterior_parameter) /
        pmax(posterior_parameter, 1)
    )
    omega = omega_new
    posterior_parameter = posterior_parameter_new
    if (relative_change < tol) break
  }

  list(
    omega = omega,
    posterior_parameter = posterior_parameter,
    expected_log_pi = dirichlet_expected_log(posterior_parameter),
    posterior_mean_pi = posterior_parameter / sum(posterior_parameter),
    allocation_mass = colMeans(omega),
    iterations = iter
  )
}

dirichlet_mixture_elbo <- function(component_elbo, omega,
                                   posterior_parameter,
                                   prior_parameter,
                                   eps = 1e-12) {
  N = nrow(component_elbo)
  S = ncol(component_elbo)
  expected_log_pi = dirichlet_expected_log(posterior_parameter)

  expected_component_fit = sum(omega * component_elbo)
  expected_cluster_log_probability = sum(
    omega * matrix(expected_log_pi, nrow = N, ncol = S, byrow = TRUE)
  )
  cluster_entropy = -sum(omega * log(pmax(omega, eps)))

  expected_component_fit + expected_cluster_log_probability +
    cluster_entropy - dirichlet_kl(
      parameter = posterior_parameter,
      prior_parameter = prior_parameter
    )
}

miso_fixed_gamma_dirichlet_pi <- function(
    Y, F, gamma_bar, max_iters = 50,
    n_inner = 5, update_prior = TRUE,
    prior_shape = 1, prior_beta = 1,
    alpha0_init = NULL, beta0_init = NULL,
    pi_init = NULL, omega_init = NULL,
    dirichlet_prior = 1,
    mixture_max_iters = 20, mixture_tol = 1e-10,
    rho_prior = NULL,
    tol = 1e-5, min_iters = 5,
    patience = 3, update_F = FALSE,
    update_gamma = TRUE,
    update_slot_scale = FALSE,
    slot_scale_init = NULL,
    slot_scale_step_init = 0.2,
    slot_scale_step_ramp = 10,
    min_slot_scale = 1e-3,
    max_slot_scale = 1e3,
    gamma_step_init = 0.5,
    gamma_step_ramp = 10,
    F_step_init = 0.2,
    F_step_ramp = 20,
    F_pseudocount = .Machine$double.eps,
    block_size = 100, eps = 1e-12) {
  N = nrow(Y)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = nrow(F)

  if (length(dirichlet_prior) == 1) {
    dirichlet_prior = rep(dirichlet_prior, S)
  }
  if (length(dirichlet_prior) != S ||
      any(!is.finite(dirichlet_prior)) || any(dirichlet_prior <= 0)) {
    stop("dirichlet_prior must contain one positive value per cluster.")
  }

  F = normalize_rows(F, eps)
  gamma_bar = normalize_gamma_bar(gamma_bar, eps)
  if (is.null(rho_prior)) rho_prior = array(1 / K, dim = c(S, D, K))
  if (is.null(slot_scale_init)) {
    slot_scale = matrix(1, S, D)
  } else {
    slot_scale = slot_scale_init
  }

  if (xor(is.null(alpha0_init), is.null(beta0_init))) {
    stop("alpha0_init and beta0_init must be supplied together.")
  }
  if (is.null(alpha0_init)) {
    alpha0 = matrix(prior_shape, nrow = S, ncol = D)
    beta0 = matrix(prior_beta, nrow = S, ncol = D)
  } else {
    if (!all(dim(alpha0_init) == c(S, D)) ||
        !all(dim(beta0_init) == c(S, D))) {
      stop("alpha0_init and beta0_init must both be S by D matrices.")
    }
    if (any(!is.finite(alpha0_init)) || any(alpha0_init <= 0) ||
        any(!is.finite(beta0_init)) || any(beta0_init <= 0)) {
      stop("Initial Gamma shapes and rates must be finite and positive.")
    }
    alpha0 = alpha0_init
    beta0 = beta0_init
  }
  alpha0_initial = alpha0
  beta0_initial = beta0

  if (is.null(omega_init)) {
    if (is.null(pi_init)) {
      omega = matrix(1 / S, nrow = N, ncol = S)
    } else {
      if (length(pi_init) != S || any(!is.finite(pi_init)) ||
          any(pi_init < 0) || sum(pi_init) <= 0) {
        stop("pi_init must be a nonnegative vector of length S.")
      }
      pi_init = pi_init / sum(pi_init)
      omega = matrix(pi_init, nrow = N, ncol = S, byrow = TRUE)
    }
  } else {
    if (!all(dim(omega_init) == c(N, S)) ||
        any(!is.finite(omega_init)) || any(omega_init < 0) ||
        any(rowSums(omega_init) <= 0)) {
      stop("omega_init must be a nonnegative N by S responsibility matrix.")
    }
    omega = omega_init / rowSums(omega_init)
  }

  omega_initial = omega
  dirichlet_posterior = dirichlet_prior + colSums(omega)
  dirichlet_posterior_initial = dirichlet_posterior
  pi_initial = dirichlet_posterior / sum(dirichlet_posterior)
  alpha = NULL
  elbo = rep(NA_real_, max_iters)
  gamma_history = vector("list", max_iters)
  small_improve_count = 0
  converged = FALSE
  n_iter = max_iters
  component_fit = NULL

  for (iter in seq_len(max_iters)) {
    gamma_step = min(1, gamma_step_init + (1 - gamma_step_init) *
                       (iter - 1) / max(gamma_step_ramp - 1, 1))
    F_step = min(1, F_step_init + (1 - F_step_init) *
                   (iter - 1) / max(F_step_ramp - 1, 1))
    slot_scale_step = min(1, slot_scale_step_init +
                            (1 - slot_scale_step_init) * (iter - 1) /
                            max(slot_scale_step_ramp - 1, 1))

    component_fit = fit_miso_components(
      Y = Y,
      F = F,
      gamma_bar = gamma_bar,
      alpha0 = alpha0,
      beta0 = beta0,
      slot_scale = slot_scale,
      alpha_init = alpha,
      n_inner = n_inner,
      block_size = block_size,
      omega = omega,
      compute_C = TRUE,
      compute_component_elbo = TRUE,
      eps = eps
    )
    alpha = component_fit$alpha

    mixture_fit = update_dirichlet_mixture(
      component_elbo = component_fit$component_elbo,
      prior_parameter = dirichlet_prior,
      posterior_parameter_init = dirichlet_posterior,
      max_iters = mixture_max_iters,
      tol = mixture_tol
    )
    omega = mixture_fit$omega
    dirichlet_posterior = mixture_fit$posterior_parameter

    component_fit = fit_miso_components(
      Y = Y,
      F = F,
      gamma_bar = gamma_bar,
      alpha0 = alpha0,
      beta0 = beta0,
      slot_scale = slot_scale,
      alpha_init = alpha,
      n_inner = 1,
      block_size = block_size,
      omega = omega,
      compute_C = TRUE,
      compute_component_elbo = TRUE,
      eps = eps
    )
    alpha = component_fit$alpha

    if (update_prior) {
      prior_fit = update_miso_priors(
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        alpha0 = alpha0,
        beta0 = beta0,
        eps = eps
      )
      alpha0 = prior_fit$alpha0
      beta0 = prior_fit$beta0
    }

    if (update_gamma) {
      gamma_bar = update_miso_gamma(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        slot_scale = slot_scale,
        rho_prior = rho_prior,
        gamma_step = gamma_step,
        eps = eps
      )
    }

    if (update_slot_scale) {
      slot_scale = update_miso_slot_scale(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        slot_scale = slot_scale,
        slot_scale_step = slot_scale_step,
        min_scale = min_slot_scale,
        max_scale = max_slot_scale,
        eps = eps
      )
    }

    if (update_F) {
      F = update_miso_F(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        F_step = F_step,
        F_pseudocount = F_pseudocount,
        eps = eps
      )
    }

    elbo[iter] = dirichlet_mixture_elbo(
      component_elbo = component_fit$component_elbo,
      omega = omega,
      posterior_parameter = dirichlet_posterior,
      prior_parameter = dirichlet_prior,
      eps = eps
    ) - miso_gamma_kl(gamma_bar, rho_prior = rho_prior, eps = eps)
    gamma_history[[iter]] = gamma_bar

    if (iter > 1) {
      rel_improve = (elbo[iter] - elbo[iter - 1]) /
        (abs(elbo[iter - 1]) + 1)
      if (iter >= min_iters && rel_improve < tol) {
        small_improve_count = small_improve_count + 1
      } else {
        small_improve_count = 0
      }
      if (iter >= min_iters && small_improve_count >= patience) {
        converged = TRUE
        n_iter = iter
        break
      }
    }
  }

  component_fit = fit_miso_components(
    Y = Y,
    F = F,
    gamma_bar = gamma_bar,
    alpha0 = alpha0,
    beta0 = beta0,
    slot_scale = slot_scale,
    alpha_init = alpha,
    n_inner = n_inner,
    block_size = block_size,
    omega = omega,
    compute_C = TRUE,
    compute_component_elbo = TRUE,
    eps = eps
  )
  mixture_fit = update_dirichlet_mixture(
    component_elbo = component_fit$component_elbo,
    prior_parameter = dirichlet_prior,
    posterior_parameter_init = dirichlet_posterior,
    max_iters = mixture_max_iters,
    tol = mixture_tol
  )
  omega = mixture_fit$omega
  dirichlet_posterior = mixture_fit$posterior_parameter
  pi = mixture_fit$posterior_mean_pi

  list(
    omega = omega,
    z_hat = max.col(omega),
    pi = pi,
    allocation_mass = mixture_fit$allocation_mass,
    expected_log_pi = mixture_fit$expected_log_pi,
    dirichlet_prior = dirichlet_prior,
    dirichlet_posterior = dirichlet_posterior,
    dirichlet_posterior_initial = dirichlet_posterior_initial,
    gamma_bar = gamma_bar,
    alpha = component_fit$alpha,
    beta_sd = component_fit$beta_sd,
    slot_scale = slot_scale,
    alpha0 = alpha0,
    beta0 = beta0,
    alpha0_initial = alpha0_initial,
    beta0_initial = beta0_initial,
    omega_initial = omega_initial,
    pi_initial = pi_initial,
    component_elbo = component_fit$component_elbo,
    C = component_fit$C,
    elbo = elbo,
    converged = converged,
    n_iter = n_iter,
    gamma_history = gamma_history[seq_len(n_iter)],
    F = F
  )
}
