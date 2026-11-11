library(LaMa)

# AR(1) transition matrix
ar1_transition <- function(phi, sigma_s, grid) {
  m <- length(grid)
  boundaries <- c(-Inf, (grid[-1] + grid[-m]) / 2, Inf)
  gamma <- matrix(0, m, m)
  
  for (i in seq_len(m)) {
    gamma[i, ] <- diff(
      pnorm(boundaries, mean = phi * grid[i], sd = sigma_s)
    )
  }
  
  gamma
}

# Observation densities
observation_densities <- function(y, grid, eta, sigma, pi_zero) {
  allprobs <- matrix(0, length(y), length(grid))
  positive <- y > 0
  
  for (j in seq_along(grid)) {
    allprobs[!positive, j] <- pi_zero[!positive]
    allprobs[positive, j] <- (1 - pi_zero[positive]) *
      dnorm(log(y[positive]), eta[positive] + grid[j], sigma) /
      y[positive]
  }
  
  allprobs
}

# General negative log-likelihood
nll_ssm_general <- function(par, y, grid, x = NULL, z = NULL) {
  if (!all(is.finite(par))) return(1e12)
  
  mu <- par[2]
  sigma <- exp(par[3])
  phi <- tanh(par[4])
  sigma_s <- exp(par[5])
  
  if (sigma < 0.01 || sigma_s < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) return(1e12)
  
  n_beta <- if (is.null(x)) 0 else ncol(as.matrix(x))
  n_alpha <- if (is.null(z)) 0 else ncol(as.matrix(z))
  
  beta <- if (n_beta > 0) par[6:(5 + n_beta)] else numeric(0)
  alpha <- if (n_alpha > 0) {
    par[(6 + n_beta):(5 + n_beta + n_alpha)]
  } else numeric(0)
  
  eta <- rep(mu, length(y))
  if (n_beta > 0) eta <- eta + as.vector(as.matrix(x) %*% beta)
  
  logit_pi <- rep(par[1], length(y))
  if (n_alpha > 0) {
    logit_pi <- logit_pi + as.vector(as.matrix(z) %*% alpha)
  }
  
  pi_zero <- plogis(logit_pi)
  
  nll <- tryCatch({
    gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(gamma)
    
    allprobs <- observation_densities(
      y, grid, eta, sigma, pi_zero
    )
    
    -LaMa::forward(
      delta = delta,
      Gamma = gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(nll)) return(1e12)
  nll
}

# Baseline model
nll_ssm <- function(par, y, grid) {
  nll_ssm_general(par, y, grid)
}

# Model with volume covariates
nll_ssm_cov <- function(par, y, x, grid) {
  nll_ssm_general(par, y, grid, x = x)
}

# Model with volume and zero-probability covariates
nll_ssm_zero_cov <- function(par, y, x, z, grid) {
  nll_ssm_general(par, y, grid, x = x, z = z)
}
