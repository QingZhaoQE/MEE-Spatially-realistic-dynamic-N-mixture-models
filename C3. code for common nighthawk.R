code <- nimbleCode({
      
  # Priors
  ### Process
  for (k in 1:2) {
    beta0[k] ~ dnorm(0, sd = 2)
  } # k

  for (k in 1:4) {
    beta_rho[k] ~ dnorm(0, sd = 2)
  } # k

  log_kappa ~ dnorm(0, sd = 2)
  kappa <- exp(log_kappa)

  beta_eta ~ dnorm(0, sd = 2)

  log_delta ~ dnorm(0, sd = 2)
  delta <- exp(log_delta)

  ### IMBCR data
  sigma_imbcr ~ dgamma(0.01, 0.01)
  theta_imbcr ~ dunif(0, 1)

  ### BBS data
  chi_route_sd ~ dgamma(0.01, 0.01)
  chi_obser_sd ~ dgamma(0.01, 0.01)
  chi_stop_beta[1] ~ dnorm(0, sd=2)
  chi_stop_beta[2] ~ T(dnorm(0, sd=2), 0, )
  chi_mu ~ dnorm(0, sd=2)
  chi_new ~ dnorm(0, sd=2)
  chi_delta <- 0.001

  # Ecological processes
  for (i in 1:nsite) {
    lambda[i,1] <- exp(beta0[1] + beta0[2] * x[i,1])
    N[i,1] ~ dpois(lambda[i,1])

    for (t in 2:nyear) {
      rho[i,t-1] <- exp(
        beta_rho[1] + 
        beta_rho[2] * (N[i,t-1] - exp(beta0[1])) / exp(beta0[1]) + 
        beta_rho[3] * x[i,t] + 
        beta_rho[4] * (N[i,t-1] - exp(beta0[1])) / exp(beta0[1]) * x[i,t] )

      for (j in 1:nsite) {
        eta[i,j,t-1] <- exp(
          -1 * kappa * dist[i,j] + 
          beta_eta * (x[j,t] - x[i,t]) ) 
      } # j
      theta[i,1:nsite,t-1] <- eta[i,1:nsite,t-1] / sum(eta[i,1:nsite,t-1])

      lambda[i,t] <- sum(N[1:nsite,t-1] * rho[1:nsite,t-1] * theta[1:nsite,i,t-1]) + delta
      N[i,t] ~ dpois(lambda[i,t])
    } # t
  } # i

  # Observational processes
  for (d in 1:ndist_imbcr) {
    intval_imbcr[d] <- sigma_imbcr^2 * (exp(-1*(breaks_imbcr[d]^2)/(2*sigma_imbcr^2)) - exp(-1*(breaks_imbcr[d+1]^2)/(2*sigma_imbcr^2)))
    psi_imbcr[d] <- 2 * intval_imbcr[d] / (cutoff_imbcr ^ 2)
  } # d
  psi_imbcr_sum <- sum(psi_imbcr[1:ndist_imbcr])
  for(d in 1:ndist_imbcr) {
    psi_imbcr_prop[d] <- psi_imbcr[d] / psi_imbcr_sum
  } # d

  for (r in 1:ntime_imbcr) {
    phi_imbcr[r] <- (1 - theta_imbcr) ^ (r - 1) * theta_imbcr
  } # r
  phi_imbcr_sum <- sum(phi_imbcr[1:ntime_imbcr])
  for (r in 1:ntime_imbcr) {
    phi_imbcr_prop[r] <- phi_imbcr[r] / phi_imbcr_sum
  } # r

  for (d in 1:ndist_imbcr) {
    for (r in 1:ntime_imbcr) {
      pi_imbcr[(r-1)*ndist_imbcr+d] <- psi_imbcr_prop[d] * phi_imbcr_prop[r]
    } # r
  } # d

  for(k in 1:nobs_imbcr) {
    imbcr_count_sum[k] ~ dbinom(0.0490874 * npoint_imbcr[k] * psi_imbcr_sum * phi_imbcr_sum, N[sites_imbcr[k], years_imbcr[k]])
    imbcr_count[k,1:(ndist_imbcr*ntime_imbcr)] ~ dmultinom(pi_imbcr[1:(ndist_imbcr*ntime_imbcr)], imbcr_count_sum[k])
  } # k

  # BBS data
  for (k in 1:nroute_bbs) {
    chi_route_epsilon[k] ~ dnorm(0, sd=chi_route_sd)
  } # k
  for (k in 1:nobser_bbs) {
    chi_obser_epsilon[k] ~ dnorm(0, sd=chi_obser_sd)
  } # k
  for (j in 1:nstop_bbs) {
    chi_stop_sd[j] <- exp(chi_stop_beta[1] + chi_stop_beta[2] * stops_bbs[j])
  } # j
  for (k in 1:nobs_bbs) {
    log_chi_mu[k] <- 
      chi_mu + 
      chi_route_epsilon[route_bbs[k]] + 
      chi_obser_epsilon[obser_bbs[k]] + 
      chi_new * new_bbs[k]
    for (j in 1:nstop_bbs) {
      log_chi[k,j] ~ dnorm(log_chi_mu[k], chi_stop_sd[j])
      chi[k,j] <- exp(log_chi[k,j])
      bbs_count[k,j] ~ dpois(N[sites_bbs[k], years_bbs[k]] * chi[k,j] + chi_delta)
    } # j
  } # k

})
