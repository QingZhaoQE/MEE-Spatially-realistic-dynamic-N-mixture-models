code <- nimbleCode({
      
  # Priors
  for (k in 1:2) {
    beta0[k] ~ dnorm(0, sd = 2)
  } # k

  for (k in 1:2) {
    beta_rho_int[k] ~ dnorm(0, sd = 2)
  } # k
  beta_rho_ddp ~ dnorm(0, sd = 2)

  for (k in 1:2) {
    log_kappa[k] ~ dnorm(0, sd = 2)
    kappa[k] <- exp(log_kappa[k])
  } # k
  beta_eta_ddp ~ dnorm(0, sd = 2)

  alpha_int ~ dnorm(0, sd = 2)
  alpha_tmp ~ dnorm(0, sd = 2)
  alpha_wnd ~ dnorm(0, sd = 2)

  # Ecological processes
  for (i in 1:nsite) {
    lambda[i,1] <- exp(beta0[x[i]])
    N[i,1] ~ dpois(lambda[i,1])

    for (t in 2:nyear) {
      rho[i,t-1] <- exp(
        beta_rho_int[x[i]] + 
        beta_rho_ddp * (N[i,t-1] - mean(exp(beta0[1:2]))) / mean(exp(beta0[1:2])) / a[i] )

      for (j in 1:nsite) {
        eta[i,j,t-1] <- exp(
          -1 * kappa[x[i]] * d[i,j] + 
          beta_eta_ddp * (N[i,t-1] - mean(exp(beta0[1:2]))) / mean(exp(beta0[1:2])) / a[i] ) 
      } # j
      theta[i,1:nsite,t-1] <- eta[i,1:nsite,t-1] / sum(eta[i,1:nsite,t-1])

      lambda[i,t] <- sum(N[1:nsite,t-1] * rho[1:nsite,t-1] * theta[1:nsite,i,t-1])
      N[i,t] ~ dpois(lambda[i,t])
    } # t
  } # i

  # Observational processes
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        p[i,t,j] <- ilogit(alpha_int + alpha_tmp * w1[i,t,j] + alpha_wnd * w2[i,t,j])
        y1[i,t,j] ~ dbinom(p[i,t,j], N[i,t])
        y2[i,t,j] ~ dbinom(p[i,t,j], N[i,t])
      } # j
    } # t
  } # i

})
