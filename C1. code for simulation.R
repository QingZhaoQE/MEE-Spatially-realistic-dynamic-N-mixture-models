library(boot)
library(nimble)
library(nimbleHMC)

#===============
# Simulate data
#===============
# Basic values
nsite <- 100  # number of sites
nyear <- 10   # number of years, 20, 10, 5
nreps <- 3    # number of within-season replicates
pmiss <- 0.50 # proportion of unsurveyed sites, 0, 0.50, 0.67, 0.90

beta0 <- c(1.8, 1)               # intercept and slopes for initial abundance
beta_rho <- c(0, -0.8, 0.6, 0.2) # intercept and slopes for population growth
kappa <- 2                       # distance effect on movement
beta_eta <- c(-0.8, 0.6, 0.2)    # slopes for movement
alpha <- c(-0.9, -0.5)            # intercept and slopes for detection probability

# Environmental covariates
x <- matrix(rnorm(nsite * nyear, 0, 1), nsite, nyear)

# Locations of sites
lon <- runif(nsite, -10, 10) # longitude or easting
lat <- runif(nsite, - 8,  8) # latitude or northing
d <- matrix(, nsite, nsite) # distance between sites
for (i in 1:nsite) {
  for (j in 1:nsite) {
    d[i,j] <- sqrt((lon[i] - lon[j]) ^ 2 + (lat[i] - lat[j]) ^ 2)
  } # j
} # i

# Ecological processes
lambda <- matrix(, nsite, nyear) # expectation of abundance
N <- matrix(, nsite, nyear) # abundance

lambda[,1] <- exp(cbind(1,x[,1]) %*% beta0) # expectation of initial abundance
N[,1] <- rpois(nsite, lambda[,1]) # initial abundance

rho <- matrix(0, nsite, nyear-1) # population growth rate
eta <- array(, dim=c(nsite, nsite, nyear-1)) # unstandardized colonization rate
theta <- array(, dim=c(nsite, nsite, nyear-1)) # standardized colonization rate
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t], (N[,t-1]-lambda[,1])/lambda[,1]*x[,t]) %*% beta_rho)
  for (i in 1:nsite) {
    eta[i,,t-1] <- exp(-1 * kappa * d[i,] + cbind((N[,t-1]-N[i,t-1])/exp(beta0[1]), x[,t]-x[i,t], (N[,t-1]-N[i,t-1])/exp(beta0[1])*(x[,t]-x[i,t])) %*% beta_eta)
  } # i
  theta[,,t-1] <- eta[,,t-1] / rowSums(eta[,,t-1])
  lambda[,t] <- (N[,t-1] * rho[,t-1]) %*% theta[,,t-1]
  N[,t] <- rpois(nsite, lambda[,t])
} # t

# Observational processes
w <- array(rnorm(nsite * nyear * nreps, 0, 1), dim=c(nsite, nyear, nreps)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    p[i,t,] <- inv.logit(cbind(1,w[i,t,]) %*% alpha)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps)) # count data
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y[i,t,] <- rbinom(nreps, N[i,t], p[i,t,])
  } # t
} # i

smiss <- sort(sample(1:nsite, nsite*pmiss, replace=F)) # unsurveyed sites
y[smiss,,] <- NA # create missing values for unsurveyed sites

#========================
# Define model in Nimble
#========================
code <- nimbleCode({
      
  # Priors
  for (k in 1:2) {
    beta0[k] ~ dnorm(0, sd = 2)
  } # k

  for (k in 1:4) {
    beta_rho[k] ~ dnorm(0, sd = 2)
  } # k

  log_kappa ~ dnorm(0, sd = 2)
  kappa <- exp(log_kappa)

  for (k in 1:3) {
    beta_eta[k] ~ dnorm(0, sd = 2)
  } # k

  for (k in 1:2) {
    alpha[k] ~ dnorm(0, sd = 2)
  } # k

  # Ecological processes
  for (i in 1:nsite) {
    lambda[i,1] <- exp(beta0[1] + beta0[2] * x[i,1])
    N[i,1] ~ dpois(lambda[i,1])

    for (t in 2:nyear) {
      rho[i,t-1] <- exp(
        beta_rho[1] + 
        beta_rho[2] * (N[i,t-1]-lambda[i,1])/lambda[i,1] + 
        beta_rho[3] * x[i,t] + 
        beta_rho[4] * (N[i,t-1]-lambda[i,1])/lambda[i,1] * x[i,t] )

      for (j in 1:nsite) {
        eta[i,j,t-1] <- exp(
          -1 * kappa * d[i,j] + 
          beta_eta[1] * (N[j,t-1]-N[i,t-1])/exp(beta0[1]) + 
          beta_eta[2] * (x[j,t]-x[i,t]) + 
          beta_eta[3] * (N[j,t-1]-N[i,t-1])/exp(beta0[1]) * (x[j,t]-x[i,t]) ) 
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
        p[i,t,j] <- ilogit(alpha[1] + alpha[2] * w[i,t,j])
        y[i,t,j] ~ dbinom(p[i,t,j], N[i,t])
      } # j
    } # t
  } # i

}) # nimbleCode

#============
# Run Nimble
#============
# Data
constants <- list(
  nsite=nsite, nyear=nyear, nreps=nreps, 
  x=x, d=d, w=w
)

data <- list(
  y=y
)

# Initial values
Ni <- apply(y, 1:2, max) + 5
Ni[which(is.na(Ni))] <- 5
inits <- list(beta0=rep(0,2), beta_rho=rep(0,4), log_kappa=0, beta_eta=rep(0,3), alpha=rep(0,2), N=Ni)

model <- nimbleModel(code, constants=constants, data=data, inits=inits, buildDerivs=TRUE)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers(c('beta0', 'beta_rho', 'log_kappa', 'beta_eta', 'alpha', 'N[1,1]', 'N[10,5]'))

mcmcConf$removeSamplers(c('beta0', 'beta_rho', 'log_kappa', 'beta_eta', 'alpha'))
mcmcConf$printSamplers(c('beta0', 'beta_rho', 'log_kappa', 'beta_eta', 'alpha', 'N[1,1]', 'N[10,5]'))

mcmcConf$addSampler(target = c('beta0'), type = "NUTS")
mcmcConf$addSampler(target = c('beta_rho'), type = "NUTS")
mcmcConf$addSampler(target = c('log_kappa', 'beta_eta'), type = "NUTS")
mcmcConf$addSampler(target = c('alpha'), type = "NUTS")
mcmcConf$printSamplers(c('beta0', 'beta_rho', 'log_kappa', 'beta_eta', 'alpha', 'N[1,1]', 'N[10,5]'))

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

chain <- 3
nmcmc <- 2000
nburn <- 1000
nthin <- 1
fit <- runMCMC(compiled$mcmc, nchains = chain, niter = nmcmc, nburnin = nburn, thin=nthin)


