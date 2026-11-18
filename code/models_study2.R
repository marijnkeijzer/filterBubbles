# Model functions extracted without changes to their statistical specifications.
bayes_m0 <- function(y_cens, x1, x2, 
                     burn_in=4000, sample_size=2000, thinning=5){
  model <- "
  model {
    for ( i in 1:Ntotal ) {
      mu[i] = x1[i] + alpha * (x2[i] - x1[i])
      y_cens[i] ~ dnorm( mu[i] , 1/sigma^2 )T(0,1)
      a1[i] <- pnorm(0, mu[i], 1/sigma^2)
      a2[i] <- pnorm(1, mu[i], 1/sigma^2)
      a[i] <- a2[i]-a1[i]
      loglik[i] <- logdensity.norm(y_cens[i], mu[i], 1/sigma^2) - log(a[i])  
    }
    alpha ~ dnorm( mu_alpha , 1 / (sigma_alpha)^2 )
    sigma ~ dunif( sdY/1000 , sdY*1000 )
  }
"
  
  Ntotal = length(y_cens)
  
  # Specify the data in a list, for later shipment to JAGS:
  dat = dump.format(
    list(
      y_cens = y_cens ,
      x1 = x1 ,
      x2 = x2 ,
      Ntotal = Ntotal ,
      mu_alpha = 0 ,
      sigma_alpha = 1 ,
      sdY = sd(y_cens,na.rm=TRUE)
    )
  )
  
  # Initial values of MCMC chains based on data:
  alpha = 0
  sigma = sd(y_cens,na.rm=TRUE) 
  
  # # initial values of censored data:
  # yInit = rep( 0.5 , length(y_cens) )
  
  initsList = list( alpha = alpha , sigma = sigma ) # , y=yInit )
  #combine with different random number generators and seeds
  inits1 <- dump.format(c(initsList, list(.RNG.name="base::Super-Duper", .RNG.seed=99999)))
  inits2 <- dump.format(c(initsList, list(.RNG.name="base::Wichmann-Hill", .RNG.seed=1234)))
  inits3 <- dump.format(c(initsList, list(.RNG.name="base::Mersenne-Twister", .RNG.seed=6666)))
  
  monitor = c( "alpha" , "sigma" , "loglik")
  results <- run.jags(model=model,#"CODE/ANALYSIS/models/weightedLinearInfluence.txt", 
                      monitor=monitor, data=dat, 
                      n.chains=3, inits=c(inits1, inits2, inits3), 
                      plots = FALSE, burnin=burn_in, sample=sample_size, thin=thinning)
  
  return(results)
}


bayes_m1 <- function(y_cens, x1, x2, 
                     burn_in=4000, sample_size=2000, thinning=5){
  model <- "
  model {
    for ( i in 1:Ntotal ) {
      mu[i] = x1[i] + alpha * (x2[i] - x1[i]) - gamma * abs(x2[i] - x1[i]) * (x2[i] - x1[i])
      y_cens[i] ~ dnorm( mu[i] , 1/sigma^2 )T(0,1)
      a1[i] <- pnorm(0, mu[i], 1/sigma^2)
      a2[i] <- pnorm(1, mu[i], 1/sigma^2)
      a[i] <- a2[i]-a1[i]
      loglik[i] <- logdensity.norm(y_cens[i], mu[i], 1/sigma^2) - log(a[i])  
    }
    alpha ~ dnorm( mu_alpha , 1 / (sigma_alpha)^2 )
    gamma ~ dnorm( mu_gamma , 1 / (sigma_gamma)^2 )
    sigma ~ dunif( sdY/1000 , sdY*1000 )
  }
"
  
  Ntotal = length(y_cens)
  
  # Specify the data in a list, for later shipment to JAGS:
  dat = dump.format(
    list(
      y_cens = y_cens ,
      x1 = x1 ,
      x2 = x2 ,
      Ntotal = Ntotal ,
      mu_alpha = 0 ,
      sigma_alpha = 1 ,
      mu_gamma = 0 ,
      sigma_gamma = 1 ,
      sdY = sd(y_cens,na.rm=TRUE)
    )
  )
  
  # Initial values of MCMC chains based on data:
  alpha = 0
  gamma = 0 
  sigma = sd(y_cens,na.rm=TRUE) 
  
  # # initial values of censored data:
  # yInit = rep( 0.5 , length(y_cens) )
  
  initsList = list( alpha = alpha , gamma = gamma , sigma = sigma ) # , y=yInit )
  #combine with different random number generators and seeds
  inits1 <- dump.format(c(initsList, list(.RNG.name="base::Super-Duper", .RNG.seed=99999)))
  inits2 <- dump.format(c(initsList, list(.RNG.name="base::Wichmann-Hill", .RNG.seed=1234)))
  inits3 <- dump.format(c(initsList, list(.RNG.name="base::Mersenne-Twister", .RNG.seed=6666)))
  
  monitor = c( "alpha" , "gamma" , "sigma" , "loglik")
  results <- run.jags(model=model,#"CODE/ANALYSIS/models/weightedLinearInfluence.txt", 
                      monitor=monitor, data=dat, 
                      n.chains=3, inits=c(inits1, inits2, inits3), 
                      plots = FALSE, burnin=burn_in, sample=sample_size, thin=thinning)
  
  return(results)
}

bayes_m2 <- function(y_cens, x1, x2, g1, 
                           burn_in=4000, sample_size=2000, thinning=5){
  model <- "
  model {
    for ( i in 1:Ntotal ) {
      mu[i] = x1[i] + alpha * (x2[i] - x1[i]) 
                    - gamma * abs(x2[i] - x1[i]) * (x2[i] - x1[i]) 
                    - zeta * g1[i] * abs(x2[i] - x1[i]) * (x2[i] - x1[i])
      y_cens[i] ~ dnorm( mu[i] , 1/sigma^2 )T(0,1)
      a1[i] <- pnorm(0, mu[i], 1/sigma^2)
      a2[i] <- pnorm(1, mu[i], 1/sigma^2)
      a[i] <- a2[i]-a1[i]
      loglik[i] <- logdensity.norm(y_cens[i], mu[i], 1/sigma^2) - log(a[i])  
    }
    alpha ~ dnorm( mu_alpha , 1 / (sigma_alpha)^2 )
    gamma ~ dnorm( mu_gamma , 1 / (sigma_gamma)^2 )
    zeta ~  dnorm( mu_zeta  , 1 / (sigma_zeta)^2  )
    sigma ~ dunif( sdY/1000 , sdY*1000 )
  }
"
  
  Ntotal = length(y_cens)
  
  # Specify the data in a list, for later shipment to JAGS:
  dat = dump.format(
    list(
      y_cens = y_cens ,
      x1 = x1 ,
      x2 = x2 ,
      g1 = g1 ,
      Ntotal = Ntotal ,
      mu_alpha = 0 ,
      sigma_alpha = 1 ,
      mu_gamma = 0 ,
      sigma_gamma = 1 ,
      mu_zeta = 0 ,
      sigma_zeta = 1 ,
      sdY = sd(y_cens,na.rm=TRUE)
    )
  )
  
  # Initial values of MCMC chains based on data:
  alpha = 0
  gamma = 0 
  zeta = 0
  sigma = sd(y_cens,na.rm=TRUE) 
  
  # # initial values of censored data:
  # yInit = rep( 0.5 , length(y_cens) )
  
  initsList = list( alpha = alpha , gamma = gamma , sigma = sigma , zeta = zeta) # , y=yInit )
  #combine with different random number generators and seeds
  inits1 <- dump.format(c(initsList, list(.RNG.name="base::Super-Duper", .RNG.seed=99999)))
  inits2 <- dump.format(c(initsList, list(.RNG.name="base::Wichmann-Hill", .RNG.seed=1234)))
  inits3 <- dump.format(c(initsList, list(.RNG.name="base::Mersenne-Twister", .RNG.seed=6666)))
  
  monitor = c( "alpha" , "gamma" , "sigma" , "zeta" , "loglik")
  results <- run.jags(model=model,#"CODE/ANALYSIS/models/weightedLinearInfluence.txt", 
                      monitor=monitor, data=dat, 
                      n.chains=3, inits=c(inits1, inits2, inits3), 
                      plots = FALSE, burnin=burn_in, sample=sample_size, thin=thinning)
  
  return(results)
}

