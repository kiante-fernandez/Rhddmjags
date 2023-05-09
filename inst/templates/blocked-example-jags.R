#blocked_exp_conds.R - Code to fit a HDDM with fixed start point from a blocked condition experiment using R2JAGS in R
# Copyright (C) 2022 Kianté Fernandez, <kiantefernan@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This R code was generated using Michael D. Nunez's `simple_test.py`
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 20/08/2022    Kianté Fernandez                        Original code
# 20/08/2022    Kianté Fernandez                        Plot and summary code

# Libraries
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
require(R2jags) # jags.parallel is part of R2jags

source(here("R", "Rhddmjagsutils.R")) #make sure you have correct installs

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
# Note you could remove this if statement and replace with loading your own data "gendata"

if (!file.exists(here("data", "blocked_genparam_test.RData"))) {
  
  # Number of simulated participants
  nparts <- 50
  
  # Number of conditions
  nconds <- 4
  
  # Number of trials per participant and condition
  ntrials <- 100
  
  # Number of total trials in each simulation
  N <- ntrials * nparts * nconds
  
  # Set random seed
  set.seed(2022)
  
  ter <- runif(n = nparts, min = .15, max = .6) # Uniform from .15 to .6 seconds
  alpha <- runif(nparts, .8, 1.4) # Uniform from .8 to 1.4 evidence units
  beta <- runif(nparts, .3, .7) # Uniform from .3 to .7 * alpha
  delta <- matrix(runif(nparts * nconds, -4, 4), nrow = nparts, ncol = nconds) # Uniform from -4 to 4 evidence units per second
  tertrialrange <- runif(n = nparts, 0, .1) # Uniform from 0 to .1 seconds
  deltatrialsd <- runif(n = nparts, 0, 2) # Uniform from 0 to 2 evidence units per second
  prob_lapse <- runif(n = nparts, 0, 10) # From 0 to 10 percent of trials
  y <- rep(0, N)
  rt <- rep(0, N)
  acc <- rep(0, N)
  participant <- rep(0, N) # Participant index
  condition <- rep(0, N) # Condition index
  indextrack <- seq_len(ntrials)
  for (p in seq_len(nparts)) {
    for (k in seq_len(nconds)) {
      #note difference ask MN if intentional in Py code
      tempout <- simulratcliff(
        N = ntrials, Alpha = alpha[[p]], Tau = ter[[p]], Beta = beta[[p]],
        Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = tertrialrange[[p]]
      )
      # tempout <- simulratcliff(
      #   N = ntrials, Alpha = alpha[[p, k]], Tau = ter[[p, k]],Beta = beta[[p]],
      #   Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = tertrialrange[[p]]
      # )
      tempx <- sign(Re(tempout))
      tempt <- abs(Re(tempout))
      mindwanderx <- sample(0:2, ntrials, replace = T) * 2 - 1
      mindwandert <- runif(n = ntrials, 0, 2) # Randomly distributed from 0 to 2 seconds
      
      mindwander_trials <- sample(1:ntrials, size = as.integer(round(ntrials * (prob_lapse[[p]] / 100))), replace = F)
      tempx[mindwander_trials] <- mindwanderx[mindwander_trials]
      tempt[mindwander_trials] <- mindwandert[mindwander_trials]
      y[indextrack] <- tempx * tempt
      rt[indextrack] <- tempt
      acc[indextrack] <- (tempx) / 2
      participant[indextrack] <- p
      condition[indextrack] <- k
      indextrack <- indextrack + ntrials
    }
  }
  
  genparam <- vector(mode = "list")
  genparam$ter <- ter
  genparam$beta <- beta
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$tertrialrange <- tertrialrange
  genparam$deltatrialsd <- deltatrialsd
  genparam$prob_lapse <- prob_lapse
  genparam$rt <- rt
  genparam$acc <- acc
  genparam$y <- y
  genparam$participant <- participant
  genparam$condition <- condition
  genparam$nparts <- nparts
  genparam$nconds <- nconds
  genparam$ntrials <- ntrials
  genparam$N <- N
  save(genparam, file = here("data", "blocked_genparam_test.RData"))
} else {
  # load dataset
  load(here("data", "blocked_genparam_test.RData"))
}


N <- genparam$N

# Input for mixture modeling
Ones <- rep(1, N)
Constant <- 20

# Fit model to data
y <- genparam$y
rt <- genparam$rt
participant <- genparam$participant
condition <- genparam$condition
nparts <- genparam$nparts
nconds <- genparam$nconds
ntrials <- genparam$ntrials

minrt <- matrix(rep(0, nparts * nconds), nrow = nparts, ncol = nconds)
for (p in seq_len(nparts)) {
  for (k in seq_len(nconds)) {
    minrt[[p, k]] <- min(rt[(participant == p) & (condition == k)])
  }
}

datalist <- list(
  y <- y,
  N <- N,
  nparts <- nparts,
  nconds <- nconds,
  condition <- condition,
  participant <- participant,
  Ones <- Ones,
  Constant <- Constant
)
# get names for the list
names(datalist) <- c("y", "N", "nparts", "nconds", "condition", "participant", "Ones", "Constant")

# JAGS code

# Set random seed
set.seed(2022)

tojags <- "model {
    
    ##########
    #Between-condition variability priors
    ##########

    #Between-condition variability in drift rate to correct
    deltasdcond ~ dgamma(1,1)

    #Between-condition variability in non-decision time
    tersdcond ~ dgamma(.3,1)

    #Between-condition variability in speed-accuracy trade-off
    alphasdcond ~ dgamma(1,1)

    ##########
    #Between-participant variability priors
    ##########

    #Between-participant variability in drift rate to correct
    deltasd ~ dgamma(1,1)

    #Between-participant variability in non-decision time
    tersd ~ dgamma(.3,1)

    #Between-participant variability in Speed-accuracy trade-off
    alphasd ~ dgamma(1,1)

    #Between-participant variability in lapse trial probability
    problapsesd ~ dgamma(.3,1)

    ##########
    #Hierarchical DDM parameter priors
    ##########

    #Hierarchical drift rate to correct
    deltahier ~ dnorm(0, pow(2, -2))

    #Hierarchical Non-decision time
    terhier ~ dnorm(.5, pow(.25,-2))

    #Hierarchical boundary parameter (speed-accuracy tradeoff)
    alphahier ~ dnorm(1, pow(.5,-2))

    #Hierarchical lapse trial probability
    problapsehier ~ dnorm(.3, pow(.15,-2))

    ##########
    #Participant-level DDM parameter priors
    ##########
    for (p in 1:nparts) {

        #Participant-level drift rate to correct
        deltapart[p] ~ dnorm(deltahier, pow(deltasd, -2))

        #Participant-level non-decision time
        terpart[p] ~ dnorm(terhier, pow(tersd,-2))

        #Participant-level boundary parameter (speed-accuracy tradeoff)
        alphapart[p] ~ dnorm(alphahier, pow(alphasd,-2))

        #Probability of a lapse trial
        problapse[p] ~ dnorm(problapsehier, pow(problapsesd,-2))T(0, 1)
        probDDM[p] <- 1 - problapse[p]


        ##########
        #Condition-level DDM parameter priors
        ##########
        for (c in 1:nconds) {

            #Drift rate to correct
            delta[p,c] ~ dnorm(deltapart[p], pow(deltasdcond, -2))

            #Non-decision time
            ter[p,c] ~ dnorm(terpart[p], pow(tersdcond,-2))T(0, 1)

            #Boundary parameter (speed-accuracy tradeoff)
            alpha[p,c] ~ dnorm(alphapart[p], pow(alphasdcond,-2))T(0, 3)

        }

    }

    ##########
    # Wiener likelihood and uniform mixture using Ones trick
    for (i in 1:N) {

        # Log density for DDM process of rightward/leftward RT
        ld_comp[i, 1] <- dlogwiener(y[i], alpha[participant[i],condition[i]], ter[participant[i],condition[i]], .5, delta[participant[i],condition[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -3, 3)

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        selected_density[i] <- exp(ld_comp[i, DDMorLapse[i]] - Constant)
        
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(selected_density[i])

        # Probability of mind wandering trials (lapse trials)
        DDMorLapse[i] ~ dcat( c(probDDM[participant[i]], problapse[participant[i]]) )
    }
}"

# Rjags code
load.module("wiener")
load.module("dic")
list.modules()

writeLines(tojags, here("jagscode", "blocked_exp_conds.jags"))

nchains <- 6
burnin <- 2000
nsamps <- 10000

modelfile <- here("jagscode", "blocked_exp_conds.jags")

# Track these variables
jags_params <- c('deltasdcond', 'tersdcond', 'alphasdcond',
  'deltasd', 'tersd', 'alphasd', 'problapsesd',
  'deltahier', 'terhier', 'alphahier', 'problapsehier',
  'deltapart', 'terpart', 'alphapart',
  'delta', 'ter', 'alpha', 'problapse', 'DDMorLapse')

initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$deltasdcond <- runif(n = 1, .1, 3)
    chaininit$tersdcond <- runif(n = 1, .01, .2)
    chaininit$alphasdcond <- runif(1, .01, 1)
    chaininit$deltasd <- runif(1, .1, 3)
    chaininit$tersd <- runif(1, .01, .2)
    chaininit$alphasd <- runif(1, .01, 1.)
    chaininit$problapsesd <- runif(1, .01, .5)
    chaininit$deltahier <- runif(1, -4., 4.)
    chaininit$terhier <- runif(1, .1, .5)
    chaininit$alphahier <- runif(1, .5, 2.)
    chaininit$problapsehier <- runif(1, .01, .1)
    chaininit$deltapart <- runif(nparts, -4., 4.)
    chaininit$terpart <- runif(nparts, .1, .5)
    chaininit$alphapart <- runif(nparts, .5, 2.)
    chaininit$problapse <- runif(nparts, .01, .1)
    chaininit$delta <- matrix(runif(nparts * nconds, -4., 4.), nrow = nparts, ncol = nconds)
    chaininit$ter <- matrix(runif(nparts * nconds, .1, .5), nrow = nparts, ncol = nconds)
    chaininit$alpha <- matrix(runif(nparts * nconds, .5, 2.), nrow = nparts, ncol = nconds)
    #chaininit$ter <- runif(nparts, .1, .5)
    #chaininit$alpha <- runif(nparts, .5, 2.)
    
    for (p in seq_len(nparts)) {
      for (k in seq_len(nconds)) {
        chaininit$ter[[p, k]] <- runif(1, 0, minrt[[p, k]] / 2)
      }
    }
    
    return(chaininit)
  }
  initials[[c]] <- initsList()
}

print(paste0("Fitting model 2 ..."))
samples <- jags(
  model.file = modelfile,
  data = datalist, inits = initials, jags_params,
  n.iter = nsamps,
  n.chains = nchains,
  n.thin = 10,
  n.burnin = burnin, 
  jags.module = "wiener"
)
savestring <- here("modelfits", "genparam_test_model2.Rdata")
print(paste0("Saving results to: ", savestring))
save(samples, file = savestring)

#Diagnostics
diags <- diagnostic(samples, exclude = c("DDMorLapse"))

#Parameter estimates
print(paste0('The median posterior drift-rate (evidence units / sec) for participant 1 in condition 1 is ', round(median(samples$BUGSoutput$sims.list$delta[,1,1]),4) ,      ' with 95%% credible interval (',round(quantile(samples$BUGSoutput$sims.list$delta[,1,1], c(2.5)/100),4), ',',round(quantile(samples$BUGSoutput$sims.list$delta[,1,1], c(97.5)/100),4), ')'))
print(paste0('The median posterior non-decision time (sec) for participant 1 in condition 1 is ', round(median(samples$BUGSoutput$sims.list$ter[,1,1]),4) ,      ' with 95%% credible interval (',round(quantile(samples$BUGSoutput$sims.list$ter[,1,1], c(2.5)/100),4), ',',round(quantile(samples$BUGSoutput$sims.list$ter[,1,1], c(97.5)/100),4), ')'))
print(paste0('The median posteriorboundary (evidence units) for participant 1 in condition 1 is ', round(median(samples$BUGSoutput$sims.list$alpha[,1,1]),4) ,      ' with 95%% credible interval (',round(quantile(samples$BUGSoutput$sims.list$alpha[,1,1], c(2.5)/100),4), ',',round(quantile(samples$BUGSoutput$sims.list$alpha[,1,1], c(97.5)/100),4), ')'))

#Posterior distributions
jellyfish(samples, "delta",filename = "figures/delta_posteriors_blocked_exp.png")

jellyfish(samples, "ter",filename = "figures/ter_posteriors_blocked_exp.png")

jellyfish(samples, "alpha",filename = "figures/alpha_posteriors_blocked_exp.png")

#Find posterior probability of each trial being from a lapse process



