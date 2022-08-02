# recovery_test.R - Testing JAGS fits of HDDM models in JAGS using R2jags in R
#                    This hierarchical model is usuaully appropriate for trials with intermixed experimental conditions
#
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
# This R code was generated using Michael D. Nunez's `recovery_test.py`
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 01/08/2022    Kianté Fernandez                        Original code generation



# Libraries
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
library(R2jags) # jags.parallel is part of R2jags

source(here("R", "Rhddmjagsutils.R"))

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
# Note you could remove this if statement and replace with loading your own data and nameing it `genparam`

if (!file.exists(here("data", "genparam_test.RData"))) {
  
  # Number of simulated participants
  nparts <- 10
  
  # Number of conditions
  nconds <- 3
  
  # Number of trials per participant and condition
  ntrials <- 50
  
  # Number of total trials in each simulation
  N <- ntrials * nparts * nconds
  
  # Set random seed
  set.seed(2022)
  
  ndt <- runif(n = nparts, min = .15, max = .6) # Uniform from .15 to .6 seconds
  alpha <- runif(nparts, .8, 1.4) # Uniform from .8 to 1.4 evidence units
  beta <- runif(nparts, .3, .7) # Uniform from .3 to .7 * alpha
  delta <- matrix(runif(nparts * nconds, -4, 4), nrow = nparts, ncol = nconds) # Uniform from -4 to 4 evidence units per second
  ndttrialrange <- runif(n = nparts, 0, .1) # Uniform from 0 to .1 seconds
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
      tempout <- simulratcliff(
        N = ntrials, Alpha = alpha[[p]], Tau = ndt[[p]], Beta = beta[[p]],
        Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = ndttrialrange[[p]]
      )
      tempx <- sign(Re(tempout))
      tempt <- abs(Re(tempout))
      mindwanderx <- sample(0:2, ntrials, replace = T) * 2 - 1
      mindwandert <- runif(n = ntrials, 0, 2) # Randomly distributed from 0 to 2 seconds
      
      mindwander_trials <- sample(1:ntrials, prob = rep(prob_lapse[[p]] / 100, ntrials), replace = F)
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
  genparam$ndt <- ndt
  genparam$beta = beta
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$ndttrialrange <- ndttrialrange
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
  save(genparam, file = here("data", "genparam_test.RData"))
} else {
  # load dataset
  load(here("data", "genparam_test.RData"))
}


# Set random seed
set.seed(2022)

# JAGS code
tojags <- "model {
    
    ##########
    #Between-condition variability parameters priors
    ##########

    #Between-condition variability in drift rate to choice A
    deltasdcond ~ dgamma(1,1)

    ##########
    #Between-participant variability parameters priors
    ##########

    #Between-participant variability in non-decision time
    tersd ~ dgamma(.3,1)

    #Between-participant variability in Speed-accuracy trade-off
    alphasd ~ dgamma(1,1)

    #Between-participant variability in choice A start point bias
    betasd ~ dgamma(.3,1)

    #Between-participant variability in lapse trial probability
    problapsesd ~ dgamma(.3,1)

    #Between-participant variability in drift rate to choice A
    deltasd ~ dgamma(1,1)

    ##########
    #Hierarchical DDM parameter priors
    ##########

    #Hierarchical Non-decision time
    terhier ~ dnorm(.5, pow(.25,-2))

    #Hierarchical boundary parameter (speed-accuracy tradeoff)
    alphahier ~ dnorm(1, pow(.5,-2))

    #Hierarchical start point bias towards choice A
    betahier ~ dnorm(.5, pow(.25,-2))

    #Hierarchical lapse trial probability
    problapsehier ~ dnorm(.3, pow(.15,-2))

    #Hierarchical drift rate to choice A
    deltahier ~ dnorm(0, pow(2, -2))

    ##########
    #Participant-level DDM parameter priors
    ##########
    for (p in 1:nparts) {

        #Non-decision time
        ter[p] ~ dnorm(terhier, pow(tersd,-2))T(0, 1)

        #Boundary parameter (speed-accuracy tradeoff)
        alpha[p] ~ dnorm(alphahier, pow(alphasd,-2))T(0, 3)

        #Rightward start point bias towards choice A
        beta[p] ~ dnorm(betahier, pow(betasd,-2))T(0, 1)

        #Probability of a lapse trial
        problapse[p] ~ dnorm(problapsehier, pow(problapsesd,-2))T(0, 1)
        probDDM[p] <- 1 - problapse[p]

        #Participant-level drift rate to choice A
        deltapart[p] ~ dnorm(deltahier, pow(deltasd, -2))

        for (c in 1:nconds) {

            #Participant-level drift rate to choice A
            delta[p,c] ~ dnorm(deltapart[p], pow(deltasdcond, -2))

        }

    }

    ##########
    # Wiener likelihood and uniform mixture using Ones trick
    for (i in 1:N) {

        # Log density for DDM process of rightward/leftward RT
        ld_comp[i, 1] <- dlogwiener(y[i], alpha[participant[i]], ter[participant[i]], beta[participant[i]], delta[participant[i],condition[i]])

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

writeLines(tojags, here("jagscode", "recovery_test.jags"))

nchains <- 3
burnin <- 200
nsamps <- 1000

modelfile <- here("jagscode", "recovery_test.jags")

# Track these variables
jags_params <- c("deltasdcond", 
                 "tersd", "alphasd","betasd","problapsesd","deltasd",
                 "terhier","alphahier","betahier","problapsehier","deltahier",
                 "ter","alpha","beta","problapse","deltapart",
                "delta","DDMorLapse")

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

minrt <- rep(0, nparts)
for (p in seq_len(nparts)) {
  minrt[[p]] <- min(rt[(participant == p)])
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
names(datalist) <- c("y","N","nparts","nconds","condition","participant","Ones","Constant")


initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$deltasdcond <- runif(1, .1, 3)
    chaininit$tersd <- runif(1, .01, .2)
    chaininit$alphasd <- runif(1, .01, 1.)
    chaininit$betasd <- runif(1, .01, .2)
    chaininit$problapsesd <- runif(1, .01, .5)
    chaininit$deltasd <- runif(1, .1, 3.)
    chaininit$deltapart <- runif(nparts, -4., 4.)
    chaininit$delta <- matrix(runif(nparts * nconds, -4., 4.), nrow = nparts, ncol = nconds)
    chaininit$ter <- runif(nparts, .1, .5)
    chaininit$alpha <- runif(nparts, .5, 2.)
    chaininit$beta <- runif(nparts, .2, .8)
    chaininit$problapse <- runif(nparts, .01, .1)
    chaininit$deltahier <- runif(1, -4., 4.)
    chaininit$terhier <- runif(1, .1, .5)
    chaininit$alphahier <- runif(1, .5, 2.)
    chaininit$betahier <- runif(1, .2, .8)
    chaininit$problapsehier <- runif(1, .01, .1)
    for (p in seq_len(nparts)) {
      chaininit$ter[[p]] <- runif(1, 0, minrt[[p]] / 2)
    }
    return(chaininit)
  }
  initials[[c]] <- initsList()
}

print(paste0("Fitting model ..."))

samples <- jags(
  model.file = modelfile,
  data = datalist, inits = initials, jags_params,
  n.iter = nsamps,
  n.chains = nchains,
  n.burnin = burnin, 
  jags.module = "wiener"
)
savestring <- here("modelfits", "genparam_test_model.Rdata")
print(paste0("Saving results to: ", savestring))
save(samples, file = savestring)

# Diagnostics
# for now just call the jags object Diagnostics() function soon to come!
samples

# Posterior distributions

# jellyfish(samples, "delta")
# # 
# jellyfish(samples, "ter")
# # 
# jellyfish(samples, "alpha")
# 
# # Recovery
# recovery(samples, genparam["delta"])
# 
# recovery(jagsfit, genparam["ndt"])
# 
# recovery(jagsfit, genparam["alpha"])
# 
# recovery(samples, genparam["delta_int"])
# 
# recovery(jagsfit, genparam["ndt_int"])
# 
# # recovery(samples, genparam["alpha_int"])
# 
# recovery(jagsfit, genparam["delta_gamma"])
# 
# recovery(samples, genparam["ndt_gamma"])
# 
# recovery(samples, genparam["alpha_gamma"])


