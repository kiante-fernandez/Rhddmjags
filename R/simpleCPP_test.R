# simpleCPP_test.R- Testing JAGS fits of a non-hierarchical Neural-DDM model
# in JAGS using Rjags in R assumes CPP slopes are generated from drift rates
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
# This R code was generated using Michael D. Nunez's `simpleCPP_test.py`
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 20/07/2022    Kianté Fernandez                        Original code generation


# Libraries
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
library(R2jags) # jags.parallel is part of R2jags

source(here("R", "Rhddmjagsutils.R"))

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
# Note you could remove this if statement and replace with loading your own data to dictionary "gendata"

if (!file.exists(here("data", "genparam_test.RData"))) {
  
  # Number of simulated participants
  nparts <- 10
  
  # Number of trials per participant and condition
  ntrials <- 100
  
  # Number of total trials in each simulation
  N <- ntrials * nparts
  
  # Set random seed
  set.seed(2022)
  
  ndt <- runif(n = nparts, min = .15, max = .6) # Uniform from .15 to .6 seconds
  alpha <- runif(nparts, .8, 1.4) # Uniform from .8 to 1.4 evidence units
  beta <- runif(nparts, .3, .7) # Uniform from .3 to .7 * alpha
  delta <- runif(nparts, -4, 4) # Uniform from -4 to 4 evidence units per second
  deltatrialsd <- runif(nparts, 0, 2) # Uniform from 0 to 2 evidence units per second
  CPPnoise = runif(nparts, 0, 1) # Uniform from 0 to 1 evidence units per second
  y <- rep(0, N)
  rt <- rep(0, N)
  acc <- rep(0, N)
  CPP <-  rep(0, N)
  participant <- rep(0, N) # Participant index
  indextrack <- seq_len(ntrials)
  for (p in seq_len(nparts)) {
    for (k in seq_len(nconds)) {
      tempout <- simulratcliff(
        N = ntrials, Alpha = alpha[[p]], Tau = ndt[[p]], Beta = beta[[p]],
        Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = ndttrialrange[[p]]
      )
      tempx <- sign(Re(tempout))
      tempt <- abs(Re(tempout))
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
  genparam$beta <- beta
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$deltatrialsd <- deltatrialsd
  genparam$CPPnoise <- CPPnoise
  genparam$CPP <- CPP
  genparam$rt <- rt
  genparam$acc <- acc
  genparam$y <- y
  genparam$participant <- participant
  genparam$nparts <- nparts
  genparam$ntrials <- ntrials
  genparam$N <- N
  save(genparam, file = here("data", "genparam_test.RData"))
} else {
  # load dataset
  load(here("data", "genparam_test.RData"))
}

# JAGS code

# Set random seed
set.seed(2022)

tojags = "
model {
    
    ##########
    #Simple NDDM parameter priors
    ##########
    for (p in 1:nparts) {
    
        #Boundary parameter (speed-accuracy tradeoff) per participant
        alpha[p] ~ dnorm(1, pow(.5,-2))T(0, 3)

        #Non-decision time per participant
        ndt[p] ~ dnorm(.5, pow(.25,-2))T(0, 1)

        #Start point bias towards choice A per participant
        beta[p] ~ dnorm(.5, pow(.25,-2))T(0, 1)

        #Drift rate to choice A per participant
        delta[p] ~ dnorm(0, pow(2, -2))

        #Noise in observed EEG measure, the CentroParietal Positivity (CPP) slope per participant
        CPPnoise[p] ~ dnorm(1, pow(.5,-2))T(0, 3)

    }

    ##########
    # Wiener likelihood
    for (i in 1:N) {

        # Observations of accuracy*RT for DDM process of rightward/leftward RT
        y[i] ~ dwiener(alpha[participant[i]], ndt[participant[i]], beta[participant[i]], delta[participant[i]])

        # Observations of CentroParietal Positivity (CPP) slope per trial
        CPP[i] ~ dnorm(delta[participant[i]],pow(CPPnoise[participant[i]],-2))

    }
}
"

# Rjags code

load.module("wiener")
load.module("dic")
list.modules()

writeLines(tojags, here("jagscode", "simpleCPP_test.jags"))

nchains <- 3
burnin <- 250
nsamps <- 5000

modelfile <- here("jagscode", "simpleCPP_test.jags")

# Track these variables
jags_params <- c(
  "alpha", "ndt", "beta", "delta", "CPPnoise"
)

# Fit model to data
N <- genparam$N

y <- genparam$y
rt <- genparam$rt
CPP <- genparam$CPP
participant <- genparam$participant
nparts <- genparam$nparts
ntrials <- genparam$ntrials

minrt <- rep(0, nparts)

datalist <- list(
  y <- y,#note are you feeding the correct thing here?
  N <- N,
  CPP <- CPP,
  nparts <- nparts,
  participant <- participant
)
for (p in seq_len(nparts)) {
  minrt[[p]] <- min(rt[(participant == p)])
}
# get names for the list
names(datalist) <- c("y", "N","CPP", "nparts", "participant")

#initialize initial values
initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$alpha <- runif(nparts, .5, 2.)
    chaininit$ndt <- runif(nparts, .1, .5)
    chaininit$beta <- runif(nparts, .2, .8)
    chaininit$delta <- runif(nparts, -4., 4.)
    chaininit$CPPnoise <- runif(nparts, .5, 2.)
    for (p in seq_len(nparts)) {
      chaininit$ndt[[p]] <- runif(1, 0, minrt[[p]] / 2)
    }
    return(chaininit)
  }
  initials[[c]] <- initsList()
}

print(paste0("Fitting ", "simpleEEG", " model ..."))

jagsfit <- R2jags::jags(
  model.file = modelfile,
  data = datalist, inits = initials, jags_params,
  n.iter = nsamps,
  n.chains = nchains,
  n.burnin = burnin, jags.module = "wiener"
)

samples <- update(jagsfit, n.iter = nsamps)

savestring <- here("modelfits", "simpleCPP_test.Rdata")
print(paste0("Saving results to: ", savestring))

save(samps, file = savestring)

# Diagnostics

# Posterior distributions


