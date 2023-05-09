# simple_test.R - Testing JAGS fits of a non-hierarchical DDM model without lapse process in JAGS using R2jags in R
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
# This R code was generated using Michael D. Nunez's `simple_test.py`
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
# Note you could remove this if statement and replace with loading your own data to dictionary "gendata"

if (!file.exists(here("data", "simpleparam_test.RData"))) {
  
  # Number of simulated participants
  nparts <- 100
  
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
  y <- rep(0, N)
  rt <- rep(0, N)
  acc <- rep(0, N)
  participant <- rep(0, N) # Participant index
  indextrack <- seq_len(ntrials)
  for (p in seq_len(nparts)) {
    tempout <- simulratcliff(
      N = ntrials, Alpha = alpha[[p]], Tau = ndt[[p]], Beta = beta[[p]],
      Nu = delta[[p]], Eta = deltatrialsd[[p]]
    )
    tempx <- sign(Re(tempout))
    tempt <- abs(Re(tempout))
    y[indextrack] <- tempx * tempt
    rt[indextrack] <- tempt
    acc[indextrack] <- (tempx) / 2
    participant[indextrack] <- p
    indextrack <- indextrack + ntrials
  }
  
  genparam <- vector(mode = "list")
  genparam$ndt <- ndt
  genparam$beta <- beta
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$deltatrialsd <- deltatrialsd
  genparam$rt <- rt
  genparam$acc <- acc
  genparam$y <- y
  genparam$participant <- participant
  genparam$nparts <- nparts
  genparam$ntrials <- ntrials
  genparam$N <- N
  save(genparam, file = here("data", "simpleparam_test.RData"))
} else {
  # load dataset
  load(here("data", "simpleparam_test.RData"))
}

# JAGS code

# Set random seed
set.seed(2022)

tojags <- "
model {
    
    ##########
    #Simple DDM parameter priors
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

    }

    ##########
    # Wiener likelihood
    ##########
    for (i in 1:N) {

        # Observations of accuracy*RT for DDM process of rightward/leftward RT
        y[i] ~ dwiener(alpha[participant[i]], ndt[participant[i]], beta[participant[i]], delta[participant[i]])

    }
}
"

# Rjags code

load.module("wiener")
load.module("dic")
list.modules()

writeLines(tojags, here("jagscode", "simple_test.jags"))

nchains <- 6
burnin <- 2000
nsamps <- 10000

modelfile <- here("jagscode", "simple_test.jags")

# Track these variables
jags_params <- c("alpha", "ndt", "beta", "delta")

# Fit model to data
N <- genparam$N

y <- genparam$y
rt <- genparam$rt
participant <- genparam$participant
nparts <- genparam$nparts
ntrials <- genparam$ntrials

minrt <- rep(0, nparts)

datalist <- list(
  y <- y,
  N <- N,
  nparts <- nparts,
  participant <- participant
)
for (p in seq_len(nparts)) {
  minrt[[p]] <- min(rt[(participant == p)])
}
# get names for the list
names(datalist) <- c("y", "N","nparts","participant")

# initialize initial values
initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$alpha <- runif(nparts, .5, 2.)
    chaininit$ndt <- runif(nparts, .1, .5)
    chaininit$beta <- runif(nparts, .2, .8)
    chaininit$delta <- runif(nparts, -4., 4.)
    for (p in seq_len(nparts)) {
      chaininit$ndt[[p]] <- runif(1, 0, minrt[[p]] / 2)
    }
    return(chaininit)
  }
  initials[[c]] <- initsList()
}
print(paste0("Fitting ", "simple", " model ..."))
jagsfit <- R2jags::jags(
  model.file = modelfile,
  data = datalist, inits = initials, jags_params,
  n.iter = nsamps,
  n.chains = nchains,
  n.burnin = burnin, jags.module = "wiener"
)
samples <- update(jagsfit, n.iter = nsamps)
savestring <- here("modelfits", "simple_test_simple.Rdata")
print(paste0("Saving results to: ", savestring))
save(samples, file = savestring)

# Diagnostics
diags <- diagnostic(samples)

# Posterior distributions
jellyfish(samples, "alpha",filename = "figures/alpha_posteriors_simple.png")

jellyfish(samples, "ndt","figures/ndt_posteriors_simple.png")

jellyfish(samples, "beta","figures/beta_posteriors_simple.png")

jellyfish(samples, "delta","figures/delta_posteriors_simple.png")

# Recovery
recovery(samples, genparam["alpha"], "figures/alpha_recovery_simple.png")

recovery(samples, genparam["ndt"],"figures/ndt_recovery_simple.png")

recovery(samples, genparam["beta"],"figures/beta_recovery_simple.png")

recovery(samples, genparam["delta"],"figures/delta_recovery_simple.png")

