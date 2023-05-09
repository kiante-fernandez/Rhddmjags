# regression_test.R - Testing JAGS fits of HDDM models with participant-level regressors in JAGS using Rjags in R
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
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 11/07/22      Kianté  Fernandez                      Starting coding
# 13/07/22      Kianté  Fernandez                      Completed Buggs code
# 20/07/22      Kianté  Fernandez                      added jelly and recovery plots
# 01/08/22      Kianté  Fernandez                      fixed lapse trial bug
# 02/08/22      Kianté  Fernandez                      added diagnostic


# Libraries
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
library(R2jags) # jags.parallel is part of R2jags

source(here("R", "Rhddmjagsutils.R"))

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
# Note you could remove this if statement and replace with loading your own data and nameing it `genparam`

if (!file.exists(here("data", "genparam_reg_test.RData"))) {

  # Number of simulated participants
  nparts <- 40

  # Number of conditions
  nconds <- 6

  # Number of trials per participant and condition
  ntrials <- 50

  # Number of total trials in each simulation
  N <- ntrials * nparts * nconds

  # Set random seed
  set.seed(2022)

  # Intercepts of linear regressions
  ndt_int <- matrix(rep(runif(nconds, .4, .7), nparts), nrow = nparts, ncol = nconds) # Uniform from .4 to .7 seconds
  alpha_int <- matrix(rep(runif(nconds, .8, 1.4), nparts), nrow = nparts, ncol = nconds) # Uniform from .8 to 1.4 evidence units
  delta_int <- matrix(rep(runif(nconds, -2, 2), nparts), nrow = nparts, ncol = nconds) # Uniform from -2 to 2 evidence units per second

  # Slopes of linear regressions
  ndt_gamma <- matrix(rep(runif(nconds, 0, .1), nparts), nrow = nparts, ncol = nconds)
  alpha_gamma <- matrix(rep(runif(nconds, -.1, .1), nparts), nrow = nparts, ncol = nconds)
  delta_gamma <- matrix(rep(runif(nconds, -1, 1), nparts), nrow = nparts, ncol = nconds)

  # Regressors
  regressors1 <- matrix(rnorm(nparts * nconds), nrow = nparts, ncol = nconds) # The same regressors for each parameter, from a standard normal distribution

  # True parameters
  ndt <- ndt_int + ndt_gamma * regressors1
  alpha <- alpha_int + alpha_gamma * regressors1
  delta <- delta_int + delta_gamma * regressors1


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
        N = ntrials, Alpha = alpha[[p, k]], Tau = ndt[[p, k]],
        Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = ndttrialrange[[p]]
      )
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
  genparam$ndt <- ndt
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$ndt_int <- ndt_int
  genparam$alpha_int <- alpha_int
  genparam$delta_int <- delta_int
  genparam$ndt_gamma <- ndt_gamma
  genparam$alpha_gamma <- alpha_gamma
  genparam$delta_gamma <- delta_gamma
  genparam$regressors1 <- regressors1
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
  save(genparam, file = here("data", "genparam_reg_test.RData"))
} else {
  # load dataset
  load(here("data", "genparam_reg_test.RData"))
}

# Fit model to data
y <- genparam$y
rt <- genparam$rt
participant <- genparam$participant
condition <- genparam$condition
nparts <- genparam$nparts
nconds <- genparam$nconds
regressors1 <- genparam$regressors1
ntrials <- genparam$ntrials
N <- genparam$N

minrt <- matrix(rep(0, nparts * nconds), nrow = nparts, ncol = nconds)

for (p in seq_len(nparts)) {
  for (k in seq_len(nconds)) {
    minrt[[p, k]] <- min(rt[(participant == p) & (condition == k)])
  }
}

# Input for mixture modeling
Ones <- rep(1, N)
Constant <- 20

datalist <- list(
  y <- y,
  N <- N,
  nparts <- nparts,
  nconds <- nconds,
  condition <- condition,
  participant <- participant,
  regressors1 <- regressors1,
  Ones <- Ones,
  Constant <- Constant
)

# get names for the list
names(datalist) <- c("y", "N", "nparts", "nconds", "condition", "participant", "regressors1", "Ones", "Constant")

# Set random seed
set.seed(2022)

# JAGS code
tojags <- "model {

    ##########
    #Between-condition variability priors
    ##########

    #Between-condition variability in drift rate to correct
    deltasdcond ~ dgamma(1,1)

    #Between-condition variability in non-decision time
    ndtsdcond ~ dgamma(.3,1)

    #Between-condition variability in speed-accuracy trade-off
    alphasdcond ~ dgamma(1,1)

    ##########
    #Between-participant variability priors
    ##########

    #Between-participant variability in lapse trial probability
    problapsesd ~ dgamma(.3,1)

    ##########
    #Hierarchical DDM parameter priors
    ##########

    #Hierarchical lapse trial probability
    problapsehier ~ dnorm(.3, pow(.15,-2))

    ##########
    #Condition-level DDM parameter priors
    ##########

    for (c in 1:nconds) {

        #Drift rate intercept
        delta_int[c] ~ dnorm(0, pow(6, -2))

        #Non-decision time intercept
        ndt_int[c] ~ dnorm(0, pow(2,-2))

        #Boundary parameter intercept
        alpha_int[c] ~ dnorm(0, pow(4,-2))

        #Effect of regressor1 on Drift rate
        delta_gamma[c] ~ dnorm(0, pow(3, -2))

        #Effect of regressor1 on Non-decision time
        ndt_gamma[c] ~ dnorm(0, pow(1,-2))

        #Effect of regressor1 on boundary parameter
        alpha_gamma[c] ~ dnorm(0, pow(2,-2))


    }


    ##########
    #Participant-level DDM parameter priors
    ##########
    for (p in 1:nparts) {

        #Probability of a lapse trial
        problapse[p] ~ dnorm(problapsehier, pow(problapsesd,-2))T(0, 1)
        probDDM[p] <- 1 - problapse[p]

        for (c in 1:nconds) {

            #Participant-level drift rate to correct
            delta[p,c] ~ dnorm(delta_int[c] + delta_gamma[c]*regressors1[p,c], pow(deltasdcond, -2))

            #Non-decision time
            ndt[p,c] ~ dnorm(ndt_int[c] + ndt_gamma[c]*regressors1[p,c], pow(ndtsdcond,-2))T(0, 1)

            #Boundary parameter (speed-accuracy tradeoff)
            alpha[p,c] ~ dnorm(alpha_int[c] + alpha_gamma[c]*regressors1[p,c], pow(alphasdcond,-2))T(0, 3)

        }

    }

    ##########
    # Wiener likelihood and uniform mixture using Ones trick
    for (i in 1:N) {

        # Log density for DDM process of rightward/leftward RT
        ld_comp[i, 1] <- dlogwiener(y[i], alpha[participant[i],condition[i]], ndt[participant[i],condition[i]], .5, delta[participant[i],condition[i]])

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

writeLines(tojags, here("jagscode", "regression_test.jags"))

nchains <- 6
burnin <- 4000
nsamps <- 20000

modelfile <- here("jagscode", "regression_test.jags")

# Track these variables
jags_params <- c(
  "deltasdcond", "ndtsdcond", "alphasdcond", "problapsesd",
  "problapsehier", "delta_int", "ndt_int", "alpha_int",
  "delta_gamma", "ndt_gamma", "alpha_gamma",
  "delta", "ndt", "alpha", "problapse", "DDMorLapse"
)

initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$deltasdcond <- runif(1, .1, 3)
    chaininit$ndtsdcond <- runif(1, .01, .2)
    chaininit$alphasdcond <- runif(1, .01, 1)
    chaininit$problapsesd <- runif(1, .01, .5)
    chaininit$problapsehier <- runif(1, .01, .1)
    chaininit$delta_int <- runif(nconds, -4., 4.)
    chaininit$ndt_int <- runif(nconds, .1, .5)
    chaininit$alpha_int <- runif(nconds, .5, 2)
    chaininit$delta_gamma <- runif(nconds, -1, 1)
    chaininit$ndt_gamma <- runif(nconds, -1, 1)
    chaininit$alpha_gamma <- runif(nconds, -1, 1)
    chaininit$delta <- matrix(runif(nparts * nconds, -4., 4.), nrow = nparts, ncol = nconds)
    chaininit$ndt <- matrix(runif(nparts * nconds, .1, .5), nrow = nparts, ncol = nconds)
    chaininit$alpha <- matrix(runif(nparts * nconds, .5, 2.), nrow = nparts, ncol = nconds)
    chaininit$problapse <- runif(nparts, .01, .1)

    for (p in seq_len(nparts)) {
      for (k in seq_len(nconds)) {
        chaininit$ndt[[p, k]] <- runif(1, 0, minrt[[p, k]] / 2)
      }
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
  n.thin = 10,
  n.burnin = burnin, jags.module = "wiener"
)
savestring <- here("modelfits", "genparam_regression_test.Rdata")
print(paste0("Saving results to: ", savestring))
save(samples, file = savestring)

# Diagnostics
diags <- diagnostic(samples, exclude = "DDMorLapse")

# Posterior distributions
jellyfish(samples, "delta", "figures/delta_posteriors_model.png")

jellyfish(samples, "ndt", "figures/ndt_posteriors_model.png")

jellyfish(samples, "alpha", "figures/alpha_posteriors_model.png")

# Recovery
recovery(samples, genparam["delta"], "figures/delta_recovery_model.png")

recovery(samples, genparam["ndt"], "figures/ndt_recovery_model.png")

recovery(samples, genparam["alpha"], "figures/alpha_recovery_model.png")

recovery(samples, genparam["delta_int"], "figures/delta_int_recovery_model.png")

recovery(samples, genparam["ndt_int"], "figures/ndt_int_recovery_model.png")

recovery(samples, genparam["alpha_int"], "figures/alpha_int_recovery_model.png")

recovery(samples, genparam["delta_gamma"], "figures/delta_gamma_recovery_model.png")

recovery(samples, genparam["ndt_gamma"], "figures/ndt_gamma_recovery_model.png")

recovery(samples, genparam["alpha_gamma"], "figures/alpha_gamma_recovery_model.png")
