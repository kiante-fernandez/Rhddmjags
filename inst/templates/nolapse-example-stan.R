# nolapse_test_rstan.R - Testing STAN fits of HDDM models without lapse process in STAN using RStan in R
#
# Copyright (C) 2023 Kianté Fernandez, <kiantefernan@gmail.com>
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
# 2023/02/11      Kianté  Fernandez                   intial code

library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
library(rstan)
library(bayesplot)

source(here("R", "Rhddmjagsutils.R")) #make sure you have correct installs

options(mc.cores=4)

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
# Note you could remove this if statement and replace with loading your own data "gendata"

if (!file.exists(here("data", "genparam_test.RData"))) {
  
  # Number of simulated participants
  nparts <- 10
  
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
  y <- rep(0, N)
  rt <- rep(0, N)
  acc <- rep(0, N)
  participant <- rep(0, N) # Participant index
  condition <- rep(0, N) # Condition index
  indextrack <- seq_len(ntrials)
  
  for (p in seq_len(nparts)) {
    for (k in seq_len(nconds)) {
      tempout <- simulratcliff(
        N = ntrials, Alpha = alpha[[p]], Tau = ter[[p]], Beta = beta[[p]],
        Nu = delta[[p, k]], Eta = deltatrialsd[[p]], rangeTau = tertrialrange[[p]]
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
  genparam$ter <- ter
  genparam$beta <- beta
  genparam$alpha <- alpha
  genparam$delta <- delta
  genparam$tertrialrange <- tertrialrange
  genparam$deltatrialsd <- deltatrialsd
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

# STAN code

# Set random seed
set.seed(2023)

# Running stan code
modelfile <- here("stancode", "nolapse_test.stan")

nchains <- 3
burnin <- 250
nsamps <- 5000

model = stan_model(modelfile)

# Track these variables
stan_params <- c(
  "deltasdcond",
  "tersd", "alphasd", "betasd", "deltasd",
  "terhier", "alphahier", "betahier", "deltahier",
  "ter", "alpha", "beta", "deltapart",
  "delta"
)

# Fit model to data
y <- genparam$y
rt <- genparam$rt
participant <- genparam$participant
condition <- genparam$condition
nparts <- genparam$nparts
nconds <- genparam$nconds
ntrials <- genparam$ntrials
N <- genparam$N

minrt <- rep(0, nparts)

datalist <- list(
  y <- y,
  N <- N,
  nparts <- nparts,
  nconds <- nconds,
  condition <- condition,
  participant <- participant
)
for (p in seq_len(nparts)) {
  minrt[[p]] <- min(rt[(participant == p)])
}
# get names for the list
names(datalist) <- c("y", "N", "nparts", "nconds", "condition", "participant")

initials <- vector(mode = "list")
for (c in seq_len(nchains)) {
  initsList <- function() {
    chaininit <- vector(mode = "list")
    chaininit$deltasdcond <- runif(n = 1, .1, 3)
    chaininit$tersd <- runif(1, .01, .2)
    chaininit$alphasd <- runif(1, .01, 1.)
    chaininit$betasd <- runif(1, .01, .2)
    chaininit$deltasd <- runif(1, .1, 3)
    chaininit$deltapart <- runif(nparts, -4., 4.)
    chaininit$delta <- matrix(runif(nparts * nconds, -4., 4.), nrow = nparts, ncol = nconds)
    chaininit$ter <- runif(nparts, .1, .5)
    chaininit$alpha <- runif(nparts, .5, 2.)
    chaininit$beta <- runif(nparts, .2, .8)
    chaininit$deltahier <- runif(1, -4., 4.)
    chaininit$terhier <- runif(1, .1, .5)
    chaininit$alphahier <- runif(1, .5, 2.)
    chaininit$betahier <- runif(1, .2, .8)
    for (p in seq_len(nparts)) {
      chaininit$ter[[p]] <- runif(1, 0, minrt[[p]] / 2)
    }
    return(chaininit)
  }
  initials[[c]] <- initsList()
}

fit = sampling(model,
               data = datalist,
               pars = stan_params,
               warmup = burnin,
               iter=nsamps,
               chains=nchains,
               thin=10, 
               init = initials,
               seed = 2023)


