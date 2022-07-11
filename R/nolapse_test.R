# nolapse_test.py - Testing JAGS fits of HDDM models without lapse process in JAGS using Rjags in R
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
# 06/07/22      Kianté  Fernandez                 Rewrote python code for generate data
# 11/07/22      Kianté  Fernandez                 started Jaggs code Recoding


# Libraries
library(dplyr) # A Grammar of Data Manipulation, CRAN v1.0.8
library(tidyr) # Tidy Messy Data, CRAN v1.2.0
library(readr) # Read Rectangular Text Data, CRAN v2.1.1
library(magrittr) # A Forward-Pipe Operator for R, CRAN v2.0.2
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics, CRAN v3.3.5
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
require(rjags) # Bayesian Graphical Models using MCMC, CRAN v4-12. NOTE: Must have previously installed package rjags.
source(here("R","Rhddmjagsutils.R"))

### Simulations ###

# Generate samples from the joint-model of reaction time and choice
#Note you could remove this if statement and replace with loading your own data "gendata"

if (!file.exists(here('data','genparam_test.RData'))) {
  
  # Number of simulated participants
  nparts <- 10
  
  # Number of conditions
  nconds <-  6
  
  # Number of trials per participant and condition
  ntrials <-  50
  
  # Number of total trials in each simulation
  N <-  ntrials*nparts*nconds
  
  # Set random seed
  set.seed(2022)
  
  ndt   <- runif(n = nparts, min = .15, max = .6) # Uniform from .15 to .6 seconds
  alpha <- runif(nparts, .8, 1.4) # Uniform from .8 to 1.4 evidence units
  beta  <- runif(nparts, .3, .7) # Uniform from .3 to .7 * alpha
  delta <- matrix(runif(nparts*nconds,-4, 4), nrow = nparts, ncol = nconds) # Uniform from -4 to 4 evidence units per second
  ndttrialrange <- runif(n = nparts, 0,.1) # Uniform from 0 to .1 seconds
  deltatrialsd <- runif(n = nparts, 0,2) # Uniform from 0 to 2 evidence units per second
  y <-  rep(0,N)
  rt <-   rep(0,N)
  acc <-   rep(0,N)
  participant <-  rep(0,N) #Participant index
  condition <-  rep(0,N) #Condition index
  indextrack <-  seq_len(ntrials)
  
  for (p in seq_len(nparts)){
    for (k in seq_len(nconds)){
      #tempout <- simulratcliff() # for testing
      tempout <- simulratcliff(N = ntrials, Alpha = alpha[[p]], Tau = ndt[[p]], Beta = beta[[p]],
                               Nu = delta[[p,k]], Eta = deltatrialsd[[p]], rangeTau = ndttrialrange[[p]])
      tempx <- sign(Re(tempout))
      tempt <- abs(Re(tempout))
      y[indextrack] <- tempx*tempt
      rt[indextrack] <- tempt
      acc[indextrack] <- (tempx + 1)/2 #do you need this 1 here?
      participant[indextrack] <- p
      condition[indextrack] <- k 
      indextrack <- indextrack + ntrials
      
  }
  } 
  
  genparam <- vector(mode = "list")
  genparam$ndt <- ndt
  genparam$beta <-  beta
  genparam$alpha <-  alpha
  genparam$delta <-  delta
  genparam$ndttrialrange <-  ndttrialrange
  genparam$deltatrialsd <-  deltatrialsd
  genparam$rt <- rt
  genparam$acc <-  acc
  genparam$y <-  y
  genparam$participant <-  participant
  genparam$condition <-  condition
  genparam$nparts <-  nparts
  genparam$nconds <-  nconds
  genparam$ntrials <-  ntrials
  genparam$N <-  N
  save(genparam, file = here("data", "genparam_test.RData"))
  
} else {
  #load dataset
  load(here("data", "genparam_test.RData"))
}

# JAGS code

# Set random seed
set.seed(2022)


    

