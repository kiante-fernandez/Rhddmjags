# Rhddmjagsutils.py - Functions for simulation, model diagnostics, and parameter recovery
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
# 06/07/22      Kianté Fernandez                      Rewrote Michael python code in R
# 11/07/22      Kianté Fernandez                      Rewrote Joachim's translations

# Libraries

library(dplyr) # A Grammar of Data Manipulation, CRAN v1.0.8
library(tidyr) # Tidy Messy Data, CRAN v1.2.0
library(readr) # Read Rectangular Text Data, CRAN v2.1.1
library(magrittr) # A Forward-Pipe Operator for R, CRAN v2.0.2
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics, CRAN v3.3.5
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
require(rjags) # Bayesian Graphical Models using MCMC, CRAN v4-12. NOTE: Must have previously installed package rjags.


### Functions ###

# Simulate diffusion models slowly with intrinsic trial-to-trial variability in parameters
simul_ratcliff_slow <- function(N = 100, Alpha = 1, Tau = .4, Nu = 1, Beta = .5, rangeTau = 0, rangeBeta = 0, Eta = .3, Varsigma = 1, nsteps = 300, step_length = .01) {
  # SIMUL_RATCLIFF_SLOW  Generates data according to a drift diffusion model with optional trial-to-trial variability
  #
  # Parameters
  # ----------
  # N: a integer denoting the size of the output vector
  # (defaults to 100 experimental trials)
  #
  # Alpha: the mean boundary separation across trials  in evidence units
  # (defaults to 1 evidence unit)
  #
  # Tau: the mean non-decision time across trials in seconds
  # (defaults to .4 seconds)
  #
  # Nu: the mean drift rate across trials in evidence units per second
  # (defaults to 1 evidence units per second, restricted to -5 to 5 units)
  #
  # Beta: the initial bias in the evidence process for choice A as a proportion of boundary Alpha
  # (defaults to .5 or 50% of total evidence units given by Alpha)
  #
  # rangeTau: Non-decision time across trials is generated from a uniform
  # distribution of Tau - rangeTau/2 to  Tau + rangeTau/2 across trials
  # (defaults to 0 seconds)
  #
  # rangeZeta: Bias across trials is generated from a uniform distribution
  # of Zeta - rangeZeta/2 to Zeta + rangeZeta/2 across trials
  # (defaults to 0 evidence units)
  #
  # Eta: Standard deviation of the drift rate across trials
  # (defaults to 3 evidence units per second, restricted to less than 3 evidence units)
  #
  # Varsigma: The diffusion coefficient, the standard deviation of the
  # evidence accumulation process within one trial. It is recommended that
  # this parameter be kept fixed unless you have reason to explore this parameter
  # (defaults to 1 evidence unit per second)
  #
  # Returns
  # -------
  # Vector with reaction times (in seconds) multiplied by the response vector
  # such that negative reaction times encode response B and positive reaction times
  # encode response A

  if (Nu < -5 || Nu > 5) {
    Nu <- sign(Nu) * 5
    warning(paste0("Nu is not in the range [-5 5], bounding drift rate to", Nu))
  }
  if (Eta > 3) {
    warning(paste0("Standard deviation of drift rate is out of bounds, bounding drift rate to 3"))
    Eta <- 3
  }
  if (Eta == 0) {
    Eta <- 1e-16
  }

  # Initialize output vectors
  rts <- rep(0, N)
  choice <- rep(0, N)

  for (n in seq_len(N)) {
    random_walk <- vector(mode = "numeric", length = nsteps)
    start_point <- runif(1, Beta - rangeBeta / 2, Beta + rangeBeta / 2)
    ndt <- runif(1, Tau - rangeTau / 2, Tau + rangeTau / 2)
    drift <- rnorm(1, mean = Nu, sd = Eta)
    random_walk[[1]] <- start_point * Alpha
    for (s in 2:nsteps) {
      random_walk[[s]] <- random_walk[[s - 1]] + rnorm(1, mean = drift * step_length, sd = Varsigma * sqrt(step_length))
      if (random_walk[[s]] >= Alpha) {
        random_walk[s:nsteps] <- Alpha
        rts[[n]] <- s * step_length + ndt
        choice[[n]] <- 1 # Correct choice shown with positive RTs
        break
      } else if (random_walk[[s]] <= -Alpha) {
        random_walk[s:nsteps] <- -Alpha
        rts[[n]] <- s * step_length + ndt
        choice[[n]] <- -1 # Incorrect choice shown with positive RTs
        break
      } else if (s == (nsteps - 1)) {
        rts[[n]] <- NaN
        choice[[n]] <- NaN
        break
      }
    }
  }
  result <- rts * choice
  return(result)
}
simul_ratcliff_slow() # works?

# Simulate diffusion models quickly with intrinsic trial-to-trial variability in parameters
simulratcliff <- function(N = 100, Alpha = 1, Tau = .4, Nu = 1, Beta = .5, rangeTau = 0, rangeBeta = 0, Eta = .3, Varsigma = 1) {
  # SIMULRATCLIFF  Generates data according to a drift diffusion model with optional trial-to-trial variability
  #
  #
  # Reference:
  # Tuerlinckx, F., Maris, E.,
  # Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods for
  # simulating the diffusion process. Behavior Research Methods,
  # Instruments, & Computers, 33, 443-456.
  #
  # Parameters
  # ----------
  #   N: a integer denoting the size of the output vector
  # (defaults to 100 experimental trials)
  #
  # Alpha: the mean boundary separation across trials  in evidence units
  # (defaults to 1 evidence unit)
  #
  # Tau: the mean non-decision time across trials in seconds
  # (defaults to .4 seconds)
  #
  # Nu: the mean drift rate across trials in evidence units per second
  # (defaults to 1 evidence units per second, restricted to -5 to 5 units)
  #
  # Beta: the initial bias in the evidence process for choice A as a proportion of boundary Alpha
  # (defaults to .5 or 50% of total evidence units given by Alpha)
  #
  # rangeTau: Non-decision time across trials is generated from a uniform
  # distribution of Tau - rangeTau/2 to  Tau + rangeTau/2 across trials
  # (defaults to 0 seconds)
  #
  # rangeZeta: Bias across trials is generated from a uniform distribution
  # of Zeta - rangeZeta/2 to Zeta + rangeZeta/2 across trials
  # (defaults to 0 evidence units)
  #
  # Eta: Standard deviation of the drift rate across trials
  # (defaults to 3 evidence units per second, restricted to less than 3 evidence units)
  #
  # Varsigma: The diffusion coefficient, the standard deviation of the
  # evidence accumulation process within one trial. It is recommended that
  # this parameter be kept fixed unless you have reason to explore this parameter
  # (defaults to 1 evidence unit per second)
  #
  # Returns
  # -------
  # Vector with reaction times (in seconds) multiplied by the response vector
  # such that negative reaction times encode response B and positive reaction times
  # encode response A
  #
  #
  # Converted from simuldiff.m MATLAB script by Joachim Vandekerckhove,
  # Then converted from pyhddmjags utils python  script by Kianté Fernandez
  # See also http://ppw.kuleuven.be/okp/dmatoolbox.

  if (Nu < -5 || Nu > 5) {
    Nu <- sign(Nu) * 5
    warning(paste0("Nu is not in the range [-5 5], bounding drift rate to", Nu))
  }
  if (Eta > 3) {
    warning(paste0("Standard deviation of drift rate is out of bounds, bounding drift rate to 3"))
    Eta <- 3
  }
  if (Eta == 0) {
    Eta <- 1e-16
  }

  # Initialize output vectors
  results <- rep(0, N)
  Ts <- rep(0, N)
  XX <- rep(0, N)

  # Called sigma in 2001 paper
  D <- (Varsigma^(2)) / 2

  # Program specifications
  eps <- 2.220446049250313e-16 # precision from 1.0 to next double-precision number
  delta <- eps

  for (n in seq_len(N)) {
    r1 <- rnorm(1)
    mu <- Nu + r1 * Eta
    bb <- Beta - rangeBeta / 2 + rangeBeta * runif(1)
    zz <- bb * Alpha
    finish <- 0
    totaltime <- 0
    startpos <- 0
    Aupper <- Alpha - zz
    Alower <- -zz
    radius <- min(c(abs(Aupper), abs(Alower)))
    while (finish == 0) {
      lambda_ <- 0.25 * mu^(2) / D + 0.25 * D * pi^(2) / radius^(2)
      # eq. formula (13) in 2001 paper with D = sigma^2/2 and radius = Alpha/2
      Fs <- D * pi / (radius * mu)
      Fs <- Fs^(2) / (1 + Fs^(2))
      # formula p447 in 2001 paper
      prob <- exp(radius * mu / D)
      prob <- prob / (1 + prob)
      dir_ <- 2 * (runif(1) < prob) - 1
      l <- -1
      s2 <- 0
      while (s2 > l) {
        s2 <- runif(1)
        s1 <- runif(1)
        tnew <- 0
        told <- 0
        uu <- 0
        while (abs(tnew - told) > eps || uu == 0) {
          told <- tnew
          uu <- uu + 1
          tnew <- told + (2 * uu + 1) * -1^(uu) * s1^(Fs * 2 * uu + 1^(2))
          # infinite sum in formula (16) in BRMIC,2001
        }
        l <- 1 + s1^(-Fs) * tnew
      }
      # rest of formula (16)
      t <- abs(log(s1)) / lambda_
      # is the negative of t* in (14) in BRMIC,2001
      totaltime <- totaltime + t
      dir_ <- startpos + dir_ * radius
      ndt <- Tau - rangeTau / 2 + rangeTau * runif(1)
      if ((dir_ + delta) > Aupper) {
        Ts[n] <- ndt + totaltime
        XX[n] <- 1
        finish <- 1
      } else if ((dir_ - delta) < Alower) {
        Ts[n] <- ndt + totaltime
        XX[n] <- -1
        finish <- 1
      } else {
        startpos <- dir_
        radius <- min(abs(c(Aupper - startpos, Alower - startpos)))
      }
    }
  }
  result <- Ts * XX
  return(result)
}

simulratcliff() # works?

simuldiff2ndt <- function(N = 100, Alpha = 1, Vet = .2, Rmr = .2, Nu = 1, Zeta = None, rangeVet = 0, rangeRmr = 0, rangeZeta = 0, Eta = .3, Varsigma = 1) {
  #   SIMULDIFF2NDT  Generates data according to a diffusion model with non-decision time split into two parts with independent variance
  #
  #
  #   Reference:
  #   Tuerlinckx, F., Maris, E.,
  #   Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods for
  #   simulating the diffusion process. Behavior Research Methods,
  #   Instruments, & Computers, 33, 443-456.
  #
  #   Parameters
  #   ----------
  #   N: a integer denoting the size of the output vector
  #   (defaults to 100 experimental trials)
  #
  #   Alpha: the mean boundary separation across trials  in evidence units
  #   (defaults to 1 evidence unit)
  #
  #   Vet: the mean visual encoding time across trials in seconds
  #   (defaults to .2 seconds)
  #
  #   Rmr: the mean residual motor response time across trials in seconds
  #   (defaults to .2 seconds)
  #
  #   Nu: the mean drift rate across trials in evidence units per second
  #   (defaults to 1 evidence units per second, restricted to -5 to 5 units)
  #
  #   Zeta: the initial bias in the evidence process for choice A
  #   (defaults to 50% of total evidence units given by Alpha)
  #
  #   rangeVet: Visual encoding time across trials is generated from a uniform
  #   distribution of Vet - rangeVet/2 to  Vet + rangeVet/2 across trials
  #   (defaults to 0 seconds)
  #
  #   rangeRmr: Residual motor response time across trials is generated from a uniform
  #   distribution of Rmr - rangeRmr/2 to  Rmr + rangeRmr/2 across trials
  #   (defaults to 0 seconds)
  #
  #   rangeZeta: Bias across trials is generated from a uniform distribution
  #   of Zeta - rangeZeta/2 to Zeta + rangeZeta/2 across trials
  #   (defaults to 0 evidence units)
  #
  #   Eta: Standard deviation of the drift rate across trials
  #   (defaults to 3 evidence units per second, restricted to less than 3 evidence units)
  #
  #   Varsigma: The diffusion coefficient, the standard deviation of the
  #   evidence accumulation process within one trial. It is recommended that
  #   this parameter be kept fixed unless you have reason to explore this parameter
  #   (defaults to 1 evidence unit per second)
  #
  #   Returns
  #   -------
  #   Numpy complex vector with 1) Real component ( np.real(x) ): reaction times (in seconds) multiplied by the response vector
  #   such that negative reaction times encode response B and positive reaction times
  #   encode response A  and 2) Imaginary component ( np.imag(x) ): N200 peak-latencies in seconds
  #
  #
  # Converted from simuldiff.m MATLAB script by Joachim Vandekerckhove.
  # Then, converted from pyhddmjagsutils.py Python script by Kianté Fernandez
  # See also http://ppw.kuleuven.be/okp/dmatoolbox.
}
