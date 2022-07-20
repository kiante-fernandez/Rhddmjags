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
# 13/07/22      Kianté Fernandez                      Jellyfish plot code
# 20/07/22      Kianté Fernandez                      Recovery plot

# Libraries

library(dplyr) # A Grammar of Data Manipulation, CRAN v1.0.8
library(tidyr) # Tidy Messy Data, CRAN v1.2.0
library(readr) # Read Rectangular Text Data, CRAN v2.1.1
library(magrittr) # A Forward-Pipe Operator for R, CRAN v2.0.2
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics, CRAN v3.3.5
library(here) # A Simpler Way to Find Your Files, CRAN v1.0.1
require(rjags) # Bayesian Graphical Models using MCMC, CRAN v4-12. NOTE: Must have previously installed package rjags.
library(ggstar)
library(coda)
library(gtools)


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
# simul_ratcliff_slow() # works?

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

# simulratcliff() # works?

diagnostic <- function(insamples) {
  insamples <- samps # for testing

  # Returns two versions of Rhat (measure of convergence, less is better with an approximate
  # 1.10 cutoff) and Neff, number of effective samples). Note that 'rhat' is more diagnostic than 'oldrhat' according to
  # Gelman et al. (2014).
  #
  # Reference for preferred Rhat calculation (split chains) and number of effective sample calculation:
  #     Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. & Rubin, D. B. (2014).
  #     Bayesian data analysis (Third Edition). CRC Press:
  #     Boca Raton, FL
  #
  # Reference for original Rhat calculation:
  #     Gelman, A., Carlin, J., Stern, H., & Rubin D., (2004).
  #     Bayesian Data Analysis (Second Edition). Chapman & Hall/CRC:
  #     Boca Raton, FL.
  #
  # Parameters
  # ----------
  # insamples: array
  # Sampled values of monitored variables as an array of size
  #  variable by dim_1, dim_n, iterations, chains
  # dim_1, ..., dim_n describe the shape of variable in JAGS model.
  #
  # Returns
  # -------
  # list:
  #     rhat, oldrhat, neff, posterior mean, and posterior std for each variable. Prints maximum Rhat and minimum Neff across all variables

  result <- list()

  maxrhatsold <- rep(0, length(attributes(insamples)$dimnames$var))
  maxrhatsnew <- rep(0, length(attributes(insamples)$dimnames$var))
  minneff <- rep(0, length(attributes(insamples)$dimnames$var))

  var_rec <- table(factor(gsub("\\[|\\]|[0-9]+", "", attributes(insamples)$dimnames$var)))
  keys <- names(var_rec)
  var_order <- factor(gsub("\\[|\\]|[0-9]+", "", attributes(insamples)$dimnames$var))

  var_order[var_order == keys[1]]

  for (key_idx in seq_len(length(keys))) {
    result[[key_idx]] <- list()

    possamps <- insamples[, which(var_order == keys[key_idx]), ]
    dim(possamps)
    # Number of chains
    nchains <- dim(possamps)[[3]]

    # Number of samples per chain
    nsamps <- dim(possamps)[[1]]

    # Number of variables per key
    nvars <- dim(possamps)[[2]]

    # Reshape data
    allsamps <- apply(possamps, 1, c)
    dim(allsamps)

    # Reshape data to preduce R_hatnew
    possampsnew <- array(0, c(nvars, nsamps / 2, nchains * 2))
    newc <- 1
    dim(possamps)
    dim(possampsnew)
    for (chain_idx in seq_len(nchains)) {
      # possampsnew[,,newc] = (((possamps,seq_len(nsamps/2),axis=-2),chain_idx
      # possampsnew[,,newc + ] = (((possamps,seq_len(nsamps/2),axis=-2),chain_idx)
      # this needs to be checked
      possampsnew[, , newc] <- possamps[seq_len(nsamps / 2), , 1]
      possampsnew[, , newc + 1] <- possamps[seq_len(nsamps / 2), , 1]
      newc <- newc + 1

      # Index of variables
      varindx <- which(var_order == keys[key_idx])

      # Reshape data (all aperm not working?)
      # alldata <-  aperm(possamps, c(nvars, nsamps, nchains))

      # Mean of each chain for rhat
      chainmeans <- apply(possamps, 3, mean)
      # Mean of each chain for rhatnew
      chainmeansnew <- apply(possampsnew, 3, mean)
      # Global mean of each parameter for rhat
      globalmean <- mean(chainmeans)
      globalmeannew <- mean(chainmeansnew)
      result[[1]]$mean <- globalmean
      result[[1]]$std <- apply(allsamps, 1, std)
    }
  }
}

jellyfish <- function(samples, parameter, filename = NULL) {
  #   Plots posterior distributions of given posterior samples in a jellyfish
  #   plot. Jellyfish plots are posterior distributions (mirrored over their
  #   horizontal axes) with 99% and 95% credible intervals (currently plotted
  #   from the .5% and 99.5% & 2.5% and 97.5% percentiles respectively.
  #   Also plotted are the median, mode, and mean of the posterior distributions"
  #
  # Parameters
  # ----------
  # samples the rjags object
  # parameter: a string with the parameter of interest
  # optional file location to save the plot
  #
  # for calculating the highest density point (mode) per parameter
  highestdensity <- function(v) {
    temp_idx <- which.max(density(as.numeric(v))[["y"]])
    density(as.numeric(v))[["x"]][temp_idx]
  }
  ## if sample_dat  is the model output from R2jags
  sample_dat <- as.data.frame(as.matrix(as.mcmc(samples)))

  ## name your predicted factor latent.mean, and the CI between latent.lower and latent.upper
  post_mean <- apply(sample_dat, 2, mean)
  post_median <- apply(sample_dat, 2, median)
  post_mode <- apply(sample_dat, 2, highestdensity)
  # get the intervals
  post_lower1 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.025)))
  post_upper1 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.975)))
  post_lower2 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.005)))
  post_upper2 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.995)))
  # get parameter names
  subject <- colnames(sample_dat)
  # create data frame of statistics to plot
  dat <- data.frame(post_mean, post_median, post_mode, post_lower1, post_upper1, post_lower2, post_upper2, subject)
  # remove white space and symbols for indexing
  var_order <- factor(gsub("[[:punct:]]", "", gsub("\\[|\\]|[0-9]+", "", rownames(dat))))
  # select the parameters of interest to plot:
  plt_data <- dat[which(var_order == parameter), ]
  # order the data by the posterior MAPS, for better visualization
  plt_data <- plt_data[order(plt_data$post_mean), ]
  # order the observation IDs
  plt_data$subject2 <- reorder(plt_data$subject, plt_data$post_mean)
  # get the title of the plot
  title <- paste0("Posterior distributions of ", parameter, " parameter")
  # make plot using the ggplot:
  # orange square (mode)
  # black circle (median)
  # cyan star (mean)
  jellyplot <- ggplot(plt_data, aes(x = post_mean, y = subject2)) +
    geom_segment(aes(x = post_lower2, xend = post_upper2, y = subject2, yend = subject2), color = "cyan2", size = 1) +
    geom_segment(aes(x = post_lower1, xend = post_upper1, y = subject2, yend = subject2), color = "blue", size = 2) +
    geom_point(aes(x = post_mode), color = "darkorange", shape = 15, size = 3) +
    geom_point(aes(x = post_median), color = "black", shape = 16, size = 4) +
    ggstar::geom_star(color = "cyan2", fill = "cyan2", size = 3) +
    labs(title = title, x = "", y = "") +
    theme_classic()

  if (!is.null(filename)) {
    jellyplot
    ggsave(filename, dpi = 300)
  }
  return(jellyplot)
}

recovery <- function(samples, truevals, filename = NULL) {
  # Plots true parameters versus 99% and 95% credible intervals of recovered
  # parameters. Also plotted are the median (circles) and mean (stars) of the posterior
  # distributions.
  #
  # Parameters
  # ----------
  # samples :
  #
  # truevals :List of true parameter values (the genparam list)
  #
  # filename: optional
  
  # true value on the X axis?
  true_paramname <- names(truevals)
  
  ## if sample_dat  is the model output from R2jags
  sample_dat <- as.data.frame(as.matrix(as.mcmc(samples)))
  
  ## name your predicted factor and the CI between lower and upper
  post_mean <- apply(sample_dat, 2, mean)
  post_median <- apply(sample_dat, 2, median)
  
  # get the intervals
  post_lower1 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.025)))
  post_upper1 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.975)))
  post_lower2 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.005)))
  post_upper2 <- apply(sample_dat, 2, function(x) quantile(x, probs = c(0.995)))
  
  # get parameter names
  paramname <- colnames(sample_dat)
  
  # create data frame of statistics to plot
  dat <- data.frame(post_mean, post_median, post_lower1, post_upper1, post_lower2, post_upper2, paramname)
  # remove white space and symbols for indexing
  var_order <- factor(gsub("[[:punct:]]", "", gsub("\\[|\\]|[0-9]+", "", rownames(dat))))
  # select the parameters of interest to plot:
  plt_data <- dat[which(var_order == true_paramname), ]
  plt_data <- plt_data[mixedsort(sort(plt_data$paramname)), ] # sort data in "correct order"
  plt_data$truevals <- as.vector(t(truevals[[1]]))
  
  ## order the data by the posterior MAPS, for better visualization
  plt_data <- plt_data[order(plt_data$post_mean), ]
  
  ## order the observation IDs
  plt_data$paramname2 <- reorder(plt_data$paramname, plt_data$post_mean)
  ## get y = x for plotting
  plt_data$recoverline <- seq(min(plt_data$truevals), max(plt_data$truevals), length.out = nrow(plt_data))
  
  title <- paste0("Recovery of the ", true_paramname)
  
  ## make plot using the ggplot:
  recover_plot <- ggplot(plt_data, aes(x = truevals, y = post_mean)) +
    geom_segment(aes(x = truevals, xend = truevals, y = post_lower2, yend = post_upper2), color = "cyan2", size = 1) +
    geom_segment(aes(x = truevals, xend = truevals, y = post_lower1, yend = post_upper1), color = "blue", size = 2) +
    geom_point(aes(y = post_median), color = "black", shape = 16, size = 4) +
    ggstar::geom_star(color = "cyan2", fill = "cyan2", size = 3) +
    labs(title = title, x = "", y = "") +
    theme_classic() +
    geom_line(aes(x = recoverline, y = recoverline), color = "darkorange", size = 2)
  
  if (!is.null(filename)) {
    recover_plot
    ggsave(filename, dpi = 300)
  }
  return(recover_plot)
}

rsquared_pred <- function(trueval, predval) {
  # RSQUARED_PRED  Calculates R^2_prediction for data and statistics derived from data
  divisor <- sum(is.infinite(trueval)) - 1
  # Mean squared error of prediction
  MSEP <- sum((trueval - predval)^2) / divisor
  # Variance estimate of the true values
  vartrue <- sum((trueval - mean(trueval, na.rm = T))^2) / divisor
  # R-squared definition
  rsquared <- 1 - (MSEP / vartrue)
  return(rsquared)
}
