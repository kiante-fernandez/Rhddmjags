# simpleCPP_sim.R - Directly simulates a model that
# assumes CPP slopes are generated from drift rates
#
# This R code was generated using Michael D. Nunez's `simpleCPP_sim.py`

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
# Date              Programmers                          Descriptions of Change
# ====            ================                       ======================
# 17-June-2022     Kianté Fernandez                             Original code


# Libraries
library(tidyverse)
library(patchwork)

# Simulation data and parameters
ntrials <- 1000
alpha <- 1 # Boundary parameter
ndt <- .4 # Non-decision time in seconds
delta <- 1 # Mean drift-rate across trials
beta <- .5 # Relative start point, proportion of boundary
eta <- .3 # Additional trial-to-trial variability in drift rate
varsigma <- 1 # Accumulation variance (diffusion coefficient)
sigma <- .5 # Observation noise of the CPP slope
nsteps <- 300 # nsteps*step_length is seconds after ndt
step_length <- .01


# Initialize output vectors
rts <- rep(0, ntrials)
choice <- rep(0, ntrials)
cpp_slopes <- rep(0, ntrials)
CPPs <- matrix(0, nsteps, ntrials)

# Try different random seeds
set.seed(2022)

# Set up plotting
p <- ggplot() +
  theme_classic()
p_ccp <- ggplot() +
  theme_classic()
# Direct simulation
plot_time <- seq(0, step_length * nsteps, length.out = nsteps) # Time vector
for (n in seq_len(ntrials)) {
  random_walk <- vector(mode = "numeric", length = nsteps)
  drift <- rnorm(n = 1, mean = delta, sd = eta)
  cpp_slopes[[n]] <- rnorm(n = 1, mean = delta, sd = sigma)
  CPPs[, n] <- sin(2 * pi * ((cpp_slopes[[n]] / 4) * (plot_time - ndt)))
  random_walk[[1]] <- beta * alpha
  for (s in 2:nsteps) {
    random_walk[[s]] <- random_walk[[s - 1]] + rnorm(1, mean = drift * step_length, sd = varsigma * sqrt(step_length))
    if (random_walk[[s]] >= alpha) {
      random_walk[s:nsteps] <- alpha
      rts[[n]] <- s * step_length + ndt
      choice[[n]] <- 1 # Correct choice shown with positive RTs
      break
    } else if (random_walk[[s]] <= 0) {
      random_walk[s:nsteps] <- 0
      rts[[n]] <- s * step_length + ndt
      choice[[n]] <- -1 # Incorrect choice shown with positive RTs
      break
    } else if (s == (nsteps - 1)) {
      rts[[n]] <- NaN
      choice[[n]] <- NaN
      break
    }
  }
  if (n < 50) {
    plot_data1 <- data.frame(cbind(cat = as.integer(n), random_walk, plot_time = plot_time + ndt))
    p <- p + geom_line(aes(x = plot_time, y = random_walk, color = factor(cat)),
      data = plot_data1
    )
    plot_data2 <- data.frame(cbind(cat = as.integer(n), CPPs = CPPs[, n], plot_time = plot_time))
    p_ccp <- p_ccp + geom_line(aes(x = plot_time, y = CPPs, color = factor(cat)),
      data = plot_data2
    )
  }
}

# Plotting parameters for evidence accumulation plots
p1 <- p + theme(legend.position = "none") +
  labs(
    title = "Simulated evidence accumulation paths",
    x = "Time (sec)",
    y = "Evidence for \u03bcV"
  ) +
  geom_hline(yintercept = 1, size = 1) +
  geom_hline(yintercept = 0.5, size = .5, linetype = "dashed") +
  geom_hline(yintercept = 0, size = 1) +
  scale_x_continuous(limits = c(ndt, 2))

p3 <- p_ccp + theme(legend.position = "none") +
  labs(
    title = "Simulated CPPs on single-trials",
    x = "Time (sec)",
    y = "Evidence for correct \u03bcV"
  ) +
  scale_x_continuous(limits = c(ndt, 2)) +
  scale_y_continuous(limits = c(0, 1))

# Plot estimated density of response times
p2 <- ggplot(data_frame(choice_rts = rts * choice), aes(choice_rts)) +
  geom_density(color = "blue", fill = "lightskyblue") +
  theme_classic() +
  labs(
    title = "Estimated Density of Incorrect and Correct Response Times",
    x = "Correct (+) and Incorrect (-) Response Times (secs)",
    y = "Density"
  )
# Plot estimated density of CPP slopes
p4 <- ggplot(data_frame(cpp_slopes), aes(cpp_slopes)) +
  geom_density(color = "forestgreen", fill = "springgreen3") +
  theme_classic() +
  labs(
    title = "Approximated Density of Simulated CPP Slopes",
    x = "CPP slopes (\u03bcV)",
    y = "Density"
  )

#print out the subplots
(p1 + p2) / (p3 + p4) # Combine plots

# Statistics of simulation
print(paste0("There are ", sum(is.na(rts)), " missing responses"))
print(paste0("The mean response time is ", mean(rts, na.rm = T), " seconds"))
print(paste0("The minimum response time is ", min(rts, na.rm = T), " seconds"))
print(paste0("The maximum response time is ", max(rts, na.rm = T), " seconds"))


# Plot estimate CDF of response times
data_frame(rts = rts, choice = choice) %>%
  mutate(correct = factor(ifelse(choice == 1, 1, 0), labels = c(
    "Incorrect",
    "Correct"
  ))) %>%
  ggplot(aes(rts, color = correct, group = correct)) +
  stat_function(fun = pnorm, position = position_dodge(0.05)) + #move postision
  theme_classic() +
  labs(
    title = "Cumulative Density Function of Response Times",
    x = "Response time (secs)",
    y = "Cumulative probability",
    color = ""
  ) +
  scale_color_brewer(palette = "Set1")
