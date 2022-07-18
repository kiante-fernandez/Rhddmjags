# Rhddmjags
#### (Repository version 0.1.0)

Repository for example Hierarchical Drift Diffusion Model (HDDM) code using JAGS in R

**This code is based on Michael D. Nunez's Repo [`pyhddmjags`](https://github.com/mdnunez/pyhddmjags)**

### Prerequisites

[R](https://www.r-project.org/)

[R Studio](https://www.rstudio.com/products/rstudio/download/)

[MCMC Sampling Program: JAGS](http://mcmc-jags.sourceforge.net/)

[Program: JAGS Wiener module](https://sourceforge.net/projects/jags-wiener/)

[R Packages: ]() can use the following to install the required packages:

```bash
install.packages(c("dplyr", "tidyr", "readr", "magrittr", "ggplot2","here","R2jags","ggstar","coda"))

```

### Downloading

The repository can be cloned with `git clone https://github.com/kiantefernandez/Rhddmjags.git`

The repository can also be downloaded via the Code -> _Download zip_ buttons above on this Github page.

### Getting started

At the moment each script can be run individually to simulate from hierarchical drift-diffusion models (HDDMs) and then find and recover parameter estimates from those models. The most simple HDDM lives in nolapse_test.R. See other scripts: recovery_test.R, blocked_exp_conds.R, and regression_test.R. These scripts provide useful examples for using JAGS with Rjags, the JAGS Wiener module, mixture modeling in JAGS, and Bayesian diagnostics in R. 

The script nolapse_test_Rstan.R contains Rstan and Stan code to find and recover parameters from the exact same HDDM written in JAGS within nolapse_test.R. 

### License

Rhddmjags is licensed under the GNU General Public License v3.0 and written by Kianté Fernandez from the Neruoeconomics group at the Ohio State Univeristy.

### Possible citation

### Usage examples
