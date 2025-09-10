# EpiSky

EpiSky is an R package that performs inference of time-varying reproduction number trajectories using epidemic data and genomic data jointly.
The inputs are a time series of the total number of cases (prevalence) and a dated phylogenetic tree.
The main output is a posterior distribution of the time-varying reproduction number.
EpiSky implements a [particle marginal Metropolis--Hastings algorithm](https://doi.org/10.1111/j.1467-9868.2009.00736.x) (PMMH) with the space of possible time series of latent prevalence and time-varying reproduction number is efficiently explored using sequential Monte Carlo (SMC).

For a methodological description of EpiSky, see the following pre-print:

A Gill, J Koskela, X Didelot, RG Everitt (2025), Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC, [arXiv:2311.09838v2](https://arxiv.org/abs/2311.09838v2).

## Installation

You can install EpiSky in R using the following command:

`devtools::install_github("alicia-gill/EpiSky")`

## Example

First, we simulate a dataset.

```r
library(EpiSky)

death_rate <- 0.1 #death rate of 0.1, so recovery takes 10 units of time
stop_time <- 30 #run for 30 units of time
x0 <- 1 #start with 1 index case

# Birth rate function - constant 0.3
constant <- function(t) {
  return(0.3)
}

set.seed(1)
epidemic <- epi_fn(birth_rate_fn = constant, death_rate = death_rate, stop_time = stop_time, x0 = x0)
```

Then extract the prevalence time series and the dated phylogeny from the simulated epidemic.

```r
prevalence <- prevalence(epidemic = epidemic, stop_time = stop_time)
phylogeny <- phylo_tree(epidemic = epidemic, stop_time = stop_time)

#plot the true prevalence and dated phylogeny if desired
# par(mfrow=c(1,2))
# par(mar=c(4,4,1,1))
# plot(prevalence, type = "o")
# library(ape)
# plot(phylogeny, show.tip.label = F)
# axisPhylo(backward = F)
```

To simulate the data, take a subsample of the prevalence and the phylogeny.

```r
rho <- 0.1 #10% reporting probability
set.seed(2)
prev_data <- sample_prevalence(prevalence = prevalence, reporting_prob = rho)

pi0 <- 0 #No chance of sampling a present day leaf
pi1 <- 0.1 #10% probability of sampling a historic leaf
set.seed(3)
phylo_data <- sample_phylo(ptree = phylogeny, pi0 = pi0, pi1 = pi1)

#plot the observed prevalence and dated phylogeny if desired
# par(mfrow=c(1,2))
# par(mar=c(4,4,1,1))
# plot(prev_data, type = "o")
# library(ape)
# plot(phylo_data$ptree, show.tip.label = F)
```

Choose priors for the reporting probability $\rho$, smoothness $\sigma$ and day 0 prevalence $X_{0}$.

```r
#prior on the reporting probability can be uniform or beta
pobs_prior <- "beta"
pobs_alpha <- 2
pobs_beta <- 5

#prior on the smoothness is exponential
sigma_mean <- 0.1

#prior on the time 0 prevalence can be uniform or negative binomial
x0_prior <- "nbinom"
x0_mean <- 5
x0_var <- 10
```

Set up for and run PMMH.

```r
iter <- 10000 #run for 10000 iterations
max_time <- 10*60 #run time limit of 10 minutes
target_acceptance <- 0.1 #target a 10% acceptance rate - close to optimal for PMMH
sigma0 <- 0.1 #initial mean smoothness of 0.1
reporting_prob0 <- 0.1 #initial reportion probability of 0.1
x0 <- 3 #initial time 0 prevalence of 3 cases
n_particles <- 1000 #it will be caculated automatically if left null
ess_threshold <- 500 #resampling will only happen in SMC if the ESS falls below this threshold
resampling_scheme <- "systematic" #use systematic resampling in SMC
backward_sim <- T #use backward simulation for sampling trajectories in SMC

set.seed(4)
chain <- pmmh(
  iter = iter,
  max_time = max_time,
  target_acceptance = target_acceptance,
  sigma0 = sigma0,
  reporting_prob0 = reporting_prob0,
  x0 = x0,
  death_rate = death_rate,
  ptree = phylo_data$ptree,
  ptree_lag = phylo_data$ptree_lag,
  sample_prevalence = prev_data,
  sigma_mean = sigma_mean,
  pobs_prior = pobs_prior,
  pobs_alpha = pobs_alpha,
  pobs_beta = pobs_beta,
  x0_prior = x0_prior,
  x0_mean = x0_mean,
  x0_var = x0_var,
  n_particles = n_particles,
  ess_threshold = ess_threshold,
  resampling_scheme = resampling_scheme,
  backward_sim = backward_sim
)
```

Look at the trace plots to assess convergence.

```r
par(mfrow=c(1,3))
plot(chain$reporting_prob, type="l", xlab="Iteration", ylab="Reporting probability")
plot(chain$sigma, type="l", xlab="Iteration", ylab="Smoothness")
plot(chain$x0, type="l", xlab="Iteration", ylab="Time 0 prevalence")
```

Finally, look at the inference plots.

```r
burn_in <- 2000
rt_mean <- apply(chain$birth_rate[-(1:burn_in),], 2, mean)/death_rate
rt_lcl <- matrixStats::colQuantiles(chain$birth_rate[-(1:burn_in),], probs=0.025)/death_rate
rt_ucl <- matrixStats::colQuantiles(chain$birth_rate[-(1:burn_in),], probs=0.975)/death_rate

par(mfrow=c(1,1))
plot(rt_mean, type="l", ylim=range(c(rt_lcl, rt_ucl)), xlab="Time", ylab="Reproduction number")
lines(rt_lcl, lty=2)
lines(rt_ucl, lty=2)
abline(h=3, col="red")
```

## Help

If you need help using EpiSky, you can email me at `draliciagill@gmail.com`.

