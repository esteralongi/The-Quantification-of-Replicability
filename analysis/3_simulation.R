## 2_simulation.R

## Simulation calibration of the replication metrics. One component of the
## generative model is varied at a time (effect / intercept / residual-scale
## heterogeneity) across low / medium / high levels; for each level R datasets of
## three studies are simulated, refit, and the metrics are summarised across
## replicates with a Monte Carlo standard error.
##
## Uses the empirical fit from step 1 to set the ground-truth locations.
## Slow step: R refits per level, per scenario.

library(RepliBayes)

data <- read.csv("data/synthetic_data.csv")
res  <- readRDS("results/fit_replicability.rds")

R <- 30   # replicate datasets per level

## Arguments (paper setting S = 3, consensus_level = 2):
##   res    : the fitted object from step 1 (supplies the ground-truth locations)
##   vary   : which generative component varies across low/medium/high levels
##   consensus_level = 2  : consensus level for the metrics; min_corroborating = 1 for the conditional metrics
##   R      : number of simulated datasets per level
sim_beta  <- simulate_replicability(res, data, vary = "tau_beta",  consensus_level = 2, min_corroborating = 1, R = R)  # effect heterogeneity
sim_alpha <- simulate_replicability(res, data, vary = "tau_alpha", consensus_level = 2, min_corroborating = 1, R = R)  # intercept heterogeneity
sim_sigma <- simulate_replicability(res, data, vary = "tau_sig",   consensus_level = 2, min_corroborating = 1, R = R)  # residual-scale heterogeneity

dir.create("results", showWarnings = FALSE)
saveRDS(list(beta = sim_beta, alpha = sim_alpha, sigma = sim_sigma),
        "results/simulation.rds")

cat("\n--- Effect heterogeneity (tau_beta) ---\n")
print(sim_beta, n = Inf)
