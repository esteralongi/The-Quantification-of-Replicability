## 2_simulation.R
## -----------------------------------------------------------------------------
## Simulation calibration of the replication metrics. One component of the
## generative model is varied at a time (effect / intercept / residual-scale
## heterogeneity) across low / medium / high levels; for each level R datasets of
## three studies are simulated, refit, and the metrics are summarised across
## replicates with a Monte Carlo standard error.
##
## Uses the empirical fit from step 1 to set the ground-truth locations.
## This is the slow step: R refits per level, per scenario.
## -----------------------------------------------------------------------------

library(RepliBayes)

data <- read.csv("data/synthetic_data.csv")
res  <- readRDS("results/fit_replicability.rds")

R <- 30   # replicate datasets per level

sim_beta  <- simulate_replicability(res, data, vary = "tau_a",     R = R)  # effect heterogeneity
sim_alpha <- simulate_replicability(res, data, vary = "tau_alpha", R = R)  # intercept heterogeneity
sim_sigma <- simulate_replicability(res, data, vary = "tau_sig",   R = R)  # residual-scale heterogeneity

dir.create("results", showWarnings = FALSE)
saveRDS(list(beta = sim_beta, alpha = sim_alpha, sigma = sim_sigma),
        "results/simulation.rds")

cat("\n--- Effect heterogeneity (tau_beta) ---\n")
print(sim_beta, n = Inf)
