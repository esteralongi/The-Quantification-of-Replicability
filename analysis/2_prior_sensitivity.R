## 3_prior_sensitivity.R

## Prior sensitivity of the replication metrics to the prior on the effect
## heterogeneity tau_beta. The data are held fixed while the prior on tau_beta is
## shifted across percentiles of its elicited prior; each setting is refit with
## adaptive convergence and the metrics are recomputed.
##
## Use target = "tau_alpha" or "tau_sig" for the other heterogeneities.

library(RepliBayes)

data <- read.csv("data/synthetic_data.csv")

## Arguments (paper setting S = 3, consensus_level = 2):
##   target = "tau_beta" : shift the prior on the effect heterogeneity
##                         (use "tau_alpha" or "tau_sig" for the others)
##   consensus_level = 2, min_corroborating = 1        : consensus level / conditional threshold for the metrics
##   percentiles         : prior medians at which the model is refit (default grid)
sens <- sensitivity_prior(data, target = "tau_beta", consensus_level = 2, min_corroborating = 1)

dir.create("results", showWarnings = FALSE)
saveRDS(sens, "results/sensitivity_tau_beta.rds")

cat("\n--- Per-setting diagnostics ---\n")
print(sens$grid)
cat("\n--- Metric x prior-median grid (consensus_level = 2 metrics) ---\n")
print(round(sens$p_grid, 3))
