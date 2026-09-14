## 4_threshold_free.R
## -----------------------------------------------------------------------------
## Threshold-free summaries of the three study-specific effects (no eps needed):
## the posterior mean/sd of the maximum divergence across studies, and the 3-way
## density overlap (overlapping coefficient).
## -----------------------------------------------------------------------------

library(RepliBayes)

res  <- readRDS("results/fit_replicability.rds")
post <- rstan::extract(res$fit_hierarchical)

## list of the S study-specific effect draws
beta_list <- lapply(seq_len(res$S), function(s) post$beta[, s])
tf <- replication_measures_multi(beta_list)

cat("Maximum divergence  (mean +/- sd):",
    sprintf("%.3f +/- %.3f\n", tf$mean_max_divergence, tf$sd_max_divergence))
cat("S-way density overlap            :",
    sprintf("%.3f\n", tf$density_overlap))
