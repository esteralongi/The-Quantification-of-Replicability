## 1_fit_replicability.R
## -----------------------------------------------------------------------------
## Empirical replication analysis: fit the hierarchical and independence-limit
## models and compute all replication probabilities (study-level, generative,
## retrospective and prospective), each with a Monte Carlo standard error.
##
## Runs on the SYNTHETIC dataset; the paper uses controlled-access GTEx data, so
## the numbers here are illustrative and will not match the paper exactly.
## -----------------------------------------------------------------------------

library(RepliBayes)
set.seed(42)

data <- read.csv("data/synthetic_data.csv")   # columns: study, x, m

## Elicited empirical-Bayes priors used in the paper (edit to supply your own).
priors <- default_priors()

## One call fits both models and computes every metric.
## eps defaults to 10% of the baseline (x = 0) mean outcome.
res <- fit_replicability(data, priors = priors)

dir.create("results", showWarnings = FALSE)
saveRDS(res, "results/fit_replicability.rds")

cat(sprintf("eps = %.4g\n", res$eps))
print(res)

cat("\n--- Hierarchical model ---\n")
print(res$metrics_hierarchical,  n = Inf)
cat("\n--- Independence limit ---\n")
print(res$metrics_independence,  n = Inf)
cat("\n--- Retrospective (left-out study) ---\n")
print(res$metrics_retrospective, n = Inf)
cat("\n--- Prospective (new study) ---\n")
print(res$metrics_prospective,   n = Inf)
