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

## Illustrative priors for the synthetic data (edit to supply your own).
## The paper's elicited empirical-Bayes priors (Table 3) are derived from the
## controlled-access GTEx pool and cannot be redistributed.
priors <- default_priors()

## One call fits both models and computes every metric. Arguments used here
## reproduce the paper's setting (S = 3 studies, consensus level consensus_level = 2):
##   priors  : the hyperparameters (illustrative defaults here)
##   consensus_level = 2   : replication requires at least 2 of the 3 studies to agree
##   min_corroborating = 1   : each study's conditional metric requires >= 1 other agreeing study
##   eps     : practical-relevance threshold; the default (NULL) sets it to
##             10% of the baseline (x = 0) mean outcome
##   retrospective / prospective default to TRUE; with 3 studies the retrospective
##   analysis leaves out study 1 and uses studies 2-3 as the body of evidence,
##   and the prospective analysis draws a new study from all 3.
## (S is read from data$study, so nothing else is needed for the 3-study case.)
res <- fit_replicability(data, priors = priors, consensus_level = 2, min_corroborating = 1)

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
