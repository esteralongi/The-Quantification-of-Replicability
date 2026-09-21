# The Quantification of Replicability

Code and data to reproduce the analyses and figures in **Alongi, Altoè &
Parmigiani (2026), *The Quantification of Replicability***.

This paper offers a unified conceptual framework for discussing why, 
when and how to quantify the replicability of an effect across several studies.
It proposes to assess it through a Bayesian hierarchical model, working 
with posterior probabilities of the effects, defined relative to a 
practical-relevance threshold `eps`.
The method is implemented in the companion R package 
**[RepliBayes](https://github.com/esteralongi/RepliBayes)**; this
repository contains the scripts that apply it and produce the paper's figures 
and tables.

## Data availability

The paper analyses donor-level GTEx data — the eQTL effect of SNP
rs4731702 on *KLF14* across three ancestral groups — which are 
controlled-access and cannot be redistributed. To keep the pipeline fully
runnable, every script here uses a synthetic dataset
(`data/synthetic_data.csv`) with the same structure (`study`, `x`, `m`). The
results produced from it are therefore *illustrative* and will not reproduce the
exact numbers reported in the paper.

## Repository structure

```
The-Quantification-Of-Replicability/
├── analysis/
│   ├── 0_motivating_example.R  # two-topics motivating example (standalone; Section "The Added Value of Replicability")
│   ├── 1_fit_replicability.R   # empirical metrics: hierarchical, independence, retrospective, prospective
│   ├── 2_prior_sensitivity.R   # prior sensitivity to the prior on tau_beta
│   ├── 3_simulation.R          # simulation calibration, varying each heterogeneity (tau_beta/alpha/sigma)
│   └── 4_figures.R             # paper figure simulation calibration
├── data/
│   └── synthetic_data.csv      # synthetic stand-in for the GTEx data
├── results/                    # script outputs (.rds); created on first run
├── figures/                    # generated figures (PDF)
├── session_info.txt            # package versions (regenerate; see the file)
└── quantifying-replicability.Rproj
```

## The RepliBayes package

All the modelling lives in the companion package. Install it first:

```r
# install.packages("remotes")
remotes::install_github("esteralongi/RepliBayes")
```

`RepliBayes` uses [rstan](https://mc-stan.org/rstan/); a working C++ toolchain is
required to compile the two Stan models on first use.

## Dependencies

```r
# the framework
remotes::install_github("esteralongi/RepliBayes")

# for the analysis scripts and figures
install.packages(c("ggplot2", "tidyr", "patchwork"))
```

All scripts were developed under R >= 4.1.

## How to reproduce

All paths are relative to the repository root (the RStudio project directory).
Run the scripts in order:

```r
source("analysis/0_motivating_example.R")  # standalone; writes the two-topics figure
source("analysis/1_fit_replicability.R")   # fits both the hierarchical and the independence-limit models and reports
all replication probabilities; writes results/fit_replicability.rds
source("analysis/2_simulation.R")          # simulation calibration for the three heterogeneity components
source("analysis/3_prior_sensitivity.R")   # refits across a grid of priors on the effect heterogeneity tau_beta
source("analysis/4_figures.R")             # writes the paper figures to figures/
```

## Outputs

| script | output | Where in the paper |
|---|---|---|
| `0_motivating_example.R` | `figures/fig_two_topics.pdf`, `results/two_topics_suffstats.rds` | Section *Introduction - The Added Value of Replicability*; Fig. 1, Table 1 |
| `1_fit_replicability.R` | `results/fit_replicability.rds`; empirical, retrospective and prospective metric tables | Section *Empirical Application - Results*; Fig. 5. Supplement A 0 Section *Monte Carlo Diagnostics for the Empirical Replication Probabilities*; Table S1, S2, S3 |
| `2_prior_sensitivity.R` | `results/sensitivity_tau_beta.rds`; metric-by-prior-median grid | Section Empirical Application - Results; Table 4 |
| `3_simulation.R` | `results/simulation.rds`; per-scenario metrics with MCSE | Section *Calibrating the Metrics through Simulation - Empirically Grounded Design*; Table 5 6 7 8 |
| `4_figures.R` | `figures/simulation_calibration.pdf` | Section *Calibrating the Metrics through Simulation - Results*; Fig. 6 |

> Alongi, E., Altoè, G. & Parmigiani, G. (2026). *The Quantification of
> Replicability.*
