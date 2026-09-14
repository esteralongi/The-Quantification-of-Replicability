# The Quantification of Replicability

Code and data to reproduce the analyses and figures in **Alongi, Altoè &
Parmigiani (2026), *The Quantification of Replicability***.

The paper introduces a Bayesian hierarchical framework that quantifies the
replicability of an effect across several studies as posterior probabilities,
defined relative to a practical-relevance threshold `eps` on the scale of the
effects (rather than on p-values). The method is implemented in the companion R
package **[RepliBayes](https://github.com/esteralongi/RepliBayes)**; this
repository contains the scripts that apply it and produce the paper's figures and
tables.

## Data availability

The paper analyses donor-level **GTEx** data — the eQTL effect of SNP
rs4731702 on *KLF14* across three ancestral groups — which are
**controlled-access and cannot be redistributed**. To keep the pipeline fully
runnable, every script here uses a **synthetic dataset**
(`data/synthetic_data.csv`) with the same structure (`study`, `x`, `m`). The
results produced from it are therefore *illustrative* and will not reproduce the
exact numbers reported in the paper.

## Repository structure

```
The-Quantification-Of-Replicability/
├── analysis/
│   ├── 1_fit_replicability.R   # empirical metrics: hierarchical, independence, retrospective, prospective
│   ├── 2_simulation.R          # simulation calibration, varying each heterogeneity (tau_beta/alpha/sigma)
│   ├── 3_prior_sensitivity.R   # prior sensitivity to the prior on tau_beta
│   ├── 4_threshold_free.R      # threshold-free summaries of the study-specific effects
│   └── 5_figures.R             # paper figures (prior -> posterior movement; simulation calibration)
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
required to compile the two Stan models on first use. (Replace `esteralongi` with
your GitHub username once the package repo is created.)

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
source("analysis/1_fit_replicability.R")   # fits both models; writes results/fit_replicability.rds
source("analysis/2_simulation.R")          # simulation calibration (slow: R refits per scenario)
source("analysis/3_prior_sensitivity.R")   # refits across a grid of priors on tau_beta
source("analysis/4_threshold_free.R")      # threshold-free summaries
source("analysis/5_figures.R")             # writes the paper figures to figures/
```

Step 1 fits the hierarchical and independence-limit models and reports all
replication probabilities (each with a Monte Carlo standard error). Step 2
recomputes the simulation calibration for the three heterogeneity components.
Step 3 reruns the analysis while shifting the prior on the effect heterogeneity.
Step 5 rebuilds the figures from the saved results.

## Outputs

| script | output |
|---|---|
| `1_fit_replicability.R` | `results/fit_replicability.rds`; empirical, retrospective and prospective metric tables |
| `2_simulation.R` | `results/simulation.rds`; per-scenario metrics with MCSE |
| `3_prior_sensitivity.R` | `results/sensitivity_tau_a.rds`; metric-by-prior-median grid |
| `4_threshold_free.R` | max-divergence and 3-way density-overlap summaries |
| `5_figures.R` | `figures/simulation_calibration.pdf`, `figures/prior_posterior.pdf` |

## Citation

> Alongi, E., Altoè, G. & Parmigiani, G. (2026). *The Quantification of
> Replicability.* Statistical Science (to appear).
