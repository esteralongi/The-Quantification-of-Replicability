## 4_figures.R
## -----------------------------------------------------------------------------
## Paper figure, rebuilt from the saved results:
##   figures/simulation_calibration.pdf  simulation calibration (3 heterogeneities),
##                                        with the empirical prior and posterior
##                                        overlaid as reference lines.
##
## All metrics use the paper-style names returned by compute_replication_probs_hier()
## and simulate_replicability() (P_overall_2, P_gen_pos_2, P_c_pos_2, ...).
## Run steps 1 and 2 first (they write results/fit_replicability.rds and
## results/simulation.rds).
## -----------------------------------------------------------------------------

library(RepliBayes)
library(ggplot2)
set.seed(42)

dir.create("figures", showWarnings = FALSE)

res    <- readRDS("results/fit_replicability.rds")
priors <- res$priors
eps    <- res$eps
nu     <- priors$nu

## metric order and pretty (plotmath) labels ----------------------------------
metrics_k2 <- c("P_overall_2", "P_non_null_2", "P_pos_2", "P_null_2", "P_neg_2",
                "P_cond_O1_m1", "P_cond_O2_m1", "P_cond_O3_m1", "P_beta",
                "P_gen_overall_2", "P_gen_non_null_2", "P_gen_pos_2", "P_gen_null_2",
                "P_gen_neg_2", "P_c_pos_2", "P_c_null_2", "P_c_neg_2")

label_dict <- c(
  "P_overall_2"      = "P['overall, 2']",
  "P_non_null_2"     = "P['non-null, 2']",
  "P_pos_2"          = "P['pos, 2']",
  "P_null_2"         = "P['null, 2']",
  "P_neg_2"          = "P['neg, 2']",
  "P_cond_O1_m1"     = "P[list('cond' ~ '|' ~ O == 1, m == 1)]",
  "P_cond_O2_m1"     = "P[list('cond' ~ '|' ~ O == 2, m == 1)]",
  "P_cond_O3_m1"     = "P[list('cond' ~ '|' ~ O == 3, m == 1)]",
  "P_beta"           = "P[beta]",
  "P_gen_overall_2"  = "P['gen, overall, 2']",
  "P_gen_non_null_2" = "P['gen, non-null, 2']",
  "P_gen_pos_2"      = "P['gen, pos, 2']",
  "P_gen_null_2"     = "P['gen, null, 2']",
  "P_gen_neg_2"      = "P['gen, neg, 2']",
  "P_c_pos_2"        = "P['c, pos, 2']",
  "P_c_null_2"       = "P['c, null, 2']",
  "P_c_neg_2"        = "P['c, neg, 2']"
)

## helper: named list of probabilities -> tidy data frame (metric, value) ------
probs_to_df <- function(pr) {
  d <- data.frame(metric = names(pr), value = as.numeric(unlist(pr)),
                  stringsAsFactors = FALSE)
  d[d$metric %in% metrics_k2, ]
}

## posterior reference (from the empirical hierarchical fit) -------------------
post <- rstan::extract(res$fit_hierarchical)
probs_post <- compute_replication_probs_hier(list(post$beta[, 1], post$beta[, 2], post$beta[, 3]),
                                             post$mu_beta, eps)

## prior reference (coherent draws: mu is the SAME beta that generates a1,a2,a3)
Nd         <- 50000
tau_a_draw <- rtrunc_t_pos(Nd, nu, priors$mu_tau_beta, priors$scale_tau_beta)
beta_draw  <- priors$mu_beta + priors$scale_beta * rt(Nd, df = nu)
a1 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
a2 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
a3 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
probs_prior <- compute_replication_probs_hier(list(a1, a2, a3), beta_draw, eps)

## ============================================================================
## Simulation calibration (three heterogeneity scenarios), with the empirical
## prior and posterior overlaid as reference lines (the paper figure)
## ============================================================================
sim <- readRDS("results/simulation.rds")
lvlmap <- c(low = "Low heterogeneity", medium = "Medium heterogeneity",
            high = "High heterogeneity")

prep <- function(df, scenario) {
  df <- df[df$metric %in% metrics_k2, c("level", "metric", "value", "mcse")]
  df$level    <- lvlmap[as.character(df$level)]
  df$scenario <- scenario
  df
}
simdf <- rbind(
  prep(sim$beta,  "A: Effect heterogeneity"),
  prep(sim$alpha, "B: Intercept heterogeneity"),
  prep(sim$sigma, "C: Residual-scale heterogeneity")
)

## prior / posterior reference lines, replicated in every panel
ref2  <- rbind(
  data.frame(level = "Prior",     probs_to_df(probs_prior), mcse = 0),
  data.frame(level = "Posterior", probs_to_df(probs_post),  mcse = 0)
)
ref_all <- do.call(rbind, lapply(unique(simdf$scenario),
                                 function(sc) data.frame(ref2, scenario = sc)))

plotdf <- rbind(simdf, ref_all[, names(simdf)])
lev_order <- c("Low heterogeneity", "Medium heterogeneity", "High heterogeneity",
               "Prior", "Posterior")
plotdf$level  <- factor(plotdf$level, levels = lev_order)
plotdf$metric <- factor(plotdf$metric, levels = metrics_k2)

pal <- c("Low heterogeneity" = "#9ECAE1", "Medium heterogeneity" = "#4292C6",
         "High heterogeneity" = "#084594", "Prior" = "grey60", "Posterior" = "black")

g2 <- ggplot(plotdf, aes(metric, value, color = level, group = level)) +
  geom_ribbon(aes(ymin = pmax(0, value - 2 * mcse),
                  ymax = pmin(1, value + 2 * mcse), fill = level),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = function(x) parse(text = label_dict[x])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "", y = "Probability", color = "") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", panel.grid.major.x = element_blank())

ggsave("figures/simulation_calibration.pdf", g2, width = 7, height = 9)

cat("Wrote figures/simulation_calibration.pdf\n")
