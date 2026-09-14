## 5_figures.R
## -----------------------------------------------------------------------------
## Paper figures, rebuilt from the saved results:
##   figures/prior_posterior.pdf        prior -> posterior movement of the metrics
##   figures/simulation_calibration.pdf simulation calibration (3 heterogeneities)
##
## All metrics use the INTERNAL names of compute_replication_probs_hier(), which
## are the names returned by simulate_replicability().
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
metrics_k2 <- c("Agr_2", "P_Rep_2", "P_pos_2", "P_null_2", "P_neg_2",
                "Cond_S1", "Cond_S2", "Cond_S3", "P_beta",
                "Agr_Gen_2", "Agr_Gen_non_null_2", "Pos_Gen_2", "Null_Gen_2",
                "Neg_Gen_2", "CondPos_Gen", "CondNull_Gen", "CondNeg_Gen")

label_dict <- c(
  "Agr_2"              = "P['overall, 2']",
  "P_Rep_2"            = "P['non-null, 2']",
  "P_pos_2"            = "P['pos, 2']",
  "P_null_2"           = "P['null, 2']",
  "P_neg_2"            = "P['neg, 2']",
  "Cond_S1"            = "P[list('cond' ~ '|' ~ O == 1, m == 1)]",
  "Cond_S2"            = "P[list('cond' ~ '|' ~ O == 2, m == 1)]",
  "Cond_S3"            = "P[list('cond' ~ '|' ~ O == 3, m == 1)]",
  "P_beta"             = "P[beta]",
  "Agr_Gen_2"          = "P['gen, overall, 2']",
  "Agr_Gen_non_null_2" = "P['gen, non-null, 2']",
  "Pos_Gen_2"          = "P['gen, pos, 2']",
  "Null_Gen_2"         = "P['gen, null, 2']",
  "Neg_Gen_2"          = "P['gen, neg, 2']",
  "CondPos_Gen"        = "P['c, pos, 2']",
  "CondNull_Gen"       = "P['c, null, 2']",
  "CondNeg_Gen"        = "P['c, neg, 2']"
)

## helper: named list of probabilities -> tidy data frame (metric, value) ------
probs_to_df <- function(pr) {
  d <- data.frame(metric = names(pr), value = as.numeric(unlist(pr)),
                  stringsAsFactors = FALSE)
  d[d$metric %in% metrics_k2, ]
}

## posterior reference (from the empirical hierarchical fit) -------------------
post <- rstan::extract(res$fit_hierarchical)
probs_post <- compute_replication_probs_hier(post$a[, 1], post$a[, 2], post$a[, 3],
                                             post$mu_a, eps)

## prior reference (coherent draws: mu is the SAME beta that generates a1,a2,a3)
Nd         <- 50000
tau_a_draw <- rtrunc_t_pos(Nd, nu, priors$mu_tau_a, priors$scale_tau_a)
beta_draw  <- priors$mu_a + priors$scale_a * rt(Nd, df = nu)
a1 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
a2 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
a3 <- beta_draw + tau_a_draw * rt(Nd, df = nu)
probs_prior <- compute_replication_probs_hier(a1, a2, a3, beta_draw, eps)

## ============================================================================
## Figure 1: prior -> posterior movement
## ============================================================================
ref <- rbind(
  data.frame(type = "Prior",     probs_to_df(probs_prior)),
  data.frame(type = "Posterior", probs_to_df(probs_post))
)
ref$metric <- factor(ref$metric, levels = metrics_k2)
ref$type   <- factor(ref$type, levels = c("Prior", "Posterior"))

g1 <- ggplot(ref, aes(metric, value, group = metric)) +
  geom_line(color = "grey70", linewidth = 0.4) +
  geom_point(aes(shape = type, color = type), size = 2.3) +
  scale_shape_manual(values = c(Prior = 17, Posterior = 16)) +
  scale_color_manual(values = c(Prior = "grey55", Posterior = "black")) +
  scale_x_discrete(labels = function(x) parse(text = label_dict[x])) +
  ylim(0, 1) +
  labs(x = "", y = "Probability", shape = "", color = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", panel.grid.major.x = element_blank())

ggsave("figures/prior_posterior.pdf", g1, width = 9, height = 4.5)

## ============================================================================
## Figure 2: simulation calibration (three heterogeneity scenarios)
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

cat("Wrote figures/prior_posterior.pdf and figures/simulation_calibration.pdf\n")
