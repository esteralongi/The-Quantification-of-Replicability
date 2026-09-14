## motivating_example.R

## Motivating example (Section "The Added Value of Replicability").
## Two topics, two studies each (i = topic, j = study):
##
## study data:    ybar_ij | beta_ij        ~ N(beta_ij, sigma2_ij / n)
## study effects: beta_ij | beta_i, tau_i  ~ N(beta_i, tau_i^2)
## topic effects: beta_i                   ~ N(m0, A2)
##
## The between-study sd tau_i is unknown with prior LogNormal(log 0.5, 0.50^2),
## and the residual variances sigma2_ij are unknown with InvGamma(3, 2) priors
## (mean 1, sd 1); the sample variances s2_ij enter the sufficient statistics.
## Conditional on the variance parameters the model is conjugate normal-normal,
## so those parameters are integrated on a 3-D grid over (tau, sigma2_1,
## sigma2_2) and every beta posterior is a finite mixture of Gaussians

## All probabilities are P(beta > eps | y) with eps = 0.1
##
## Construction (as in the paper): ybar_11 = 1; ybar_12 solved so that
## P(beta_12 > eps | y) = 0.5 while P(beta_11 > eps | y) > 0.95. Topic 1 has a
## precise study (s2_11 = 0.5) and a noisy one (s2_12 = 4.0). Topic 2's common
## center and sample-variance level are solved so that E[beta_2|y] = E[beta_1|y]
## and SD[beta_2|y] = SD[beta_1|y], with the two study means straddling the
## center by +/- delta and their variances asymmetric
##
## This script is independent of the eQTL analysis: it uses its own eps = 0.1
## and n = 25, so run it in a fresh session


## fixed quantities 
n     <- 25       # sample size of every study
m0    <- 0        # prior mean of beta_i
A2    <- 1        # prior variance of beta_i
eps   <- 0.1      # minimally important effect: probabilities are P(beta > eps | y)
delta <- 0.20     # half-distance between the topic-2 study means

## mixture helpers (a posterior is a weighted mixture of N(m, V)) 
mix_mean <- function(wg, m)    sum(wg * m)
mix_sd   <- function(wg, m, V) sqrt(sum(wg * (V + m^2)) - mix_mean(wg, m)^2)
mix_prob <- function(wg, m, V) sum(wg * pnorm(eps, m, sqrt(V), lower.tail = FALSE))
mix_dens <- function(x, wg, m, V)
  vapply(x, function(xx) sum(wg * dnorm(xx, m, sqrt(V))), 0)

## variance-integration machinery
mlog <- log(0.5); slog <- 0.50          # tau_i ~ LogNormal(mlog, slog^2)
a0 <- 3; b0 <- 2                        # sigma2_ij ~ InvGamma(a0, b0): mean 1, sd 1
tau_grid <- seq(0.05, 3.0, length.out = 60)
ldinvgamma <- function(x, a, b) a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x
## log p(s2 | sigma2): (n-1) s2 / sigma2 ~ chisq(n-1)
lds2 <- function(s2, sigma2) dchisq((n - 1) * s2 / sigma2, n - 1, log = TRUE) +
  log(n - 1) - log(sigma2)
## per-study sigma2 grid from its (conjugate, s2-only) posterior quantiles
qig <- function(p, a, b) 1 / qgamma(1 - p, a, rate = b)
sig_grid <- function(s2) {
  a <- a0 + (n - 1) / 2; b <- b0 + (n - 1) * s2 / 2
  seq(qig(5e-4, a, b), qig(1 - 5e-4, a, b), length.out = 40)
}

## per topic: 3-D grid posterior over (tau, sigma2_1, sigma2_2) with conditional
## posterior moments; returns mixture components for beta_i (Eb, Vb), beta_ij
## (Em, Vm: N x 2), and sigma2_ij (sg: N x 2)
topic_post <- function(ybar, s2) {
  g  <- expand.grid(tau = tau_grid, sg1 = sig_grid(s2[1]), sg2 = sig_grid(s2[2]))
  v_1 <- g$sg1 / n; v_2 <- g$sg2 / n
  V1 <- A2 + g$tau^2 + v_1; V2 <- A2 + g$tau^2 + v_2; C <- A2
  d  <- V1 * V2 - C^2
  y1 <- ybar[1] - m0; y2 <- ybar[2] - m0
  lml <- -log(2 * pi) - 0.5 * log(d) -
    0.5 * (V2 * y1^2 - 2 * C * y1 * y2 + V1 * y2^2) / d
  lw <- lml + dlnorm(g$tau, mlog, slog, log = TRUE) +
    ldinvgamma(g$sg1, a0, b0) + lds2(s2[1], g$sg1) +
    ldinvgamma(g$sg2, a0, b0) + lds2(s2[2], g$sg2)
  wg <- exp(lw - max(lw)); wg <- wg / sum(wg)
  p1 <- 1 / (v_1 + g$tau^2); p2 <- 1 / (v_2 + g$tau^2)
  P  <- 1 / A2 + p1 + p2
  Eb <- (m0 / A2 + p1 * ybar[1] + p2 * ybar[2]) / P;  Vb <- 1 / P
  w1 <- (n / g$sg1) / (n / g$sg1 + g$tau^-2)
  w2 <- (n / g$sg2) / (n / g$sg2 + g$tau^-2)
  Em <- cbind(w1 * ybar[1] + (1 - w1) * Eb, w2 * ybar[2] + (1 - w2) * Eb)
  Vm <- cbind(1 / (n / g$sg1 + g$tau^-2) + (1 - w1)^2 * Vb,
              1 / (n / g$sg2 + g$tau^-2) + (1 - w2)^2 * Vb)
  list(wg = wg, tau = g$tau, Eb = Eb, Vb = Vb, Em = Em, Vm = Vm,
       sg = cbind(g$sg1, g$sg2))
}

## construct the sufficient statistics of the example 
ybar11  <- 1
s2_1    <- c(0.5, 4.0)       # topic-1 sample variances: precise vs noisy study
s2ratio <- c(1.4, 0.7)       # topic-2 variance asymmetry

## topic-1 second study mean: solved so P(beta_12 > eps | y) = 0.5
ybar12 <- uniroot(function(y12) {
  tp <- topic_post(c(ybar11, y12), s2_1)
  mix_prob(tp$wg, tp$Em[, 2], tp$Vm[, 2]) - 0.5
}, c(-1.0, 0.4), tol = 1e-8)$root
tp1 <- topic_post(c(ybar11, ybar12), s2_1)
E1  <- mix_mean(tp1$wg, tp1$Eb)
S1  <- mix_sd(tp1$wg, tp1$Eb, tp1$Vb)

## topic-2 center s and sample-variance level c: solved so that
## E[beta_2 | y] = E[beta_1 | y] and SD[beta_2 | y] = SD[beta_1 | y]
obj <- function(par) {
  tp <- topic_post(c(par[1] + delta, par[1] - delta), exp(par[2]) * s2ratio)
  (mix_mean(tp$wg, tp$Eb) - E1)^2 + (mix_sd(tp$wg, tp$Eb, tp$Vb) - S1)^2
}
opt <- optim(c(0.45, log(3)), obj, method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 500))
s_c <- opt$par[1]; c_lvl <- exp(opt$par[2])
tp2 <- topic_post(c(s_c + delta, s_c - delta), c_lvl * s2ratio)

suff <- data.frame(topic = rep(1:2, each = 2), study = rep(1:2, 2), n = n,
                   ybar = c(ybar11, ybar12, s_c + delta, s_c - delta),
                   s2   = c(s2_1, c_lvl * s2ratio))
post <- list(tp1, tp2)

## prior/posterior mean and sd (plus P(>eps|y) for the betas) of every
## parameter, with each study's sufficient statistics attached to its rows
param_table <- function(post, suff) {
  prior_tau_mean <- exp(mlog + slog^2 / 2)
  prior_tau_sd   <- sqrt((exp(slog^2) - 1) * exp(2 * mlog + slog^2))
  prior_bij_sd   <- sqrt(A2 + exp(2 * mlog + 2 * slog^2))   # marginal prior sd of beta_ij
  tab <- NULL
  for (i in 1:2) {
    tp <- post[[i]]; si <- suff[suff$topic == i, ]
    for (j in 1:2)
      tab <- rbind(tab, data.frame(
        param = sprintf("beta_%d%d", i, j), ybar = si$ybar[j], s2 = si$s2[j],
        prior_mean = m0, prior_sd = prior_bij_sd,
        post_mean = mix_mean(tp$wg, tp$Em[, j]),
        post_sd   = mix_sd(tp$wg, tp$Em[, j], tp$Vm[, j]),
        P_eps     = mix_prob(tp$wg, tp$Em[, j], tp$Vm[, j])))
    tab <- rbind(tab, data.frame(
      param = sprintf("beta_%d", i), ybar = NA, s2 = NA,
      prior_mean = m0, prior_sd = sqrt(A2),
      post_mean = mix_mean(tp$wg, tp$Eb), post_sd = mix_sd(tp$wg, tp$Eb, tp$Vb),
      P_eps = mix_prob(tp$wg, tp$Eb, tp$Vb)))
    tab <- rbind(tab, data.frame(
      param = sprintf("tau_%d", i), ybar = NA, s2 = NA,
      prior_mean = prior_tau_mean, prior_sd = prior_tau_sd,
      post_mean = mix_mean(tp$wg, tp$tau), post_sd = mix_sd(tp$wg, tp$tau, 0),
      P_eps = NA))
    for (j in 1:2)
      tab <- rbind(tab, data.frame(
        param = sprintf("sigma2_%d%d", i, j), ybar = NA, s2 = si$s2[j],
        prior_mean = b0 / (a0 - 1),
        prior_sd   = sqrt(b0^2 / ((a0 - 1)^2 * (a0 - 2))),
        post_mean = mix_mean(tp$wg, tp$sg[, j]),
        post_sd   = mix_sd(tp$wg, tp$sg[, j], 0), P_eps = NA))
  }
  tab
}

cat("\n--- parameter table ---\n")
print(param_table(post, suff), digits = 3, row.names = FALSE)
dir.create("results", showWarnings = FALSE)
saveRDS(suff, "results/two_topics_suffstats.rds")

## construction checks
stopifnot(
  abs(mix_prob(tp1$wg, tp1$Em[, 2], tp1$Vm[, 2]) - 0.5) < 1e-5,   # P(beta_12>eps)=0.5
  mix_prob(tp1$wg, tp1$Em[, 1], tp1$Vm[, 1]) > 0.95,              # P(beta_11>eps)>0.95
  abs(mix_mean(tp1$wg, tp1$Eb) - mix_mean(tp2$wg, tp2$Eb)) < 1e-5, # topic means matched
  abs(mix_sd(tp1$wg, tp1$Eb, tp1$Vb) - mix_sd(tp2$wg, tp2$Eb, tp2$Vb)) < 1e-4)

## figure 
## topic-level effect (blue, thick, solid); study-specific effects
## (green, thinner, distinguished by dashing)
col_meta  <- "#2a78d6"; col_study <- "#1baf7a"
COL <- c(col_study, col_study, col_meta)
LWD <- c(1.8, 1.8, 4); LTY <- c(2, 3, 1)

draw_panel <- function(tp, i, xlim, ylim, show_xlab) {
  x <- seq(xlim[1], xlim[2], length.out = 512)
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = if (show_xlab) "effect of interest" else "",
       ylab = "posterior density", main = paste("Topic", i), bty = "l", las = 1)
  abline(v = eps, col = "grey70", lty = 3, lwd = 0.8)
  lines(x, mix_dens(x, tp$wg, tp$Em[, 1], tp$Vm[, 1]), col = COL[1], lwd = LWD[1], lty = LTY[1])
  lines(x, mix_dens(x, tp$wg, tp$Em[, 2], tp$Vm[, 2]), col = COL[2], lwd = LWD[2], lty = LTY[2])
  lines(x, mix_dens(x, tp$wg, tp$Eb, tp$Vb),           col = COL[3], lwd = LWD[3], lty = LTY[3])
  
  ## label: beta_sub: P(beta_sub > eps | y) = p
  ## legend labels: beta_i1, beta_i2, beta_i
  mk <- function(sub) {
    as.expression(
      bquote(beta[.(sub)])
    )
  }
  
  legend(
    "topright",
    bty = "n",
    cex = 0.85,
    col = COL,
    lwd = LWD,
    lty = LTY,
    legend = c(
      mk(paste0(i, "1")),
      mk(paste0(i, "2")),
      mk(i)
    )
  )
}

dir.create("figures", showWarnings = FALSE)
pdf("figures/fig_two_topics.pdf", width = 5, height = 7.2)
op <- par(mfrow = c(2, 1), mar = c(3.6, 4.2, 2.2, 1), mgp = c(2.2, 0.7, 0))
xlim <- c(-1.4, 2.2)
xg <- seq(xlim[1], xlim[2], length.out = 512)
ymax <- max(vapply(post, function(tp)
  max(mix_dens(xg, tp$wg, tp$Em[, 1], tp$Vm[, 1]),
      mix_dens(xg, tp$wg, tp$Em[, 2], tp$Vm[, 2])), 0))
for (i in 1:2) draw_panel(post[[i]], i, xlim, c(0, 1.3 * ymax), i == 2)
par(op); dev.off()
cat("\nWrote: results/two_topics_suffstats.rds, figures/fig_two_topics.pdf\n")
