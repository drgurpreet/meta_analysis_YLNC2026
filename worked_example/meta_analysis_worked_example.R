# =============================================================================
# Worked example: random-effects meta-analysis in R
# Dataset : BCG vaccine trials for tuberculosis (Colditz et al., 1994)
#           13 RCTs, binary outcome (TB cases), effect measure = risk ratio
# Packages: metafor (analysis) + tidyverse (data handling, ggplot2 graphics)
# Style   : tidyverse, native pipe |>
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
# install.packages(c("metafor", "tidyverse"))
library(metafor)
library(tidyverse)

theme_set(theme_minimal(base_size = 12))

# Helper: turn metafor output (predict(), leave1out()) into a tibble with
# snake_case names, e.g. ci.lb -> ci_lb
tidy_mf <- \(x) x |>
  as.data.frame() |>
  as_tibble() |>
  rename_with(\(nm) str_replace_all(nm, "\\.", "_"))

# ---- 1. Data ----------------------------------------------------------------
# dat.bcg ships with metafor (via the metadat package).
# tpos/tneg = TB +/- in vaccinated arm; cpos/cneg = TB +/- in control arm
bcg <- dat.bcg |>
  as_tibble() |>
  mutate(
    study  = str_c(author, ", ", year),
    n_vacc = tpos + tneg,
    n_ctrl = cpos + cneg
  ) |>
  select(trial, study, year, ablat, alloc,
         tpos, n_vacc, cpos, n_ctrl, tneg, cneg)

bcg

# ---- 2. Effect sizes --------------------------------------------------------
# escalc() adds yi (log risk ratio) and vi (its sampling variance)
bcg_es <- escalc(
  measure = "RR",
  ai = tpos, bi = tneg,   # vaccinated: events, non-events
  ci = cpos, di = cneg,   # control:    events, non-events
  data = bcg,
  slab = study
)

# Per-study RR with 95% CI, back-transformed for reading
bcg_es |>
  as_tibble() |>
  mutate(
    rr    = exp(yi),
    lower = exp(yi - 1.96 * sqrt(vi)),
    upper = exp(yi + 1.96 * sqrt(vi))
  ) |>
  select(study, tpos, n_vacc, cpos, n_ctrl, rr, lower, upper) |>
  mutate(across(rr:upper, \(x) round(x, 2)))

# ---- 3. Pooled estimates: fixed (common) vs random effects ----------------
res_fe <- rma(yi, vi, data = bcg_es, method = "EE")    # common effect
res_re <- rma(yi, vi, data = bcg_es, method = "REML")  # random effects

res_re

# Tidy side-by-side comparison on the RR scale
bind_rows(
  `Common effect`  = predict(res_fe, transf = exp) |> tidy_mf(),
  `Random effects` = predict(res_re, transf = exp) |> tidy_mf(),
  .id = "model"
) |>
  select(model, rr = pred, ci_lb, ci_ub, pi_lb, pi_ub) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

# ---- 4. Forest plot ---------------------------------------------------------
forest(
  res_re,
  atransf = exp,
  at      = log(c(0.05, 0.25, 1, 4)),
  xlim    = c(-16, 6),
  ilab    = cbind(bcg$tpos, bcg$n_vacc, bcg$cpos, bcg$n_ctrl),
  ilab.xpos = c(-9.5, -8, -6, -4.5),
  header  = c("Study", "RR [95% CI]"),
  addpred = TRUE,          # show the 95% prediction interval
  shade   = TRUE,
  mlab    = "Random-effects model (REML)"
)
text(c(-9.5, -8, -6, -4.5), res_re$k + 2, c("TB+", "N", "TB+", "N"), font = 2, cex = 0.8)
text(c(-8.75, -5.25),       res_re$k + 3, c("Vaccinated", "Control"), font = 2, cex = 0.8)
# HOW TO READ
#   - Each row = one study: square = its RR, line = its 95% CI.
#   - Square size = weight in the pooled analysis.
#   - Vertical line at RR = 1 = no effect; left favours vaccination.
#   - Diamond = pooled RR (centre) and its 95% CI (width).
#   - Dotted line through the diamond = 95% prediction interval (PI).
#   - x-axis is log scale: RR 0.5 and RR 2 are equally far from 1.
# WHAT WE INFER
#   - 11 of 13 trials favour BCG; only TPT Madras (1.01) and
#     Comstock & Webster (1.56) do not.
#   - 5 trials are non-significant on their own, yet the pooled RR
#     0.49 (0.34 to 0.70) is clearly below 1: pooling adds precision.
#   - Many CIs do not overlap (e.g. Hart & Sutherland 0.24 vs TPT Madras
#     1.01): visible heterogeneity.
#   - RE weights are fairly even (3.8% to 10.2%) because tau^2 is large;
#     under the common-effect model TPT Madras alone carries 41% of the weight.
#   - PI 0.15 to 1.55 crosses 1: BCG works on average, but in some settings
#     it may give little or no protection. Report both CI and PI.


# ---- 5. Heterogeneity -------------------------------------------------------
# Q test, tau^2, I^2, H^2 come with the model; confint() gives their CIs
tibble(
  statistic = c("Q", "df", "p (Q)", "tau^2", "tau", "I^2 (%)", "H^2"),
  value     = c(res_re$QE, res_re$k - 1, res_re$QEp,
                res_re$tau2, sqrt(res_re$tau2), res_re$I2, res_re$H2)
) |>
  mutate(value = if_else(statistic == "p (Q)",
                         format.pval(value, eps = 1e-4),
                         as.character(round(value, 2))))

confint(res_re)   # CIs for tau^2, tau, I^2, H^2

# Prediction interval: where the true effect of a *new* trial is likely to lie
predict(res_re, transf = exp)

# Exploring heterogeneity: does distance from the equator explain it?
res_mr <- rma(yi, vi, mods = ~ ablat, data = bcg_es, method = "REML")
res_mr   # look at R^2 and the residual I^2

# Bubble plot of the meta-regression (ggplot2)
lat_grid <- tibble(ablat = seq(10, 60, by = 1))

lat_pred <- predict(res_mr, newmods = lat_grid$ablat, transf = exp) |>
  tidy_mf() |>
  bind_cols(lat_grid)

bcg_es |>
  as_tibble() |>
  mutate(rr = exp(yi), weight = 1 / vi) |>
  ggplot(aes(ablat, rr)) +
  geom_ribbon(data = lat_pred, aes(y = pred, ymin = ci_lb, ymax = ci_ub),
              fill = "grey85") +
  geom_line(data = lat_pred, aes(y = pred), linewidth = 0.8) +
  geom_point(aes(size = weight), shape = 21, fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  scale_size_area(max_size = 12, guide = "none") +
  labs(x = "Absolute latitude (degrees)", y = "Risk ratio (log scale)",
       title = "Meta-regression: BCG efficacy vs latitude",
       subtitle = "Point size proportional to inverse-variance weight")
# HOW TO READ
#   - Each bubble = one trial at its latitude (x) and RR (y, log scale);
#     bubble size = weight.
#   - Solid line = fitted meta-regression; grey band = its 95% CI.
#   - Dashed line at RR = 1 = no effect.
#   - A sloping line means the moderator explains heterogeneity.
# WHAT WE INFER
#   - Each extra degree of latitude multiplies the RR by 0.97
#     (about 2.9% more protection per degree, p < 0.0001).
#   - Near the equator (13-19 deg) BCG shows little effect; above 50 deg the
#     RR is around 0.2-0.3.
#   - Latitude explains ~76% of tau^2; residual I^2 falls from 92% to 68%.
#   - Caution: ecological, across-study association. Latitude is a proxy
#     (plausibly environmental mycobacteria). With k = 13, fit at most one
#     moderator (~10 studies per covariate).


# ---- 6. Funnel plot ---------------------------------------------------------
# Standard funnel plot
funnel(res_re, atransf = exp, xlab = "Risk ratio (log scale)")

# Contour-enhanced funnel: helps tell publication bias from other asymmetry
funnel(res_re, level = c(90, 95, 99),
       shade = c("white", "grey80", "grey60"),
       refline = 0, legend = TRUE,
       atransf = exp, xlab = "Risk ratio (log scale)")
# HOW TO READ
#   - x = study effect (log RR); y = standard error, inverted, so large
#     precise studies are at the top.
#   - Vertical line = pooled estimate; triangle = where 95% of studies should
#     fall if there were no bias and no heterogeneity.
#   - Expect a symmetric inverted funnel. A gap in one bottom corner
#     suggests missing small studies.
#   - Contour version: reference line at RR = 1; grey zones = p < .10, .05,
#     .01. A gap in the white (non-significant) zone points to publication
#     bias; a gap in a shaded zone points to other causes.
# WHAT WE INFER
#   - Large trials at the top are spread widely left to right: that is
#     heterogeneity, not publication bias.
#   - Many points lie outside the triangle, as expected with I^2 = 92%.
#   - Small trials mostly lie to the left (strong benefit), one to the right
#     (Comstock & Webster): only a mild hint of asymmetry.
#   - Small trials sit in both shaded and white zones: no clear sign that
#     non-significant studies are missing.


# ---- 7. Egger's regression test --------------------------------------------
egger <- regtest(res_re, model = "lm", predictor = "sei")  # classic Egger
egger

# Rank correlation (Begg) for comparison
ranktest(res_re)

# Trim-and-fill (sensitivity analysis for small-study effects)
res_tf <- trimfill(res_re)
res_tf
funnel(res_tf, atransf = exp, legend = TRUE, xlab = "Risk ratio (log scale)")
# HOW TO READ
#   - Egger: regression of effect on its SE; a slope different from 0 means
#     small studies differ systematically from large ones. Begg: same idea
#     with rank correlation.
#   - Trim-and-fill funnel: open circles = imputed "missing" studies;
#     vertical line = adjusted pooled estimate.
# WHAT WE INFER
#   - Egger p = 0.19, Begg p = 0.95: no strong evidence of small-study effects.
#   - Trim-and-fill imputes 1 study on the right (towards no effect); pooled
#     RR moves only from 0.49 to 0.52 (0.37 to 0.74). Conclusion unchanged.


# Caution for teaching: with k = 13 these tests have low power.
# Rule of thumb (Cochrane): do not test for asymmetry with fewer than 10 studies.

# ---- 8. Leave-one-out sensitivity analysis ---------------------------------
loo <- leave1out(res_re, transf = exp) |>
  tidy_mf() |>
  mutate(study = bcg$study) |>
  select(omitted = study, rr = estimate, ci_lb, ci_ub, I2, tau2)

loo |> mutate(across(where(is.numeric), \(x) round(x, 2)))

pooled_rr <- exp(res_re$b[1])

loo |>
  mutate(omitted = fct_reorder(omitted, rr)) |>
  ggplot(aes(rr, omitted)) +
  annotate("rect", xmin = exp(res_re$ci.lb), xmax = exp(res_re$ci.ub),
           ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = pooled_rr, linetype = "dashed") +
  geom_vline(xintercept = 1, colour = "firebrick") +
  geom_pointrange(aes(xmin = ci_lb, xmax = ci_ub)) +
  geom_text(aes(x = ci_ub, label = str_c("I² = ", round(I2), "%")),
            hjust = -0.15, size = 3, colour = "grey30") +
  scale_x_log10(breaks = c(0.3, 0.4, 0.5, 0.7, 1),
                expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Pooled RR when study is omitted (log scale)", y = "Study omitted",
       title = "Leave-one-out sensitivity analysis",
       subtitle = "Dashed line and band = overall pooled RR and 95% CI")
# HOW TO READ
#   - Each row = the meta-analysis re-run without the named study.
#   - Dashed line + grey band = original pooled RR and CI; red line = RR 1.
#   - Rows sorted by new RR: top/bottom rows are the most influential.
#   - Right-hand labels = I^2 after removal; a big drop marks a source of
#     heterogeneity.
#   - Look for: any estimate leaving the band, or any CI crossing 1.
# WHAT WE INFER
#   - Pooled RR ranges only 0.45 to 0.53; all CIs stay clear of 1: robust.
#   - Dropping TPT Madras gives the most favourable RR (0.45) and the lowest
#     I^2 (87%); dropping Hart & Sutherland gives the least favourable (0.53).
#   - I^2 stays above 85% whichever study is removed: heterogeneity is not
#     due to one outlier.


# Influence diagnostics (Cook's distance, DFFITS, etc.) and Baujat plot
influence(res_re) |> plot()
# HOW TO READ
#   - One point per study in 8 panels. rstudent beyond +/-2 = outlier;
#     dffits / cook.d = change in pooled estimate when removed;
#     cov.r < 1 = removal improves precision; tau2.del / QE.del = tau^2 and Q
#     after removal; hat / weight = leverage and weight.
#   - Red points = flagged as influential (Viechtbauer & Cheung, 2010).
# WHAT WE INFER
#   - No study is flagged; no residual exceeds +/-2.
#   - Highest Cook's distance: Hart & Sutherland (0.23), TPT Madras (0.20).
#     They also cut Q the most (to 97 and 68, from 152).

baujat(res_re)
# HOW TO READ
#   - x = contribution to overall heterogeneity (Q); y = influence on the
#     pooled result. Top-right = studies that do both; inspect them.
# WHAT WE INFER
#   - Hart & Sutherland (4) and TPT Madras (8) are top-right, at opposite
#     latitudes (52 vs 13 deg): again pointing to latitude.
#   - Stein & Aronson and Aronson sit near the origin: consistent studies.


# ---- 9. Session info --------------------------------------------------------
sessionInfo()
