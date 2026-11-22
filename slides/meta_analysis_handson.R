# =============================================================================
# Meta-analysis: hands-on companion script
# Session: Meta-Analysis - When and How (YLNC)
# Department of Community Medicine, AFMC Pune
#
# Works through the same BCG vaccine example used in the slides.
# Style: tidyverse, native pipe |>, one step per section.
# Packages: tidyverse, metafor (install.packages(c("tidyverse", "metafor")))
# =============================================================================

library(tidyverse)
library(metafor)

# -----------------------------------------------------------------------------
# 1. Data: one row per study (tidy)
# -----------------------------------------------------------------------------
# 13 trials of BCG vaccine against tuberculosis (Colditz et al., JAMA 1994)
# tpos/tneg = TB / no TB in vaccinated arm; cpos/cneg = same in control arm
bcg <- dat.bcg |>
  as_tibble() |>
  mutate(study = str_c(author, " (", year, ")"))

bcg |> glimpse()

# -----------------------------------------------------------------------------
# 2. Stage 1: effect size and variance for each study (log risk ratio)
# -----------------------------------------------------------------------------
# By hand, so the formula is visible
bcg_es <- bcg |>
  mutate(
    risk_vacc = tpos / (tpos + tneg),
    risk_ctrl = cpos / (cpos + cneg),
    yi  = log(risk_vacc / risk_ctrl),
    vi  = 1 / tpos - 1 / (tpos + tneg) + 1 / cpos - 1 / (cpos + cneg),
    sei = sqrt(vi),
    rr  = exp(yi),
    rr_lo = exp(yi - 1.96 * sei),
    rr_hi = exp(yi + 1.96 * sei)
  )

bcg_es |> select(study, rr, rr_lo, rr_hi)

# Same thing with metafor::escalc() - check they agree
bcg |>
  escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = _) |>
  as_tibble() |>
  select(author, year, yi, vi)

# EXERCISE: change measure to "OR" or "RD". How do the effect sizes change?

# -----------------------------------------------------------------------------
# 3. Stage 2: inverse-variance pooling by hand (fixed effect)
# -----------------------------------------------------------------------------
pooled_fe <- bcg_es |>
  mutate(wi = 1 / vi) |>
  summarise(
    theta = sum(wi * yi) / sum(wi),
    se    = sqrt(1 / sum(wi)),
    Q     = sum(wi * (yi - theta)^2),
    k     = n()
  ) |>
  mutate(
    RR    = exp(theta),
    lower = exp(theta - 1.96 * se),
    upper = exp(theta + 1.96 * se),
    I2    = pmax(0, (Q - (k - 1)) / Q) * 100
  )

pooled_fe

# -----------------------------------------------------------------------------
# 4. Fixed vs random effects with metafor
# -----------------------------------------------------------------------------
fe <- rma(yi, vi, data = bcg_es, method = "FE")
re <- rma(yi, vi, data = bcg_es, method = "REML", slab = study)

fe
re

# Random effects with Hartung-Knapp adjustment (recommended with few studies)
re_hk <- rma(yi, vi, data = bcg_es, method = "REML", test = "knha", slab = study)
re_hk

# Pooled RR with confidence and prediction intervals, back-transformed
predict(re, transf = exp)

# Compare estimators of tau^2 in one tidy table
c("DL", "REML", "PM") |>
  set_names() |>
  map(\(m) rma(yi, vi, data = bcg_es, method = m)) |>
  imap(\(fit, m) tibble(
    method = m,
    RR     = exp(fit$b[1]),
    lower  = exp(fit$ci.lb),
    upper  = exp(fit$ci.ub),
    tau2   = fit$tau2,
    I2     = fit$I2
  )) |>
  list_rbind()

# Weights under each model
bcg_es |>
  mutate(w_fixed = weights(fe), w_random = weights(re)) |>
  select(study, w_fixed, w_random) |>
  arrange(desc(w_fixed))

# -----------------------------------------------------------------------------
# 5. Forest plot (ggplot2)
# -----------------------------------------------------------------------------
forest_df <- bcg_es |>
  mutate(weight = weights(re)) |>
  arrange(year) |>
  mutate(study = fct_rev(fct_inorder(study)))

pooled_row <- tibble(
  study = "Random effects (REML)",
  rr = exp(re$b[1]), rr_lo = exp(re$ci.lb), rr_hi = exp(re$ci.ub)
)

forest_df |>
  ggplot(aes(rr, study)) +
  geom_vline(xintercept = 1, colour = "grey50") +
  geom_vline(xintercept = pooled_row$rr, colour = "#C8553D", linetype = "dashed") +
  geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi), height = 0) +
  geom_point(aes(size = weight), shape = 15, colour = "#12355B") +
  scale_x_log10() +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(
    x = "Risk ratio (log scale)", y = NULL,
    title = "BCG vaccine and tuberculosis",
    subtitle = str_glue(
      "Pooled RR {round(pooled_row$rr, 2)} ",
      "(95% CI {round(pooled_row$rr_lo, 2)} to {round(pooled_row$rr_hi, 2)}); ",
      "I² = {round(re$I2)}%"
    )
  ) +
  theme_minimal(base_size = 13)

# metafor's built-in forest plot, for comparison
forest(re, atransf = exp, header = c("Study", "RR [95% CI]"), addpred = TRUE)

# -----------------------------------------------------------------------------
# 6. Explaining heterogeneity
# -----------------------------------------------------------------------------
# 6a. Subgroup analysis by allocation method
bcg_es |>
  nest(.by = alloc) |>
  mutate(
    fit = map(data, \(d) rma(yi, vi, data = d)),
    k   = map_int(fit, \(f) f$k),
    RR  = map_dbl(fit, \(f) exp(f$b[1])),
    lower = map_dbl(fit, \(f) exp(f$ci.lb)),
    upper = map_dbl(fit, \(f) exp(f$ci.ub)),
    I2  = map_dbl(fit, \(f) f$I2)
  ) |>
  select(-data, -fit)

# Formal test for subgroup differences
rma(yi, vi, mods = ~ alloc, data = bcg_es)

# 6b. Meta-regression on absolute latitude
mr <- rma(yi, vi, mods = ~ ablat, data = bcg_es)
mr
mr$R2   # % of between-study variance explained

bcg_es |>
  ggplot(aes(ablat, rr)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = 1 / vi), shape = 21, fill = "#1F6F8B", colour = "white") +
  geom_smooth(aes(weight = 1 / vi), method = "lm", formula = y ~ x,
              colour = "#C8553D", fill = "#C8553D", alpha = 0.15) +
  scale_y_log10() +
  scale_size_continuous(range = c(2, 12), guide = "none") +
  labs(x = "Absolute latitude (degrees)", y = "Risk ratio (log scale)") +
  theme_minimal(base_size = 13)
# Note: geom_smooth() gives a weighted least-squares line for display only;
# use the rma() model above for inference.

# -----------------------------------------------------------------------------
# 7. Sensitivity analysis: leave one out
# -----------------------------------------------------------------------------
leave1out(re) |>
  as_tibble() |>
  mutate(
    omitted = bcg_es$study,
    across(c(estimate, ci.lb, ci.ub), exp)
  ) |>
  select(omitted, RR = estimate, lower = ci.lb, upper = ci.ub, I2) |>
  arrange(RR)

# -----------------------------------------------------------------------------
# 8. Small-study effects / publication bias
# -----------------------------------------------------------------------------
funnel(re)
regtest(re)        # Egger's regression test
ranktest(re)       # Begg's rank correlation test
tf <- trimfill(re)
tf
funnel(tf)         # filled (imputed) studies shown as open circles

# -----------------------------------------------------------------------------
# 9. EXERCISES
# -----------------------------------------------------------------------------
# 1. Re-run the main analysis with measure = "OR". Does the conclusion change?
# 2. Use method = "DL" and test = "knha". What happens to the CI width?
# 3. Drop the trials using alternate allocation. Re-estimate RR and I^2.
# 4. Add year as a second moderator to the meta-regression. Is it useful?
# 5. Write a two-line GRADE judgement for "inconsistency" using I^2,
#    tau^2 and the prediction interval.
