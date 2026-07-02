rm(list = ls())
graphics.off()

library(tidyverse)
library(car)
library(lmtest)
library(sandwich)
library(lme4)
library(lmerTest)
library(performance)
library(boot)

# load
datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/tdcs.csv"
data <- read_csv(datapath)

# regression per stimulation site

# Frontal vs sham
model_Fr <- lm(
  FrSh_VarAve ~ FrSh_Na + FrSh_Nm + FrSh_Af,
  data = data
)

summary(model_Fr)
confint(model_Fr, level = 0.95)

# Parietal vs sham
model_Pz <- lm(
  PzSh_VarAve ~ PSh_Na + PSh_Nm + PSh_Af,
  data = data
)

summary(model_Pz)
confint(model_Pz, level = 0.95)

# function for regression diagnostics
check_lm_assumptions <- function(model, data_name = "model") {
  
  cat("\n=============================\n")
  cat("Diagnostics for:", data_name, "\n")
  cat("=============================\n")
  
  # Normality of residuals
  cat("\nShapiro-Wilk test:\n")
  print(shapiro.test(resid(model)))
  
  qqnorm(resid(model), main = paste("Q-Q plot:", data_name))
  qqline(resid(model))
  
  # Linearity / homoscedasticity plot
  plot(
    fitted(model),
    resid(model),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste("Residuals vs fitted:", data_name)
  )
  abline(h = 0, lty = 2)
  
  # Breusch-Pagan
  cat("\nBreusch-Pagan test:\n")
  print(bptest(model))
  
  # If heteroscedasticity is present, HC3 robust SE
  cat("\nHC3 robust coefficient tests:\n")
  print(coeftest(model, vcov = vcovHC(model, type = "HC3")))
  
  # Multicollinearity
  cat("\nVIF:\n")
  print(vif(model))
  
  # Cook's distance
  cooks_d <- cooks.distance(model)
  threshold <- 4 / nobs(model)
  
  plot(
    cooks_d,
    type = "h",
    ylab = "Cook's distance",
    main = paste("Cook's distance:", data_name)
  )
  abline(h = threshold, col = "red", lty = 2)
  
  cat("\nInfluential observations using Cook's D > 4/n:\n")
  print(which(cooks_d > threshold))
  
  # Standard diagnostic plots
  par(mfrow = c(2, 2))
  plot(model)
  par(mfrow = c(1, 1))
}

check_lm_assumptions(model_Fr, "Frontal vs sham")
check_lm_assumptions(model_Pz, "Parietal vs sham")

# Bootstrap R2 CIs for separate regressions

r2_boot_Fr <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(FrSh_VarAve ~ FrSh_Na + FrSh_Nm + FrSh_Af, data = d)
  summary(m)$r.squared
}

r2_boot_Pz <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(PzSh_VarAve ~ PSh_Na + PSh_Nm + PSh_Af, data = d)
  summary(m)$r.squared
}

set.seed(123)

boot_Fr <- boot(data = data, statistic = r2_boot_Fr, R = 5000)
boot_Pz <- boot(data = data, statistic = r2_boot_Pz, R = 5000)

boot.ci(boot_Fr, type = c("perc", "bca"))
boot.ci(boot_Pz, type = c("perc", "bca"))


# Long-format data for combined site model
tdcs_long <- bind_rows(
  data %>%
    transmute(
      SID = SID,
      Site = "Frontal",
      dWM = FrSh_VarAve,
      dNa = FrSh_Na,
      dNm = FrSh_Nm,
      dAf = FrSh_Af
    ),
  data %>%
    transmute(
      SID = SID,
      Site = "Parietal",
      dWM = PzSh_VarAve,
      dNa = PSh_Na,
      dNm = PSh_Nm,
      dAf = PSh_Af
    )
) %>%
  mutate(
    SID = factor(SID),
    Site = factor(Site, levels = c("Frontal", "Parietal")),
    dNa_c = as.numeric(scale(dNa, center = TRUE, scale = FALSE)),
    dNm_c = as.numeric(scale(dNm, center = TRUE, scale = FALSE)),
    dAf_c = as.numeric(scale(dAf, center = TRUE, scale = FALSE))
  )


# Mixed-effects model:
# Does PTM-WM relationship differ by stimulation site?

model_site <- lmer(
  dWM ~ Site * (dNa_c + dNm_c + dAf_c) + (1 | SID),
  data = tdcs_long,
  REML = FALSE
)

summary(model_site)
confint(model_site, method = "Wald")

# Model R2
r2(model_site)

# Multicollinearity for mixed model
check_collinearity(model_site)

# Residual normality
shapiro.test(resid(model_site))

qqnorm(resid(model_site))
qqline(resid(model_site))

# Residuals vs fitted
plot(
  fitted(model_site),
  resid(model_site),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Mixed model: residuals vs fitted"
)
abline(h = 0, lty = 2)

# Approximate heteroscedasticity check using lm version
model_site_lm <- lm(
  dWM ~ Site * (dNa_c + dNm_c + dAf_c),
  data = tdcs_long
)

bptest(model_site_lm)

coeftest(
  model_site_lm,
  vcov = vcovHC(model_site_lm, type = "HC3")
)

vif(model_site_lm)


# Focused additive-noise model
model_site_Na <- lmer(
  dWM ~ Site * dNa_c + (1 | SID),
  data = tdcs_long,
  REML = FALSE
)

summary(model_site_Na)
confint(model_site_Na, method = "Wald")
r2(model_site_Na)

# Residual diagnostics for focused model
shapiro.test(resid(model_site_Na))

qqnorm(resid(model_site_Na))
qqline(resid(model_site_Na))

plot(
  fitted(model_site_Na),
  resid(model_site_Na),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Focused dNa model: residuals vs fitted"
)
abline(h = 0, lty = 2)


# Bayesian model comparison
# Tests whether adding the Site x dNa interaction improves the model.

library(BayesFactor)

tdcs_long_bf <- tdcs_long %>%
  mutate(
    SID = factor(SID),
    Site = factor(Site)
  ) %>%
  filter(
    complete.cases(dWM, Site, dNa_c, SID)
  ) %>%
  droplevels()

# Check how many rows remain
nrow(tdcs_long_bf)
table(tdcs_long_bf$Site)

# No-interaction model
bf_no_interaction <- lmBF(
  dWM ~ Site + dNa_c + SID,
  data = tdcs_long_bf,
  whichRandom = "SID"
)

# Interaction model
bf_interaction <- lmBF(
  dWM ~ Site * dNa_c + SID,
  data = tdcs_long_bf,
  whichRandom = "SID"
)

BF10_interaction <- bf_interaction / bf_no_interaction
BF01_interaction <- 1 / BF10_interaction

BF10_interaction
BF01_interaction


# Direct bootstrap comparison of Na and Af slopes across sites
# Compares simple dNoise -> dWM slopes for Parietal vs Frontal.
# Bootstrap is done at the SUBJECT level to preserve within-subject dependence.

set.seed(123)

subject_ids <- unique(tdcs_long$SID)

boot_slope_diff_NaAf <- boot(
  data = tibble(SID = subject_ids),
  statistic = function(subject_data, indices) {
    
    sampled_ids <- subject_data$SID[indices]
    
    d <- map_dfr(sampled_ids, function(id) {
      tdcs_long %>% filter(SID == id)
    })
    
    # Na slopes
    m_fr_na <- lm(dWM ~ dNa, data = d %>% filter(Site == "Frontal"))
    m_pz_na <- lm(dWM ~ dNa, data = d %>% filter(Site == "Parietal"))
    
    slope_fr_na <- coef(m_fr_na)["dNa"]
    slope_pz_na <- coef(m_pz_na)["dNa"]
    
    # Af slopes
    m_fr_af <- lm(dWM ~ dAf, data = d %>% filter(Site == "Frontal"))
    m_pz_af <- lm(dWM ~ dAf, data = d %>% filter(Site == "Parietal"))
    
    slope_fr_af <- coef(m_fr_af)["dAf"]
    slope_pz_af <- coef(m_pz_af)["dAf"]
    
    c(
      Na_slope_Frontal = slope_fr_na,
      Na_slope_Parietal = slope_pz_na,
      Na_slope_diff_Pz_minus_Fr = slope_pz_na - slope_fr_na,
      Af_slope_Frontal = slope_fr_af,
      Af_slope_Parietal = slope_pz_af,
      Af_slope_diff_Pz_minus_Fr = slope_pz_af - slope_fr_af
    )
  },
  R = 5000
)

# Observed slopes and slope differences
obs_m_fr_na <- lm(dWM ~ dNa, data = tdcs_long %>% filter(Site == "Frontal"))
obs_m_pz_na <- lm(dWM ~ dNa, data = tdcs_long %>% filter(Site == "Parietal"))

obs_m_fr_af <- lm(dWM ~ dAf, data = tdcs_long %>% filter(Site == "Frontal"))
obs_m_pz_af <- lm(dWM ~ dAf, data = tdcs_long %>% filter(Site == "Parietal"))

obs_slope_fr_na <- coef(obs_m_fr_na)["dNa"]
obs_slope_pz_na <- coef(obs_m_pz_na)["dNa"]
obs_slope_diff_na <- obs_slope_pz_na - obs_slope_fr_na

obs_slope_fr_af <- coef(obs_m_fr_af)["dAf"]
obs_slope_pz_af <- coef(obs_m_pz_af)["dAf"]
obs_slope_diff_af <- obs_slope_pz_af - obs_slope_fr_af

observed_slopes <- tibble(
  Parameter = c("Na", "Af"),
  Frontal_slope = c(obs_slope_fr_na, obs_slope_fr_af),
  Parietal_slope = c(obs_slope_pz_na, obs_slope_pz_af),
  Difference_Pz_minus_Fr = c(obs_slope_diff_na, obs_slope_diff_af)
)

observed_slopes

# Bootstrap CIs
boot.ci(boot_slope_diff_NaAf, type = c("perc", "bca"), index = 3) # Na difference
boot.ci(boot_slope_diff_NaAf, type = c("perc", "bca"), index = 6) # Af difference

# Two-sided bootstrap p-values for slope differences
boot_diffs_na <- boot_slope_diff_NaAf$t[, 3]
boot_diffs_af <- boot_slope_diff_NaAf$t[, 6]

p_boot_slope_diff_na <- 2 * min(
  mean(boot_diffs_na <= 0, na.rm = TRUE),
  mean(boot_diffs_na >= 0, na.rm = TRUE)
)

p_boot_slope_diff_af <- 2 * min(
  mean(boot_diffs_af <= 0, na.rm = TRUE),
  mean(boot_diffs_af >= 0, na.rm = TRUE)
)

p_boot_slope_diff_na
p_boot_slope_diff_af

# Compact summary table
bootstrap_slope_summary <- tibble(
  Parameter = c("Na", "Af"),
  Observed_difference_Pz_minus_Fr = c(obs_slope_diff_na, obs_slope_diff_af),
  Bootstrap_p = c(p_boot_slope_diff_na, p_boot_slope_diff_af)
)

bootstrap_slope_summary



