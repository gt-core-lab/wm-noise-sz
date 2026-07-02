rm(list = ls())
graphics.off()

library(tidyverse)
library(afex)
library(effectsize)

# load
datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/tdcs.csv"
data <- read_csv(datapath)

# reorganize WM data to long format
wm_long <- data %>%
  select(
    SID,
    FCZ_Var2, FCZ_Var4,
    PZ_Var2, PZ_Var4,
    Sh_Var2, Sh_Var4
  ) %>%
  pivot_longer(
    cols = -SID,
    names_to = c("Loc", "setsize"),
    names_pattern = "(FCZ|PZ|Sh)_Var(2|4)",
    values_to = "var"
  ) %>%
  mutate(
    SID = factor(SID),
    Loc = factor(Loc, levels = c("Sh", "FCZ", "PZ")),
    setsize = factor(setsize, levels = c("2", "4"))
  )

head(wm_long)

# 2 x 3 repeated-measures ANOVA:
# setsize: 2, 4
# Loc: Sham, Frontal, Parietal
aov_tdcs_wm <- aov_ez(
  id = "SID",
  dv = "var",
  within = c("setsize", "Loc"),
  data = wm_long,
  type = 3,
  anova_table = list(correction = "none", es = "pes")
)

print(aov_tdcs_wm)

# Partial eta-squared with 95% CI
eta_ci_tdcs_wm <- eta_squared(
  aov_tdcs_wm,
  partial = TRUE,
  ci = 0.95
)

print(eta_ci_tdcs_wm)

# Assumption checks

# Residual normality
resid_tdcs_wm <- residuals(aov_tdcs_wm$lm)

shapiro.test(resid_tdcs_wm)

qqnorm(resid_tdcs_wm)
qqline(resid_tdcs_wm)

# Sphericity (Mauchly's tests and GG/HF corrections)
summary(aov_tdcs_wm)

# Inspect variances by condition
aggregate(var ~ Loc + setsize, data = wm_long, var)

# Outliers by condition
boxplot(
  var ~ Loc * setsize,
  data = wm_long,
  xlab = "Stimulation condition × set size",
  ylab = "WM recall variance",
  las = 2
)

