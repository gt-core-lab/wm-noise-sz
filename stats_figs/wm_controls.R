rm(list = ls())
graphics.off()

library(tidyverse)
library(afex)
library(effectsize)

# load
datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/ptm-wm.csv"
data <- read_csv(datapath)

# reorganize wide to long
df <- data %>%
  select(studyID, group, k1var, k2var, k4var) %>%
  rename(
    kvar1 = k1var,
    kvar2 = k2var,
    kvar4 = k4var
  ) %>%
  pivot_longer(
    cols = starts_with("kvar"),
    names_to = "setsize",
    values_to = "kvar"
  ) %>%
  mutate(
    setsize = str_remove(setsize, "kvar"),
    setsize = factor(setsize, levels = c("1", "2", "4"))
  )

# select CO / neurotypical data
df_CO <- df %>%
  filter(group == "CO")

# repeated-measures ANOVA
aov_CO <- aov_ez(
  id = "studyID",
  dv = "kvar",
  within = "setsize",
  data = df_CO,
  type = 3,
  anova_table = list(correction = "none", es = "pes")
)

print(aov_CO)

# partial eta-squared with 95% CI
eta_ci <- eta_squared(
  aov_CO,
  partial = TRUE,
  ci = 0.95
)

print(eta_ci)

# assumption checks 

# Normality of residuals
resid_CO <- residuals(aov_CO$lm)

shapiro.test(resid_CO)

qqnorm(resid_CO)
qqline(resid_CO)

# Sphericity: Mauchly test + GG/HF corrections
summary(aov_CO)

# Variance by set size
aggregate(kvar ~ setsize, data = df_CO, var)

# Outliers by set size
boxplot(
  kvar ~ setsize,
  data = df_CO,
  xlab = "Set size",
  ylab = "WM recall variance"
)
