rm(list = ls())
graphics.off()

library(tidyverse)
library(afex)
library(effectsize)
library(car)

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
    setsize = factor(setsize, levels = c("1", "2", "4")),
    group = factor(group),
    studyID = factor(studyID)
  )

# mixed ANOVA: group x set size
aov_group <- aov_ez(
  id = "studyID",
  dv = "kvar",
  within = "setsize",
  between = "group",
  data = df,
  type = 3,
  anova_table = list(correction = "none", es = "pes")
)

print(aov_group)

# partial eta-squared with 95% CI
eta_ci <- eta_squared(
  aov_group,
  partial = TRUE,
  ci = 0.95
)

print(eta_ci)

# Assumption checks

# Normality of residuals
resid_group <- residuals(aov_group$lm)

shapiro.test(resid_group)

qqnorm(resid_group)
qqline(resid_group)

# Sphericity (Mauchly test and GG/HF corrections)
summary(aov_group)

nice(aov_group,
     correction = "GG",
     es = "pes")

# Homogeneity of variance between groups at each set size
leveneTest(kvar ~ group * setsize, data = df)

# Levene test separately for each set size
by(df, df$setsize, function(subdat) {
  print(leveneTest(kvar ~ group, data = subdat))
})

# 4. Inspect variances by group and set size
aggregate(kvar ~ group + setsize, data = df, var)

# 5. Outliers by group and set size
boxplot(
  kvar ~ group * setsize,
  data = df,
  xlab = "Group × Set size",
  ylab = "WM recall variance",
  las = 2
)
