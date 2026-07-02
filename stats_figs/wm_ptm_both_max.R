rm(list = ls())
graphics.off()

library(tidyverse)
library(car)
library(lmtest)
library(boot)

datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/ptm-wm.csv"
data <- read_csv(datapath)

# Create max-noise predictor
# Important: this z-scores Na and Af within each group, matching your Python code.
data_max <- data %>%
  group_by(group) %>%
  mutate(
    Na_z = as.numeric(scale(`Na-fix`)),
    Af_z = as.numeric(scale(`Af-fix`)),
    max_z = pmax(Na_z, Af_z),
    Na_larger = Na_z > Af_z,
    min_z = pmin(Na_z, Af_z)
  ) %>%
  ungroup() %>%
  mutate(
    group = factor(group),
    max_z_c = as.numeric(scale(max_z, center = TRUE, scale = FALSE))
  )

# Regression: group, max-noise, and interaction
model_max <- lm(
  kavgvar ~ group * max_z_c,
  data = data_max
)

summary(model_max)
confint(model_max, level = 0.95)

# R2 and adjusted R2
summary(model_max)$r.squared
summary(model_max)$adj.r.squared

# Bootstrap 95% CI for R2
r2_boot_fun <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(kavgvar ~ group * max_z_c, data = d)
  summary(m)$r.squared
}

set.seed(123)
boot_r2_max <- boot(data = data_max, statistic = r2_boot_fun, R = 5000)
boot.ci(boot_r2_max, type = c("perc", "bca"))

# Assumption checks

# Normality of residuals
shapiro.test(resid(model_max))

qqnorm(resid(model_max))
qqline(resid(model_max))

# Linearity and homoscedasticity: residuals vs fitted
plot(
  fitted(model_max),
  resid(model_max),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs fitted"
)
abline(h = 0, lty = 2)

# Homoscedasticity: Breusch-Pagan test
bptest(model_max)

# Multicollinearity
vif(model_max)

# Influential observations
cooks_d <- cooks.distance(model_max)

plot(
  cooks_d,
  type = "h",
  ylab = "Cook's distance",
  main = "Cook's distance"
)
abline(h = 4 / nrow(data_max), col = "red", lty = 2)

which(cooks_d > 4 / nrow(data_max))

# Standard diagnostic plots
par(mfrow = c(2, 2))
plot(model_max)
par(mfrow = c(1, 1))

# HC3
library(sandwich)
library(lmtest)

coeftest(
  model_max,
  vcov = vcovHC(model_max, type = "HC3")
)
