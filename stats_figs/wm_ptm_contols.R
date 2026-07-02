rm(list = ls())
graphics.off()

library(tidyverse)
library(car)
library(lmtest)
library(effectsize)
library(boot)

# load
datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/ptm-wm.csv"
data <- read_csv(datapath)

# select controls only
data_CO <- data %>%
  filter(group == "CO") %>%
  mutate(
    Na_z = as.numeric(scale(`Na-fix`)),
    Af_z = as.numeric(scale(`Af-fix`)),
    max_z = pmax(Na_z, Af_z),
    min_z = pmin(Na_z, Af_z),
    Na_larger = Na_z > Af_z
  )

# multiple regression
# use backticks because the variable names contain hyphens
model_CO <- lm(
  kavgvar ~ `Na` + `Af` + `Nm`,
  data = data_CO
)

summary(model_CO)

# coefficient confidence intervals
confint(model_CO, level = 0.95)

# confidence intervals on effect size
# function to compute R2 for bootstrap samples
r2_boot_fun <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(kavgvar ~ `Na` + `Af` + `Nm`, data = d)
  summary(m)$r.squared
}

set.seed(123)

boot_r2 <- boot(
  data = data_CO,
  statistic = r2_boot_fun,
  R = 5000
)

# Bootstrap 95% CI for R2
boot.ci(boot_r2, type = c("perc", "bca"))

# Assumption checks
# Linearity and homoscedasticity: residuals vs fitted
plot(
  fitted(model_CO),
  resid(model_CO),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs fitted"
)
abline(h = 0, lty = 2)

# Normality of residuals
shapiro.test(resid(model_CO))

qqnorm(resid(model_CO))
qqline(resid(model_CO))

# Homoscedasticity: Breusch-Pagan test
bptest(model_CO)

# Multicollinearity
vif(model_CO)

# Influential observations
cooks_d <- cooks.distance(model_CO)

plot(
  cooks_d,
  type = "h",
  ylab = "Cook's distance",
  main = "Cook's distance"
)
abline(h = 4 / nrow(data_CO), col = "red", lty = 2)

which(cooks_d > 4 / nrow(data_CO))

# Standard diagnostic plots
par(mfrow = c(2, 2))
plot(model_CO)
par(mfrow = c(1, 1)) # resets
