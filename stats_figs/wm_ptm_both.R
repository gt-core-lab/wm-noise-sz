rm(list = ls())
graphics.off()

library(tidyverse)
library(car)
library(lmtest)
library(effectsize)
library(boot)

datapath <- "/Users/wpark78/Documents/code/ptm-wm-sz/csv/ptm-wm.csv"
data <- read_csv(datapath)

data <- data %>%
  mutate(
    group = factor(group),
    Na_c = as.numeric(scale(`Na-fix`, center = TRUE, scale = FALSE)),
    Af_c = as.numeric(scale(`Af-fix`, center = TRUE, scale = FALSE))
  )

model_combined <- lm(
  kavgvar ~ group * (Na_c + Af_c),
  data = data
)

summary(model_combined)
confint(model_combined, level = 0.95)
#standardize_parameters(model_combined)

summary(model_combined)$r.squared
summary(model_combined)$adj.r.squared

# Bootstrap 95% CI for R2
r2_boot_fun_combined <- function(data, indices) {
  d <- data[indices, ]
  m <- lm(kavgvar ~ group * (Na_c + Af_c), data = d)
  summary(m)$r.squared
}

set.seed(123)
boot_r2_combined <- boot(data = data, statistic = r2_boot_fun_combined, R = 5000)
boot.ci(boot_r2_combined, type = c("perc", "bca"))

# Assumption checks
shapiro.test(resid(model_combined))

qqnorm(resid(model_combined))
qqline(resid(model_combined))

plot(fitted(model_combined), resid(model_combined),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs fitted")
abline(h = 0, lty = 2)

bptest(model_combined)

vif(model_combined)

cooks_d <- cooks.distance(model_combined)
plot(cooks_d, type = "h",
     ylab = "Cook's distance",
     main = "Cook's distance")
abline(h = 4 / nrow(data), col = "red", lty = 2)

which(cooks_d > 4 / nrow(data))

par(mfrow = c(2, 2))
plot(model_combined)
par(mfrow = c(1, 1))