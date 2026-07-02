rm(list = ls())
graphics.off()

library(BayesFactor)

# frequentist result
t_value <- 1.65
n1 <- 30   # neurotypical / CO
n2 <- 27   # schizophrenia / SZ

# Bayesian independent-samples t-test from t statistic
bf10 <- ttest.tstat(
  t = t_value,
  n1 = n1,
  n2 = n2,
  rscale = "medium",
  simple = TRUE
)

bf01 <- 1 / bf10

bf10
bf01

library(effectsize)

F_to_eta2(
  f = 1.65^2,
  df = 1,
  df_error = 55,
  ci = 0.95
)

