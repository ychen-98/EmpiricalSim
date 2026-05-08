set.seed(1)
n <- 120
obs_time <- rweibull(n, shape = 1.2, scale = 15)
status   <- rbinom(n, 1, prob = 0.75)
hist_df  <- data.frame(
  time = obs_time,
  age  = rnorm(n, 60, 10),
  sex  = rbinom(n, 1, 0.45),
  ecog = sample(0:2, n, replace = TRUE, prob = c(0.4, 0.4, 0.2))
)
