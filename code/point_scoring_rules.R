ses <- function(obs, pred) (pred - obs)^2
sle <- function(obs, pred) (log(pred + 1) - log(obs + 1))^2
ae <- function(obs, pred) abs(pred - obs)
