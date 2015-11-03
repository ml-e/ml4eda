sim <- function(rho, disturbance, n) {
  stopifnot(all.equal(sum(as.vector(disturbance)), 0.))
  disturbance <- disturbance / sqrt(n)
  dim(disturbance) <- c(2, 2)
  
  x1 <- rbinom(n, 1, 0.5)
  p2_cond_x1 <- 0.5 * (1 + (2*x1 - 1) * rho)
  x2 <- rbinom(n, 1, p2_cond_x1)
  
  treatment_probs <- 0.5 + mapply(function(a, b) disturbance[a + 1, b + 1], x1, x2)
  treatment <- rbinom(n, 1, treatment_probs)
  df <- data.frame(treatment=treatment, x1=x1, x2=x2)
}


sim_null <- function(rho, n=200) {
  sim(rho, array(0, dim=c(2, 2)), n)
}


sim_A <- function(rho, eta, n=200) {
  disturbance <- eta * c(-1, 0, -1, 0)
  sim(rho, disturbance, n)
}


sim_B <- function(rho, eta, n=200) {
  disturbance <- eta * c(-1, 0, 0, 1)
  sim(rho, disturbance, n)
}


sim_C <- function(rho, eta, n=200) {
  disturbance <- eta * c(-1, 1, 1, -1)
  sim(rho, disturbance, n)
}


test_means <- function(df) {
  grouped <- group_by(df, treatment)
  means <- summarize(
    grouped,
    mean_x1=mean(x1),
    mean_x2=mean(x2),
    se_x1=var(x1) / n(),
    se_x2=var(x2) / n())
  means <- mutate(
    means,
    mean_x1=(2*treatment - 1) * mean_x1,
    mean_x2=(2*treatment - 1) * mean_x2)
  means <- ungroup(means)
  ts <- summarize(
    means,
    x1=sum(mean_x1) / sqrt(sum(se_x1)),
    x2=sum(mean_x2) / sqrt(sum(se_x2))
  )
}