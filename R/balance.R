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


sim_A <- function(rho, eta=1, n=200) {
  disturbance <- eta * c(-1, -1, 1, 1)
  sim(rho, disturbance, n)
}


sim_B <- function(rho, eta=1, n=200) {
  disturbance <- eta * c(-1, 0, 0, 1)
  sim(rho, disturbance, n)
}


sim_C <- function(rho, eta=1, n=200) {
  disturbance <- eta * c(-1, 1, 1, -1)
  sim(rho, disturbance, n)
}

SIMS <- list(null=sim_null, a=sim_A, b=sim_B, c=sim_C)


get_tstats <- function(df) {
  df <- arrange(df, treatment)
  m <- as.matrix(df)
  n <- nrow(m)
  idx <- n - sum(m[, 1])
  df0 <- m[1:idx, -1]
  df1 <- m[(idx+1):n, -1]
  map(1:ncol(df[, -1]), ~ t.test(df0[, .], df1[, .]))
}


test_means <- function(df) {
  tstats <- get_tstats(df)
  pvals <- map(tstats, ~.$p.value)
  min(p.adjust(pvals, method='bon'))
}


test_ols <- function(df) {
  model <- lm(treatment ~ x1 + x2, df)
  fstat <- summary(model)$fstatistic
  pf(
    fstat[1L], 
    fstat[2L],
    fstat[3L],
    lower.tail = FALSE)
}


test_permutation <- function(df, permutations=1000) {
  sum_sq <- function(df) {
    tstats <- get_tstats(df)
    ts <- map_dbl(tstats, ~.$statistic)
    mean(ts^2)
  }

  observed <- sum_sq(df)
  samples <- sapply(1:permutations, function(i) {
    df_perm <- mutate(df, treatment=sample(treatment))
    sum_sq(df_perm)
  })
  percentile <- ecdf(samples)
  pval <- percentile(observed)
  min(pval, (1-pval)) / 2
} 


cross_val <- function(df, method, k=5) {
 assignments <- sample(rep(1:k, length.out=nrow(df)))
 out <- lapply(1:k, function(i) {
   is_test <- assignments == i
   test <- df[is_test, ]
   train <- df[!is_test, ]
   model <- method(train)
   list(true=test$treatment, predicted=predict(model, test))
 })
 list(true=map(out, ~.$true), predicted=map(out, ~.$predicted))
}

test_prediction <- function(df, permutations=1000, test_frac=0.2, method=ols) {
  loss <- function(predicted, true) sum(abs((predicted > .5) - true))
  out <- cross_val(df, method)
  predictions <- out$predicted
  predictions <- flatten(predictions)
  true <- out$true
  observed_loss <- loss(predictions, flatten(true))
  
  loss_samples <- sapply(1:permutations, function(i) {
    resample <- flatmap(true, sample)
    loss(predictions, resample)
  })
  percentile <- ecdf(loss_samples)
  percentile(observed_loss)
}

ols <- function(df) {
  lm(treatment ~ x1 + x2, df)
}

test_ols_prediction <- function(df) {
  test_prediction(df, method=ols)
}

ols_ix <- function(df) {
  lm(treatment ~ x1 + x2 + x1 * x2, df)
}

test_ols_ix_prediction <- function(df) {
  test_prediction(df, method=ols_ix)
}
  
tree <- function(df) {
  rpart::rpart(treatment ~ x1 + x2, df)
}

test_tree_prediction <- function(df) {
  test_prediction(df, method=tree)
}

TESTS <- list(
  means=test_means,
  # permutation=test_permutation,
  ols=test_ols,
  ols_prediction=test_ols_prediction,
  ols_ix_prediction=test_ols_ix_prediction,
  tree_prediction=test_tree_prediction
)
RHOS <- c(0., 0.5)


run_one <- function(sim_name, test_name, rho, n, count=100) {
  sim_strategy <- SIMS[[sim_name]]
  test_strategy <- TESTS[[test_name]]
  results <- sapply(1:count, function(i) {
    ds <- sim_strategy(rho=rho, n=n)
    test_strategy(ds)
  })
  data.frame(sim=sim_name, test=test_name, rho=rho, n=n, result=results)
}


run_all <- function(sim_names=names(SIMS), test_names=names(TESTS), rhos=RHOS, ns=200, count=100) {
  grid <- expand.grid(sim_names=sim_names, test_names=test_names, rhos=rhos, ns=ns)
  results <- mcmapply(
    run_one,
    grid$sim_names,
    grid$test_names,
    grid$rhos,
    grid$ns,
    count=count,
    SIMPLIFY=FALSE,
    mc.cores=4
  )
  bind_rows(results)
}

summarize_output <- function(output) {
  g <- group_by(output, sim, test, rho, n)
  summarize(g, rejection_pct=sum(result<.05)/n())
}

suff <- function(base, suffix) paste(base, suffix, '.csv', sep='')

save_output <- function(output, suffix='') {
  write.csv(output, suff('results', suffix))
  summarized <- summarize_output(output)
  write.csv(summarized, suff('summary', suffix))
}