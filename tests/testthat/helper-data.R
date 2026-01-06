make_tiny_data <- function(n = 60, seed = 42) {
  set.seed(seed)
  
  x1 <- rnorm(n)
  x2 <- x1 * 0.95 + rnorm(n, sd = 0.05)  # correlated to x1
  x3 <- rnorm(n)
  f1 <- factor(sample(c("A","B","C"), n, replace = TRUE))
  
  # linear relation
  y <- 1 + 2*x1 - 1.5*x3 + ifelse(f1 == "B", 0.5, 0) + rnorm(n, sd = 0.2)
  
  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, f1 = f1)
}

make_tiny_categorical_data <- function(n = 60, seed = 42) {
  set.seed(seed)
  
  x1 <- rnorm(n)
  x2 <- x1 * 0.95 + rnorm(n, sd = 0.05)  # correlated to x1
  x3 <- rnorm(n)
  f1 <- factor(sample(c("A","B","C"), n, replace = TRUE))
  f2 <- factor(sample(c("high", "low", "deep", "shallow"), n, replace = TRUE))
  
  # linear relation
  y <- 1 + 2*x1 - 1.5*x3 + ifelse(f1 == "B", 0.5, 0)  + ifelse(f1 == "high", 0.5, 0) + rnorm(n, sd = 0.2)
  
  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, f1 = f1, f2 = f2)
}

make_tiny_data_with_na <- function(n = 60, seed = 42) {
  d <- make_tiny_data(n, seed)
  d$x3[sample.int(n, 5)] <- NA
  d
}

make_tiny_a_b_c_data <- function(n = 50000, seed = 42) {
  set.seed(seed)
  
  xA <- rnorm(n, sd = 1)
  e1 <- rnorm(n, sd = 1)
  e2 <- rnorm(n)

  xB <- e1 + xA # A == B
  xC <- e1 # B == C & A != C
  
  y <- 3 * xA + rnorm(n)
  
  data.frame(
    y = y,
    xA = xA,
    xB = xB,
    xC = xC
  )
}

make_tiny_a_b_c_equal_data <- function(n = 50000, seed = 42) {
  set.seed(seed)
  
  xA <- rnorm(n)
  xB <- xA * 2
  xC <- xA * 4
  
  y <- 3 * xA + rnorm(n)
  
  data.frame(
    y = y,
    xA = xA,
    xB = xB,
    xC = xC
  )
}

make_interaction_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x3 <- rnorm(n)
  x2 <- rnorm(n)
  
  # strong interaction: effect through x1*x3
  y <- 2 + 0.1*x1 + 0.1*x3 + 6*(x1*x3) + rnorm(n, sd = 0.5)
  
  data.frame(y, x1, x2, x3)
}

make_tiny_glm_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f1 <- factor(sample(c("A", "B"), n, replace = TRUE))
  
  # signal through x1 + f1
  eta <- -0.5 + 2.0 * x1 + ifelse(f1 == "B", 0.8, 0)
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n, size = 1, prob = p)
  
  data.frame(y = y, x1 = x1, x2 = x2, f1 = f1)
}

make_tiny_nlme_data <- function(n_groups = 8, n_per = 25, k_true = 0.8, seed = 42) {
  set.seed(seed)
  
  grp <- factor(rep(seq_len(n_groups), each = n_per))
  t <- rep(seq(0, 3, length.out = n_per), times = n_groups)
  
  # true parameters
  Asym_g <- rnorm(n_groups, mean = 10, sd = 1.5)  # random Asym/group
  
  Asym <- Asym_g[as.integer(grp)]
  y <- Asym * exp(-k_true * t) + rnorm(length(t), sd = 0.25)
  
  data.frame(y = y, t = t, grp = grp)
}

make_tiny_lmer_data <- function(n_groups = 12, n_per = 25, seed = 42) {
  set.seed(seed)
  
  grp <- factor(rep(seq_len(n_groups), each = n_per))
  n <- length(grp)
  
  x1 <- rnorm(n)
  x2 <- rnorm(n)  # noise
  # random intercept by group
  b0 <- rnorm(n_groups, sd = 1.0)[as.integer(grp)]
  
  # fixed + random + noise
  y <- 1.0 + 2.0 * x1 + b0 + rnorm(n, sd = 0.6)
  
  data.frame(y = y, x1 = x1, x2 = x2, grp = grp)
}