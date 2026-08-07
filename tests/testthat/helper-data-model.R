make_tiny_data <- function(n = 80, seed = 42) {
  set.seed(seed)

  x1 <- rnorm(n)
  x2 <- x1 * 0.95 + rnorm(n, sd = 0.05)  # correlated to x1
  x3 <- rnorm(n)
  f1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # linear relation
  y <- 1 + 2 * x1 - 1.5 * x3 + ifelse(f1 == "B", 0.5, 0) + rnorm(n, sd = 0.2)

  d <- data.frame(y = y,
                  x1 = x1,
                  x2 = x2,
                  x3 = x3,
                  f1 = f1,
                  exposure = runif(n, 1, 5),
                  trials = sample(5:20, n, TRUE))

  d$rate <- exp(0.2 + 0.5 * d$x1) * d$exposure
  d$count <- rpois(n, d$rate)
  d$prop <- rbinom(n, d$trials, plogis(-0.2 + 0.5 * d$x1)) / d$trials

  d
}

make_significant_factors_data <- function(seed = 42, n_per_group = 30) {
  set.seed(seed)
  test_data <- tidyr::expand_grid(
    f1 = c("A", "B", "C"),
    f2 = c("shallow", "medium", "deep"),
    observation = seq_len(n_per_group)
  ) |>
    dplyr::mutate(
      f1 = factor(.data$f1, levels = c("A", "B", "C")),
      f2 = factor(.data$f2, levels = c("shallow", "medium", "deep")),
      x1 = stats::rnorm(dplyr::n(), mean = 10, sd = 2),
      x2 = stats::runif(dplyr::n(), min = 0, max = 5),
      x3 = stats::runif(dplyr::n(), min = 0, max = 5),
      f1_effect = dplyr::recode(.data$f1, A = 0, B = 6, C = 12),
      f2_effect = dplyr::recode(.data$f2, shallow = 0, medium = 5, deep = 10),
      y = 20 + .data$x1 + 1.5 * .data$x1 - 2 * .data$x2 +
        .data$f1_effect + .data$f2_effect +
        stats::rnorm(dplyr::n(), mean = 0, sd = 1.5)
    ) |>
    dplyr::select(-c("f1_effect", "f2_effect", "observation"))

  test_data
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

make_grouped_data <- function(n = 80, seed = 42) {
  set.seed(42)

  x <- rnorm(n)
  trials <- sample(5:20, n, TRUE)
  p <- plogis(-0.2 + 0.8 * x)
  success <- rbinom(n, trials, p)
  failure <- trials - success

  data.frame(success = success,
             failure = failure,
             x = x)
}

make_lmer_data <- function(n_groups = 20, n_per_group = 20, seed = 42) {
  set.seed(seed)

  n <- n_groups * n_per_group

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  grp <- factor(rep(seq_len(n_groups),
                    each = n_per_group))
  group_effect <- rnorm(n_groups,
                        mean = 0,
                        sd = 0.6)

  eta <- 0.5 +
    0.4 * x1 -
    0.3 * x2 +
    0.2 * x3 +
    group_effect[grp]

  # poisson response
  y <- rpois(n,
             lambda = exp(eta))

  # binomial response
  p <- plogis(eta)
  y_binom <- rbinom(n, size = 1, prob = p)

  data.frame(
    y = y,
    y_binom = y_binom,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    grp = grp
  )
}

make_nlme_data <- function(n_groups = 8,
                           n_per = 25,
                           k_true = 0.8,
                           seed = 42) {
  set.seed(seed)

  grp <- factor(rep(seq_len(n_groups), each = n_per))
  t <- rep(seq(0, 3, length.out = n_per), times = n_groups)

  # true parameters
  asym_g <- rnorm(n_groups, mean = 10, sd = 1.5)  # random Asym/group

  asym <- asym_g[as.integer(grp)]
  y <- asym * exp(-k_true * t) + rnorm(length(t), sd = 0.25)

  data.frame(y = y,
             t = t,
             grp = grp)
}

make_nls_data <- function(n = 80, seed = 42) {
  set.seed(seed)

  x <- seq(0, 10, length.out = n)

  asym <- 5
  k <- 0.6

  y <- asym * (1 - exp(-k * x)) + rnorm(n, sd = 0.15)

  data.frame(
    y = y,
    x = x
  )
}

make_gam_data <- function(n = 150, seed = 42) {
  set.seed(seed)

  x1 <- runif(n, 0, 10)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  f1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # non-linear effect
  y <- sin(x1) + 0.4 * x2 + rnorm(n, sd = 0.25)

  data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    f1 = f1
  )
}

make_unformatted_binary_data <- function(n = 40, seed = 42) {
  set.seed(seed)

  x <- rnorm(n)
  z <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.2 + 1.2 * x - 0.8 * z))

  data.frame(
    y = y,
    x = x,
    z = z
  )
}

make_tiny_poisson_data <- function(n = 40,
                                   seed = 42,
                                   beta0 = 0.5,
                                   beta1 = 0.4) {
  set.seed(seed)

  x <- runif(n, min = 0, max = 3)
  lambda <- exp(beta0 + beta1 * x)

  data.frame(
    x = x,
    y = rpois(n, lambda = lambda)
  )
}

make_tiny_data_with_na <- function(n = 60, seed = 42) {
  set.seed(seed)
  d <- make_tiny_data(n, seed)
  d$x3[sample.int(n, 5)] <- NA
  d
}

make_tiny_proportions_data <- function(n = 60, seed = 42) {
  set.seed(seed)
  
  x1 <- rnorm(n)
  x2 <- x1 * 0.95 + rnorm(n, sd = 0.05)  # correlated to x1
  x3 <- rnorm(n)
  f1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  
  # linear relation
  y <- 1 + 2 * x1 - 1.5 * x3 + ifelse(f1 == "B", 0.5, 0) + rnorm(n, sd = 0.2)
  y <- (y - min(y)) / (max(y) - min(y))
  
  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, f1 = f1)
}
