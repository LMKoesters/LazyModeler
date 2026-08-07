make_autocor_data <- function(n = 80, seed = 42) {
  set.seed(seed)

  x1 <- rnorm(n)
  x2 <- x1 * 0.9 + rnorm(n) * .7  # correlated to x1
  x3 <- rnorm(n) * .7 + x2 # correlated to x2, but not to x1
  x4 <- rnorm(n) + x3
  x5 <- x4 * .5 + rnorm(n) * 0.1
  x6 <- rnorm(n)

  # linear relation
  y <- 1 + 2 * x1 - 1.5 * x3 + rnorm(n, sd = 0.2)

  d <- data.frame(y = y,
                  x1 = x1,
                  x2 = x2,
                  x3 = x3,
                  x4 = x4,
                  x5 = x5,
                  x6 = x6)

  d
}

make_tiny_a_b_c_data <- function(n = 50000, seed = 42) {
  set.seed(seed)

  x_a <- rnorm(n, sd = 1)
  e1 <- rnorm(n, sd = 1)

  x_b <- e1 + x_a
  x_c <- e1

  y <- 3 * x_a + rnorm(n)

  data.frame(
    y = y,
    xA = x_a,
    xB = x_b,
    xC = x_c
  )
}

make_tiny_a_b_c_equal_data <- function(n = 50000, seed = 42) {
  set.seed(seed)

  x_a <- rnorm(n)
  x_b <- x_a * 2
  x_c <- x_a * 4

  y <- 3 * x_a + rnorm(n)

  data.frame(
    y = y,
    xA = x_a,
    xB = x_b,
    xC = x_c
  )
}
