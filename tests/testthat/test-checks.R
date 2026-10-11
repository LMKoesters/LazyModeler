test_that("Transformed response is handled correctly", {
  n <- 50
  d <- data.frame(y = abs(rnorm(n)),
                  x1 = rnorm(n),
                  logy = factor(sample(c("A", "B", "C"), n, replace = TRUE)))

  (family <- check_model_family(
    log(y) ~ x1,
    "glm",
    NULL,
    automatic = TRUE,
    data = d
  )) |>
    expect_message("Gaussian is suggested")

  expect_equal(family, "gaussian")
})

test_that("gam is handled correctly", {
  d <- make_gam_data()

  (family <- check_model_family(
    y ~ s(x1) + x2 + x3,
    "gam",
    NULL,
    automatic = TRUE,
    data = d
  )) |>
    expect_message("Gaussian is suggested")

  expect_equal(family, "gaussian")
})

test_that("Unsupported model type is rejected (string)", {
  expect_error(
    check_model_type("some model"),
    regex = "not supported by LazyModeler"
  )
})

test_that("Support model type is accepted", {
  expect_no_error(
    check_model_type("glm")
  )
})

test_that("Binomial distribution is correctly identified", {
  d <- make_unformatted_binary_data()

  (family <- check_model_family(
    y ~ x,
    "glm",
    NULL,
    automatic = TRUE,
    data = d
  )) |>
    expect_message("Continuing with binomial distribution")

  expect_equal(family, "binomial")
})

test_that("Gaussian distribution is correctly identified", {
  d <- make_tiny_data()

  (family <- check_model_family(
    y ~ x1,
    "glm",
    NULL,
    automatic = TRUE,
    data = d
  )) |>
    expect_message()

  expect_equal(family, "gaussian")
})

test_that("Quasibinomial distribution is correctly identified", {
  d <- make_tiny_proportions_data()

  check_model_family(
    y ~ x1,
    "glm",
    NULL,
    automatic = TRUE,
    data = d
  ) |>
    expect_error(regex = "following distributions: gaussian and quasibinomial")
})

test_that("Poisson distribution is correctly identified", {
  d <- make_tiny_poisson_data()

  (family <- check_model_family(
    y ~ x,
    "glm",
    NULL,
    automatic = TRUE,
    data = d
  )) |>
    expect_message(regexp = "Poisson assumes non-negative")

  expect_equal(family, "poisson")
})

test_that("Warning on non-automatic incorrect family", {
  d <- make_tiny_poisson_data()

  (family <- check_model_family(
    y ~ x,
    "glm",
    "binomial",
    automatic = FALSE,
    data = d
  )) |>
    expect_warning(
      regexp = "Chosen distribution 'binomial' does not match response values"
    )

  expect_equal(family, "binomial")
})

test_that("Family as closure is accepted", {
  d <- make_tiny_poisson_data()

  family <- check_model_family(
    y ~ x,
    "glm",
    poisson,
    automatic = FALSE,
    data = d
  )

  expect_equal(family, poisson)
})

test_that("Binomial cbind is accepted", {
  d <- make_grouped_data(n = 50)

  (family <- check_model_family(
    cbind(success, failure) ~ x,
    "glm",
    "binomial",
    automatic = FALSE,
    data = d
  ))

  expect_equal(family, "binomial")
})

test_that("Non-binomial cbind is rejected", {
  d <- make_grouped_data(n = 50)
  d$x1 <- rnorm(50)

  (check_model_family(
    cbind(success, x1) ~ x,
    "glm",
    "binomial",
    automatic = TRUE,
    data = d
  )) |>
    expect_error(regexp = "Quitting.")
})
