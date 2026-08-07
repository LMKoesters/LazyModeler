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
    NULL,
    automatic = TRUE,
    data = d,
    lhs = quote(y)
  )) |>
    expect_message() |>
    expect_warning(regexp = "formatted as numeric")

  expect_equal(family, "binomial")
})

test_that("Gaussian distribution is correctly identified", {
  d <- make_tiny_data()

  (family <- check_model_family(
    NULL,
    automatic = TRUE,
    data = d,
    lhs = quote(x1)
  )) |>
    expect_message()

  expect_equal(family, "gaussian")
})

test_that("Quasibinomial distribution is correctly identified", {
  d <- make_tiny_proportions_data()

  check_model_family(
    NULL,
    automatic = TRUE,
    data = d,
    lhs = quote(y)
  ) |>
    expect_error(regex = "following distributions: gaussian and quasibinomial")
})

test_that("Poisson distribution is correctly identified", {
  d <- make_tiny_poisson_data()

  (family <- check_model_family(
    NULL,
    automatic = TRUE,
    data = d,
    lhs = quote(y)
  )) |>
    expect_warning(regexp = "formatted as numeric") |>
    expect_message()

  expect_equal(family, "poisson")
})

test_that("Warning on non-automatic incorrect family", {
  d <- make_tiny_poisson_data()

  (family <- check_model_family(
    "binomial",
    automatic = FALSE,
    data = d,
    lhs = quote(y)
  )) |>
    expect_warning(regexp = "formatted as numeric") |>
    expect_warning(regexp = "does not match response values")

  expect_equal(family, "binomial")
})

test_that("Family as closure is accepted", {
  d <- make_tiny_poisson_data()

  (family <- check_model_family(
    poisson,
    automatic = FALSE,
    data = d,
    lhs = quote(y)
  )) |>
    expect_warning(regexp = "formatted as numeric")

  expect_equal(family, poisson)
})

test_that("Response is transformed in check_model_family", {
  d <- make_tiny_data()

  (family <- check_model_family(
    NULL,
    automatic = TRUE,
    data = d,
    lhs = quote(I(x1^2))
  )) |>
    expect_message(regex = "Transformed response") |>
    expect_message(regex = "Gaussian is suggested")

  expect_equal(family, "gaussian")
})

test_that("Binomial cbind is accepted", {
  d <- make_grouped_data()


  (family <- check_model_family(
    "binomial",
    automatic = FALSE,
    data = d,
    lhs = quote(cbind(success, failure))
  )) |>
    expect_no_error() |>
    expect_warning(regex = "formatted as numeric") |>
    expect_warning(regex = "formatted as numeric")

  expect_equal(family, "binomial")
})
