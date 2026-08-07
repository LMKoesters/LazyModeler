test_that("family = gaussian() -> glm computes", {
  d <- make_tiny_data()

  (res <- create_model(
    formula = count ~ x1 + x3,
    data = d,
    model_type = "glm",
    family = gaussian()
  )) |>
    expect_no_error()

  expect_identical(as.character(stats::family(res)$family), "gaussian")
})

test_that("family = gaussian -> glm computes", {
  d <- make_tiny_data()

  (res <- create_model(
    formula = quote(count ~ x1 + x3),
    data = d,
    model_type = "glm",
    family = gaussian
  )) |>
    expect_no_error()

  expect_identical(as.character(res$family$family), "gaussian")
})

test_that("weights = quote(trials) -> glm computes", {
  d <- make_tiny_data()

  (res <- create_model(
    formula = prop ~ x1 + x3,
    data = d,
    model_type = "glm",
    family = "quasibinomial",
    model_args = list(weights = quote(trials))
  )) |>
    expect_no_error()

  res_frame <- model.frame(res)
  expect_true(any(stats::model.weights(res_frame) != 1))
})

test_that("offset = quote(log(exposure)) -> glm computes", {
  d <- make_tiny_data()

  (res <- create_model(
    formula = count ~ x1 + x3,
    data = d,
    model_type = "glm",
    family = "poisson",
    model_args = list(offset = quote(log(exposure)))
  )) |>
    expect_no_error()

  res_frame <- model.frame(res)
  expect_false(is.null(stats::model.offset(res_frame)))
})
