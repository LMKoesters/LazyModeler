test_that("Unknown arguments throw error", {
  d <- make_autocor_data()

  optimize_model(
    formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
    data = d,
    model_type = "glm",
    family = gaussian,
    some_unknown_arg = 1
  ) |>
    expect_error(regexp = "Unknown or partially matched")
})

test_that("Partially matching arguments throw error", {
  d <- make_autocor_data()

  optimize_model(
    formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
    data = d,
    model_type = "glm",
    family = gaussian,
    evaluation_method = 1
  ) |>
    expect_error(regexp = "Unknown or partially matched")
})

test_that("optimize_model returns expected structure", {
  d <- make_autocor_data()

  res <- optimize_model(
    formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
    data = d,
    model_type = "glm",
    family = gaussian,
    directions = c("backward", "forward"),
    detect_autocors = TRUE,
    remove_autocors = TRUE,
    base_formula = y ~ 1
  )

  expect_type(res, "list")
  expect_named(res, c("autocorrelation_result", "models_with_info"))
  expect_named(
    res$models_with_info,
    c("backward", "forward")
  )
  # TODO
})
