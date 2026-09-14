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
