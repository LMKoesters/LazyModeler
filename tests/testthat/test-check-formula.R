test_that("check formula: y ~ .", {
  d <- make_significant_factors_data()

  (formula <- check_formula(y ~ ., d)) |>
    expect_no_error()

  expect_identical(y ~ f1 + f2 + x1 + x2 + x3, formula)
})

test_that("check main effects are added", {
  d <- make_significant_factors_data()

  (check_formula(y ~ x1 + x2:x3 + x2:f1, d)) |>
    expect_no_error() |>
    expect_warning(regexp = paste("We will add the following main effects to",
                                  "your formula: x2, x3, f1",
                                  collapse = " "))
})
