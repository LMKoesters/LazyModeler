test_that("remove_autocorrelations returns expected structure", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 * x2 + (x3 + x4 + x5)^3,
    data = d,
    cols = c("x1", "x2", "x3", "x4", "x5"),
    remove = TRUE,
    threshold = 0.8
  )

  # check output type
  expect_type(res, "list")
  # check output structure
  expect_true(all(c("removed_predictors",
                    "autocorrelations_info",
                    "formula") %in% names(res)))
  # check output type
  expect_true(is.character(res$removed_predictors))
  # check output type
  expect_true(is.data.frame(res$autocorrelations))
  # check autocorrelations columns
  expect_true(all(c("coefficientA",
                    "coefficientB",
                    "correlation",
                    "p_value",
                    "note") %in% names(res$autocorrelations_info)))
})

test_that("no autocorrelations -> removed_predictors is empty", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x5 + x6,
    data = d,
    cols = c("x1", "x5", "x6"),
    remove = TRUE,
    threshold = 0.8
  )

  expect_null(res$autocorrelations_info)
  expect_null(res$removed_predictors)
})

test_that("handle_autocorrelations recognizes priority order of variables", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x2 + x1 + x3,
    data = d,
    cols = c("x2", "x1", "x3"),
    remove = TRUE,
    threshold = 0.8
  )

  # check priority: x2 before x1 means that x1 should be removed
  expect_true("x1" %in% res$removed_predictors)
  expect_false("x2" %in% res$removed_predictors)
})

test_that("handle_autocorrelations recognizes A-B-C constellation", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x2 + x3 + x4 + x5,
    data = d,
    cols = c("x1", "x2", "x3", "x4", "x5"),
    remove = TRUE,
    threshold = 0.8
  )

  # check that x2 was removed (a-b-c test)
  expect_true("x2" %in% res$removed_predictors)
})

test_that("invalid inputs throw informative errors", {
  d <- make_autocor_data()

  expect_error(
    handle_autocorrelations(formula = y ~ x1, data = d,
                            cols = c("does_not_exist")),
    regexp = "undefined columns selected"
  )

  expect_error(
    handle_autocorrelations(formula = y ~ x1 + x2 + x3, data = d,
                            cols = c("x1", "x2", "x3"), threshold = 1.5),
    regexp = "outside of range \\(0,1\\)"
  )
})

test_that("error if autocorrelations but no automatic removal is requested", {
  d <- make_tiny_data()

  expect_warning(
    handle_autocorrelations(
      formula = y ~ x1 + x2 + x3,
      data = d,
      cols = c("x1", "x2", "x3"),
      remove = FALSE,
      threshold = 0.8,
      cor_args = list(method = "spearman",
                      use = "complete.obs")
    ),
    regexp = "Some of your variables are autocorrelated"
  )
})

test_that("formula related columns are recognized", {
  d <- make_significant_factors_data()

  (cols <- formula_related_cols(
    y ~ x1 + x2 + f1 + I(x1^2) + x3:x2,
    d
  )) |>
    expect_no_error()

  expect_equal(c("x2", "x3", "x1"), cols)
})
