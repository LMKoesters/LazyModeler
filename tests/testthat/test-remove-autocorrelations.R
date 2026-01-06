test_that("remove_autocorrelations returns expected structure", {
  d <- make_tiny_data()

  res <- remove_autocorrelations(
    df = d,
    coefficients = c("x1", "x2", "x3"),
    use = "complete.obs",
    autocorrelation_threshold = 0.8,
    correlation_method = "spearman",
    automatic_removal = TRUE
  )

  # check output type
  expect_type(res, "list")
  # check output structure
  expect_true(all(c("removed_predictor_cols", "autocorrelations") %in% names(res)))
  # check output type
  expect_true(is.character(res$removed_predictor_cols) || is.null(res$removed_predictor_cols))
  # check output type
  expect_true(is.data.frame(res$autocorrelations))
})

test_that("no autocorrelations -> removed_predictors empty", {
  d <- make_tiny_data()
  # remove correlations from data
  d$x2 <- rnorm(nrow(d))

  res <- remove_autocorrelations(
    df = d,
    coefficients = c("x1", "x2", "x3"),
    autocorrelation_threshold = 0.8,
    automatic_removal = TRUE
  )

  expect_true(length(res$removed_predictor_cols) == 0)
})

test_that("automatic removal removes later variable according to priority order", {
  d <- make_tiny_data()

  res <- remove_autocorrelations(
    df = d,
    coefficients = c("x1", "x2", "x3"),
    autocorrelation_threshold = 0.8,
    automatic_removal = TRUE
  )

  # check priority: x1 before x2 means that x2 should be removed
  expect_true("x2" %in% res$removed_predictor_cols)
  expect_false("x1" %in% res$removed_predictor_cols)
})

test_that("invalid inputs throw informative errors", {
  d <- make_tiny_data()

  expect_error(
    remove_autocorrelations(df=d, coefficients = c("does_not_exist")),
    regexp = "does_not_exist|not found|unknown|column"
  )

  expect_error(
    remove_autocorrelations(df=d, coefficients = c("x1", "x2", "x3"), autocorrelation_threshold = 1.5),
    regexp = "threshold|between|0|1"
  )
})

test_that("error if autocorrelations but no automatic removal is requested", {
  d <- make_tiny_data()

  expect_error(
    remove_autocorrelations(
      df = d,
      coefficients = c("x1", "x2", "x3"),
      use = "complete.obs",
      autocorrelation_threshold = 0.8,
      correlation_method = "spearman",
      automatic_removal = FALSE
    ),
    regexp = "Some of your variables are autocorrelated"
  )
})

test_that("check correct A == B; B == C; A != C removal", {
  d <- make_tiny_a_b_c_data()

  res <- remove_autocorrelations(
    df = d,
    coefficients = c("xA", "xB", "xC"),
    use = "complete.obs",
    autocorrelation_threshold = 0.7,
    correlation_method = "pearson",
    automatic_removal = TRUE
  )$removed_predictor_cols

  expect_true("xB" %in% res)
})

test_that("check correct A == B; B == C; A == C removal", {
  d <- make_tiny_a_b_c_equal_data()

  res <- remove_autocorrelations(
    df = d,
    coefficients = c("xA", "xB", "xC"),
    use = "complete.obs",
    autocorrelation_threshold = 0.7,
    correlation_method = "pearson",
    automatic_removal = TRUE
  )$removed_predictor_cols

  expect_true("xC" %in% res)
})
