test_that("optimize_model returns expected top-level structure", {
  d <- make_tiny_data()
  
  res <- optimize_model(
    df = d,
    term = quote(y ~ x1 + x2 + x3 + f1),
    autocorrelation_cols = c("x1", "x2", "x3"),
    automatic_removal = TRUE,
    autocorrelation_threshold = 0.8,
    correlation_method = "spearman",
    model_type = "lm",
    evaluation_methods = c("aic"),
    simplification_direction = "backward",
    omit.na = "overall",
    scale_predictor = FALSE,
    plot_quality_assessment = "baseR",
    plot_relationships = FALSE,
    trace = FALSE
  )
  
  expect_type(res, "list")
  expect_true(any(grepl("model|final|summary|plots|autocorr", names(res), ignore.case = TRUE)))
})

test_that("optimize_model updates formula when autocorrelations removed", {
  d <- make_tiny_data()

  res <- optimize_model(
    df = d,
    term = quote(y ~ x1 + x2 + x3),
    autocorrelation_cols = c("x1", "x2", "x3"),
    automatic_removal = TRUE,
    autocorrelation_threshold = 0.8,
    model_type = "lm",
    evaluation_methods = c("aic"),
    simplification_direction = "backward",
    trace = FALSE,
    plot_relationships = FALSE
  )

  # x2 should have been removed (correlated to x1)
  expect_true("removed x2" %in% res$autocorrelations$note)

  m <- res$models_with_info$backward$final_model
  if (!is.null(m)) {
    ftxt <- paste(deparse(stats::formula(m)), collapse = " ")
    expect_false(grepl("\\bx2\\b", ftxt))
  }
})

test_that("optimize_model fails loudly on invalid model_type", {
  d <- make_tiny_data()

  expect_error(
    optimize_model(
      d,
      quote(y ~ x1 + x3),
      model_type = "definitely_not_a_model",
      trace = FALSE
    ),
    regexp = "Unknown model type"
  )
})
