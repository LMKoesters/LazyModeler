test_that("simplify_model fits a model and returns expected components", {
  d <- make_tiny_data()

  res <- simplify_model(
    df = d,
    model_type = "lm",
    term = y ~ x1 + x2 + x3 + f1,
    evaluation_methods = c("aic", "aicc", "bic"),
    direction = "both",
    categorical_vars=c('f1'),
    backward_simplify_model = TRUE,
    trace = TRUE,
    omit.na = "overall"
  )

  expect_type(res, "list")
  expect_type(res$backward, "list")
  expect_type(res$forward, "list")
  
  expect_true(any(c("final_model", "significant_variables", "marginally_significant_variables", "history") %in% names(res$backward)))
  expect_true(any(c("final_model", "significant_variables", "marginally_significant_variables", "history") %in% names(res$forward)))
})

test_that("simplification removes non-informative terms", {
  d <- make_tiny_data()
  # function removes x2 as it correlates to x1
  res <- simplify_model(
    df = d,
    model_type = "lm",
    term = y ~ x1 + x2 + x3 + f1,
    evaluation_methods = c("aic"),
    direction = "backward",
    categorical_vars=c('f1'),
    backward_simplify_model = TRUE,
    trace = FALSE
  )$backward

  m <- res$final_model
  expect_s3_class(m, "lm")
  expect_false(grepl("\\bx2\\b", paste(deparse(stats::formula(m)), collapse = " ")))
})

test_that("interaction terms keep main effects", {
  d <- make_interaction_data()

  res <- simplify_model(
    df = d,
    model_type = "lm",
    term = y ~ x1 * x3 + x2,
    evaluation_methods = c("aic"),
    direction = "backward",
    backward_simplify_model = FALSE,
    trace = FALSE
  )$backward

  m <- res$final_model
  ftxt <- paste(deparse(stats::formula(m)), collapse = " ")

  expect_true(grepl("x1 \\* x3", ftxt))
  # x1 and x3 (as main effects) should not be removed
  expect_true(grepl("\\bx1\\b", ftxt))
  expect_true(grepl("\\bx3\\b", ftxt))
})

test_that("interaction terms keep main effects", {
  d <- make_interaction_data()

  res <- simplify_model(
    df = d,
    model_type = "lm",
    term = y ~ x1:x3 + x2,
    evaluation_methods = c("aic"),
    direction = "backward",
    backward_simplify_model = FALSE,
    trace = FALSE
  )$backward

  m <- res$final_model
  ftxt <- paste(deparse(stats::formula(m)), collapse = " ")

  expect_true(grepl("x1:x3", ftxt))
  # x1 and x3 (as main effects) should not be removed
  expect_true(grepl("\\bx1\\b", ftxt))
  expect_true(grepl("\\bx3\\b", ftxt))
})

test_that("NA handling works", {
  d <- make_tiny_data_with_na()

  expect_error(
    simplify_model(
      df = d,
      model_type = "lm",
      term = y ~ x1 + x3,
      evaluation_methods = c("aic"),
      direction = "backward",
      backward_simplify_model = TRUE,
      trace = FALSE,
      omit.na = "none"
    ),
    regexp = "Unknown|option|omit"
  )

  expect_no_error(
    simplify_model(
      df = d,
      model_type = "lm",
      term = y ~ x1 + x3,
      evaluation_methods = c("aic"),
      direction = "backward",
      backward_simplify_model = TRUE,
      trace = FALSE,
      omit.na = "overall"
    )
  )
})

test_that("simplify_model works for glm (binomial) and drops a noise term", {
  d <- make_tiny_glm_data()

  res <- simplify_model(
    df = d,
    model_type = "glm",
    term = y ~ x1 + x2 + f1,
    evaluation_methods = c("aic"),
    direction = "backward",
    categorical_vars = "f1",
    backward_simplify_model = TRUE,
    trace = FALSE
  )$backward

  m <- res$final_model %||% res$model %||% res$fit

  expect_s3_class(m, "glm")
  expect_identical(stats::family(m)$family, "binomial")

  # check term
  tt <- attr(stats::terms(m), "term.labels")

  # signal should remain
  expect_true("x1" %in% tt)
  expect_true("f1" %in% tt)

  # noise should be removed
  expect_false("x2" %in% tt)
})

test_that("simplify_model runs for nlme and returns an nlme object", {
  d <- make_tiny_nlme_data()

  # nlmer needs starting values
  start <- c(Asym = 9, k = 0.7)

  res <- simplify_model(
    df = d,
    model_type = "nlme",
    term = quote(y ~ Asym * exp(-k * t)),
    start = start,
    non_linear = quote(y ~ Asym * exp(-k * t)),
    random = quote(Asym ~ 1 | grp),
    fixed = Asym + k ~ 1,
    evaluation_methods = c("aic"),
    direction = "backward",
    backward_simplify_model = TRUE,
    trace = TRUE
  )$nlme

  m <- res$final_model

  expect_s3_class(m, c("nlme", "lme"))

  # basic sanity: fitted values exist
  fv <- stats::fitted(m)
  expect_true(is.numeric(fv))
  expect_equal(length(fv), nrow(d))

  ll <- stats::logLik(m)
  expect_true(is.finite(as.numeric(ll)))
})

test_that("simplify_model works for lmer and returns a mixed model", {
  d <- make_tiny_lmer_data()

  res <- simplify_model(
    df = d,
    model_type = "lmer",
    term = y ~ x1 + x2 + x1:x2 + (1 | grp),
    evaluation_methods = "aic",
    direction = "both",
    backward_simplify_model = TRUE,
    trace = FALSE
  )$backward

  m <- res$final_model

  # lmerTest returns lmerModLmerTest; inherits merMod too
  expect_true(inherits(m, c("lmerModLmerTest", "lmerMod", "merMod")))

  # fitted values exist and have right length
  fv <- stats::fitted(m)
  expect_true(is.numeric(fv))
  expect_equal(length(fv), nrow(d))

  # fixed effect x1 should be present
  tt <- attr(stats::terms(m), "term.labels")
  expect_true("x1" %in% tt)

  # noise term should be dropped
  expect_false("x2" %in% tt)
})
