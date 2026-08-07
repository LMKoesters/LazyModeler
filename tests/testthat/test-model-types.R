test_that("optimize_model() fails on invalid model_type", {
  d <- make_tiny_data()

  expect_error(
    optimize_model(
      formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
      data = d,
      model_type = "some unknown model",
      family = gaussian
    ),
    regexp = "some unknown model is not supported by LazyModeler"
  )
})

test_that("simple glm computes", {
  d <- make_tiny_data()

  expect_no_error(
    res <- create_model(
      formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
      data = d,
      model_type = "glm",
      family = gaussian
    )
  )

  expect_s3_class(res, "glm")
  expect_named(stats::coef(res), c("(Intercept)",
                                   "x1",
                                   "x2",
                                   "x3",
                                   "f1B",
                                   "f1C",
                                   "I(x1^2)"))
})

test_that("simple lm computes", {
  d <- make_tiny_data()

  expect_no_error(
    res <- create_model(
      formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
      data = d,
      model_type = "lm"
    )
  )

  expect_s3_class(res, "lm")
  expect_named(stats::coef(res), c("(Intercept)",
                                   "x1",
                                   "x2",
                                   "x3",
                                   "f1B",
                                   "f1C",
                                   "I(x1^2)"))
})

test_that("simple glmer computes", {
  d <- make_lmer_data()

  expect_no_error(
    res <- create_model(
      formula = y ~ x1 + x2 + x3 + (1 | grp),
      data = d,
      model_type = "glmer",
      family = poisson
    )
  )

  expect_s4_class(res, "glmerMod")
  expect_named(lme4::fixef(res), c("(Intercept)", "x1", "x2", "x3"))
  expect_named(stats::coef(res), "grp")
})

test_that("simple lmer computes", {
  d <- make_lmer_data()

  expect_no_error({
    res <- create_model(
      formula = y ~ x1 + x2 + x1:x2 + (1 | grp),
      data = d,
      model_type = "lmer"
    )
  })

  expect_s4_class(res, "lmerMod")
  expect_named(lme4::fixef(res), c("(Intercept)", "x1", "x2", "x1:x2"))
  expect_named(stats::coef(res), "grp")
})

test_that("simple gam computes", {
  d <- make_gam_data()

  expect_no_error({
    res <- create_model(
      formula = y ~ s(x1) + x2 + x3,
      data = d,
      model_type = "gam",
      family = gaussian
    )
  })

  expect_s3_class(res, "gam")
})

test_that("simple nlme computes", {
  d <- make_nlme_data()
  start <- c(Asym = 9, k = 0.7)

  expect_no_error({
    res <- create_model(
      formula = y ~ Asym * exp(-k * t),
      data = d,
      model_type = "nlme",
      model_args = list(
        start = start,
        random = quote(Asym ~ 1 | grp),
        fixed = quote(Asym + k ~ 1)
      )
    )
  })

  expect_s3_class(res, c("nlme", "lme"))
  expect_named(stats::coef(res), c("Asym", "k"))
})

test_that("simple nls computes", {
  d <- make_nls_data()
  start <- c(Asym = 5, k = 0.6)

  expect_no_error({
    res <- create_model(
      formula = y ~ Asym * (1 - exp(-k * x)),
      data = d,
      model_type = "nls",
      model_args = list(
        start = start
      )
    )
  })

  expect_s3_class(res, "nls")
  expect_named(stats::coef(res), c("Asym", "k"))
})
