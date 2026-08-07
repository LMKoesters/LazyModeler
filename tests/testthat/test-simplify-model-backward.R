test_that("glm is optimized (backward simplification)", {
  d <- make_significant_factors_data()

  (m <- simplify_model(
    formula = y ~ x1 + I(x1^2) + f2 + f1:x1 + x1 * x2 + x1:x3,
    data = d,
    model_type = "glm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward",
    family = gaussian
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }

    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})

test_that("lm is optimized (backward simplification)", {
  d <- make_tiny_data()

  (m <- simplify_model(
    formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
    data = d,
    model_type = "lm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward"
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }
    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})

test_that("glmer is optimized (backward simplification)", {
  d <- make_lmer_data()

  (m <- simplify_model(
    formula = y ~ x1 + x2 + x3 + x1:x3 + (1 | grp),
    data = d,
    model_type = "glmer",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward",
    family = poisson
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }
    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})

test_that("lmer is optimized (backward simplification)", {
  d <- make_lmer_data()

  (m <- simplify_model(
    formula = y ~ x1 + x2 + x1:x2 + x1:x3 + (1 | grp),
    data = d,
    model_type = "lmer",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward"
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }
    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})

test_that("gam is optimized (backward simplification)", {
  d <- make_gam_data()

  (m <- simplify_model(
    formula = y ~ s(x1) + x2 + x3 + f1,
    data = d,
    model_type = "gam",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward",
    family = gaussian
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }
    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})

test_that("nls is returned as is", {
  d <- make_nls_data()
  start <- c(Asym = 5, k = 0.6, offset = .1, slope = .5)

  (m <- simplify_model(
    formula = y ~ offset + Asym * (1 - exp(-k * x)) + slope * x,
    data = d,
    model_type = "nls",
    model_args = list(
      start = start
    ),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward"
  )) |>
    expect_no_error()

  expect_named(m, c("assessments", "final_model", "p_values"))
  expect_s3_class(m$final_model, "nls")
})

test_that("nlme is returned as is", {
  d <- make_nlme_data()
  start <- c(Asym = 9, k = 0.7)

  (m <- simplify_model(
    formula = y ~ Asym * exp(-k * t),
    data = d,
    model_type = "nlme",
    model_args = list(
      start = start,
      random = quote(Asym ~ 1 | grp),
      fixed = quote(Asym + k ~ 1)
    ),
    evaluation_methods = c("aic", "aicc", "bic", "anova")
  )) |>
    expect_no_error()

  expect_named(m, c("assessments", "final_model", "p_values"))
  expect_s3_class(m$final_model, c("nlme", "lme"))
})

test_that("omit NAs works", {
  d <- make_tiny_data_with_na()

  (m <- simplify_model(
    formula = y ~ x1 + I(x1^2) + x1 * x2 + x1:x3,
    data = d,
    model_type = "glm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward",
    family = gaussian
  )) |>
    expect_no_error()

  d <- stats::model.frame(y ~ x1 + I(x1^2) + x1 * x2 + x1:x3,
                          data = d,
                          na.action = na.omit)

  expect_equal(nrow(d), nrow(stats::model.frame(m$final_model)))
})

test_that("Shorthand y ~ . is accepted (backward)", {
  d <- make_significant_factors_data()

  (m <- simplify_model(
    formula = y ~ .,
    data = d,
    model_type = "glm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "backward",
    family = gaussian
  )) |>
    expect_no_error()

  history <- m$history
  for (i in 1:(length(history) - 1)) {
    m1 <- history[[names(history)[[i]]]]
    m2 <- history[[names(history)[[i + 1]]]]
    if (((length(m2$assessments) == 0) &&
           (i == (length(history) - 1))) ||
          ((i == (length(history) - 1)) &&
             (!identical(m2$model, m$final_model)))) {
      break
    }

    expect_false(identical(m1$model, m2$model))
    for (eval_m in c("aic", "aicc", "bic")) {
      expect_true((m2$assessments[[eval_m]] - (m1$assessments[[eval_m]]) < 2))
    }
  }
})
