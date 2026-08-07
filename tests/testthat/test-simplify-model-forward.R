test_that("glm is optimized (forward selection)", {
  d <- make_significant_factors_data()

  (m <- simplify_model(
    formula = y ~ x1 + I(x1^2) + f2 + f1:x1 + x1 * x2 + x1:x3,
    base_formula = y ~ 1,
    data = d,
    model_type = "glm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward",
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

test_that("lm is optimized (forward selection)", {
  d <- make_tiny_data()

  (m <- simplify_model(
    formula = y ~ x1 + x2 + x3 + f1 + I(x1^2),
    base_formula = y ~ 1,
    data = d,
    model_type = "lm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward"
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

test_that("glmer is optimized (forward selection)", {
  d <- make_lmer_data()

  (m <- simplify_model(
    formula = y ~ x1 + x2 + x3 + x1:x3 + (1 | grp),
    base_formula = y ~ (1 | grp),
    data = d,
    model_type = "glmer",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward",
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
    base_formula = y ~ (1 | grp),
    data = d,
    model_type = "lmer",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward"
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

test_that("gam is optimized (forward simplification)", {
  d <- make_gam_data()

  (m <- simplify_model(
    formula = y ~ s(x1) + x2 + x3 + f1,
    base_formula = y ~ 1,
    data = d,
    model_type = "gam",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward",
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

test_that("Shorthand y ~ . is accepted (forward)", {
  d <- make_significant_factors_data()

  (m <- simplify_model(
    formula = y ~ .,
    base_formula = y ~ 1,
    data = d,
    model_type = "glm",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic", "anova"),
    direction = "forward",
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
