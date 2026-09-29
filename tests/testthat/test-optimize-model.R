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

test_that("optimize_model (glm) returns expected structure", {
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
  expect_named(res, c("autocorrelation_result", "backward", "forward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model", "plots")
  )
  expect_named(
    res$forward,
    c("model_selection_result", "final_model", "plots")
  )
  expect_named(
    res$backward$plots,
    c("quality_check",
      "estimates",
      "effect_sizes",
      "categorical_variables",
      "relationships")
  )
  expect_named(
    res$forward$plots,
    c("quality_check",
      "estimates",
      "effect_sizes",
      "categorical_variables",
      "relationships")
  )
})

test_that("optimize_model + lm", {
  d <- make_significant_factors_data()

  res <- optimize_model(
    formula = y ~ I(x1^2) + f2 + f1 * x1 + x1 * x2 + x1 * x3,
    data = d,
    model_type = "lm",
    model_args = list(),
    evaluation_methods = c("anova"),
    directions = c("backward")
  )

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model", "plots")
  )
  expect_named(
    res$backward$plots,
    c("quality_check",
      "estimates",
      "effect_sizes",
      "categorical_variables",
      "relationships")
  )
})

test_that("optimize_model + glmer", {
  d <- make_lmer_data()

  res <- optimize_model(
    formula = y ~ x1 + x2 + x3 + x1:x3 + (1 | grp),
    data = d,
    model_type = "glmer",
    model_args = list(),
    evaluation_methods = c("anova"),
    directions = c("backward"),
    family = poisson
  )

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model", "plots")
  )
  expect_named(
    res$backward$plots,
    c("quality_check",
      "estimates",
      "effect_sizes",
      "categorical_variables",
      "relationships")
  )
})

test_that("optimize_model + lmer", {
  d <- make_lmer_data()

  res <- optimize_model(
    formula = y ~ x1 + x2 + x3 + x1:x3 + (1 | grp),
    data = d,
    model_type = "lmer",
    model_args = list(),
    evaluation_methods = c("anova"),
    directions = c("backward")
  )

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model", "plots")
  )
  expect_named(
    res$backward$plots,
    c("quality_check",
      "estimates",
      "effect_sizes",
      "categorical_variables",
      "relationships")
  )
})

test_that("optimize_model + gam", {
  d <- make_gam_data(n = 1000)

  (res <- optimize_model(
    formula = y ~ s(x1) + x2 + x3 + f1 + f1:x3,
    data = d,
    model_type = "gam",
    model_args = list(),
    evaluation_methods = c("aic", "aicc", "bic"),
    directions = c("backward"),
    family = gaussian
  )) |>
    expect_warning(regexp = "plotting for gam models will be added")

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model")
  )
})

test_that("optimize_model + nls", {
  d <- make_nls_data()
  start <- c(Asym = 5, k = 0.6, offset = .1, slope = .5)

  (res <- optimize_model(
    formula = y ~ offset + Asym * (1 - exp(-k * x)) + slope * x,
    data = d,
    model_type = "nls",
    model_args = list(
      start = start
    ),
    evaluation_methods = c("anova"),
    directions = c("backward"),
    scale_predictors = FALSE
  )) |> expect_warning(regexp = "do not allow for plotting of your model type")

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model")
  )
})

test_that("optimize_model + nlme", {
  d <- make_nlme_data()
  start <- c(Asym = 9, k = 0.7)

  (res <- optimize_model(
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
    expect_warning(regexp = "do not allow for plotting of your model type")

  expect_type(res, "list")
  expect_named(res, c("backward"))
  expect_named(
    res$backward,
    c("model_selection_result", "final_model")
  )
})
