test_that("plotting of formula with 0 interactions successful", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + f2,
    data = d,
    model_type = "glm",
    model_args = list()
  )
  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR"))
  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
})

test_that("plotting of formula with 0 categorical variables successful", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + x1:x2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR"))
  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
  expect_all_true(
    vapply(p$categorical_variables,
           function(x) {
             length(x) == 0
           },
           logical(1))
  )
})

test_that("plotting of formula with 0 numerical variables successful", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ f1 + f2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR"))

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
  expect_true(length(p$relationships) == 0)
})

test_that("plotting of interactions between >2 vars throws warning", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + f2 + x1:f2 + x1:x2:f2 + x1:x2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR")) |>
    expect_warning(regexp = "interactions with more than two variables")
})

test_that("plotting of interactions between factors throws warning", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + f2 + x1:f2 + f1:f2 + x1:x2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR")) |>
    expect_warning(regexp = "interactions between two factors")
})

test_that("plot output formatted correctly", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + f2 + x1:f2 + x1:x2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "baseR"))

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
  expect_equal(names(p$categorical_variables), c("plots", "stat_results"))
  expect_equal(names(p$categorical_variables$plots), c("f2"))
  expect_equal(names(p$categorical_variables$stat_results), c("f2"))
  expect_equal(names(p$relationships), c("x1", "I(x1^2)", "x1:x2", "x1:f2"))
  expect_all_true(
    vapply(p$categorical_variables$plots, inherits, logical(1), what = "ggplot")
  )
  expect_all_true(
    vapply(p$relationships, inherits, logical(1), what = "ggplot")
  )
})

test_that("glm with performance successful", {
  d <- make_significant_factors_data()

  final_model <- create_model(
    formula = y ~ x1 + I(x1^2) + f2 + x1:f2 + x1:x2,
    data = d,
    model_type = "glm",
    model_args = list()
  )
  (p <- plot_model(model = final_model,
                   model_type = "glm",
                   quality_assessment = "performance")) |>
    expect_no_error()

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
})

test_that("gam plots successfully", {
  d <- make_gam_data()

  final_model <- create_model(
    formula = y ~ s(x1) + x2 + x3 + f1,
    data = d,
    model_type = "gam",
    model_args = list(),
    family = gaussian
  )

  (p <- plot_model(model = final_model,
                   model_type = "gam",
                   quality_assessment = "baseR")) |>
    expect_no_error()

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
})

test_that("glmer plots successfully", {
  d <- make_lmer_data()

  final_model <- create_model(
    formula = y ~ x1 + x2 + x3 + (1 | grp),
    data = d,
    model_type = "glmer",
    family = poisson
  )

  (p <- plot_model(model = final_model,
                   model_type = "glmer",
                   quality_assessment = "baseR")) |>
    expect_no_error()

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
})

test_that("lmer plots successfully", {
  d <- make_lmer_data()

  final_model <- create_model(
    formula = y ~ x1 + x2 + x3 + (1 | grp),
    data = d,
    model_type = "lmer"
  )

  (p <- plot_model(model = final_model,
                   model_type = "lmer",
                   quality_assessment = "baseR")) |>
    expect_no_error()

  expect_equal(names(p), c("quality_check", "estimates", "effect_sizes",
                           "categorical_variables", "relationships"))
})
