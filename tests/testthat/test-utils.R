test_that("3-way interactions are extracted", {
  terms <- attr(stats::terms.formula(y ~ f1:f2:x1), "term.labels")
  interactions <- extract_interactions(terms)

  expect_true(length(unique(interactions$main_effect)) == 3)
})

test_that("All columns related to formula are extracted", {
  d <- make_significant_factors_data()

  cols <- formula_related_cols(y ~ x1 + I(x1^2) + f2 + x1:f2 + x1:x2:f2 + x1:x2,
                               d)

  expect_equal(cols, c("x1", "x2"))
})

test_that("3-way interactions are transformed into trait interactions", {
  d <- make_significant_factors_data()

  m <- create_model(
    formula = y ~ x1 + f2 + x3 + (x1:x3:f2)^3,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  model_overview <- stats::coef(summary(m)) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "predictor")

  c(categorical_vars,
    interactions,
    numeric_vars) %<-% get_term_factors(model = m,
                                        model_type = "glm",
                                        model_overview = model_overview)
  expect_length(setdiff(categorical_vars$predictor,
                        unique(paste0("f2", d$f2))), 0)
  expect_length(setdiff(interactions$interaction,
                        unique(paste0("x1:f2", d$f2, ":x3"))), 0)
})

test_that("factor-factor interactions transformed into trait interactions", {
  d <- make_significant_factors_data()

  m <- create_model(
    formula = y ~ x1 + f1:f2,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  model_overview <- stats::coef(summary(m)) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "predictor")

  c(categorical_vars,
    interactions,
    numeric_vars) %<-% get_term_factors(model = m,
                                        model_type = "glm",
                                        model_overview = model_overview)
  combinations <- expand.grid(
    list(paste0("f1", d$f1), paste0(":f2", d$f2)),
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )
  all_combis <- unique(paste0(combinations$Var1, combinations$Var2))
  expect_true(
    nrow(interactions[!interactions$interaction %in% all_combis, ]) == 0
  )
})
