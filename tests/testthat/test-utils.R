test_that("3-way interactions are extracted", {
  terms <- attr(stats::terms.formula(y ~ f1:f2:x1), "term.labels")
  interactions <- extract_interactions(terms)

  expect_true(length(unique(interactions$main_effect)) == 3)
})

test_that("3-way interactions are transformed into trait interactions", {
  d <- make_significant_factors_data()

  m <- create_model(
    formula = y ~ x1 + f2 + (x1 * x3 * f2)^3,
    data = d,
    model_type = "glm",
    model_args = list()
  )

  c(plot_data,
    categorical_vars,
    numeric_vars,
    mixed_interactions,
    term_map) %<-% extract_terms_per_category(
    d,
    m,
    "glm",
    "gaussian",
    list()
  )

  expect_length(setdiff(term_map$column,
                        c("x1", "f2medium", "f2deep", "x3",
                          "x1:x3", "x1:f2medium", "x1:f2deep",
                          "f2medium:x3", "f2deep:x3", "x1:f2medium:x3",
                          "x1:f2deep:x3")),
                0)
})

test_that("NAs are omitted from data", {
  d <- make_tiny_data_with_na()

  res <- omit_na_from_model_data(
    formula = prop ~ x1 + x2 + x3,
    data = d,
    model_type = "glm",
    family = quasibinomial,
    model_args = list(weights = quote(trials))
  )

  expect_shape(res, nrow = 55)
  expect_equal(colnames(res),
               colnames(d))
})

test_that("NAs are omitted from data when subset is used", {
  d <- make_tiny_data_with_na()

  res <- omit_na_from_model_data(
    formula = prop ~ x1 + x2 + x3,
    data = d,
    model_type = "glm",
    family = quasibinomial,
    model_args = list(weights = quote(trials),
                      subset = quote(x3 <= 0.))
  )

  expect_shape(res, nrow = 35)
  expect_equal(colnames(res),
               colnames(d))
})
