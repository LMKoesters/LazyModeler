test_that("remove_autocorrelations returns expected structure", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 * x2 + (x3 + x4 + x5)^3,
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("x1", "x2", "x3", "x4", "x5"),
    remove = TRUE,
    threshold = 0.8
  )

  # check output type
  expect_type(res, "list")
  # check output structure
  expect_true(all(c("removed_predictors",
                    "autocorrelations_info",
                    "formula") %in% names(res)))
  # check output type
  expect_true(is.character(res$removed_predictors))
  # check output type
  expect_true(is.data.frame(res$autocorrelations_info))
  # check autocorrelations columns
  expect_true(all(c("coefficientA",
                    "coefficientB",
                    "correlation",
                    "p_value",
                    "note") %in% names(res$autocorrelations_info)))
})

test_that("no autocorrelations -> removed_predictors is empty", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x5 + x6,
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("x1", "x5", "x6"),
    remove = TRUE,
    threshold = 0.8
  )

  expect_null(res$autocorrelations_info)
  expect_length(res$removed_predictors, 0)
})

test_that("handle_autocorrelations recognizes priority order of variables", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x2 + x1 + x3,
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("x2", "x1", "x3"),
    remove = TRUE,
    threshold = 0.8
  )

  # check priority: x2 before x1 means that x1 should be removed
  expect_true("x1" %in% res$removed_predictors)
  expect_false("x2" %in% res$removed_predictors)
})

test_that("handle_autocorrelations recognizes A-B-C constellation", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x2 + x3 + x4 + x5,
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("x1", "x2", "x3", "x4", "x5"),
    remove = TRUE,
    threshold = 0.8
  )

  # check that x2 was removed (a-b-c test)
  expect_true("x2" %in% res$removed_predictors)
})

test_that("handle_autocorrelations respects main effects", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x2 + x3 + x4 + x5 + x1:x2 + x1:x3 + f1 +
      f1:log(x1) + log(x1),
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("f1", "x1:x2", "x1:x3", "x1", "x4", "x2", "x5", "x3",
             "f1:log(x1)", "log(x1)"),
    remove = TRUE,
    threshold = 0.8
  )

  expect_equal(res$removed_predictors,
               c("log(x1)", "x2", "x1:x2", "x3", "x1:x3", "x5", "f1:log(x1)"))
  expect_equal(res$formula,
               y ~ x1 + x4 + f1)
})

test_that("handle_autocorrelations ignores extra column", {
  d <- make_autocor_data()

  res <- handle_autocorrelations(
    formula = y ~ x1 + x2 + x3 + I(x^2) + x1:x2 + log(x1),
    data = d,
    model_type = "glm",
    family = gaussian,
    cols = c("x1", "x4", "x2", "x3"),
    remove = TRUE,
    threshold = 0.8
  )

  expect_true((!"x4" %in% res$autocorrelations$coefficientA) &&
                (!"x4" %in% res$autocorrelations$coefficientB))
})

test_that("sort_term_map respects given order under (non)-interactions", {
  d <- make_autocor_data()
  formula <- y ~ x1:x2 + x1 + x2 + x3 + x5 + x4 + x1:x3 + f1

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "glm",
    family = "quasibinomial",
    model_args = list(
      subset = quote(x1 > 0),
      contrasts = list(
        f1 = "contr.sum"
      )
    )
  )

  term_map <- sort_term_map(term_map)
  expect_equal(unique(term_map$term),
               c("x1", "x2", "x3", "x5", "x4", "f1", "x1:x2", "x1:x3"))
})

test_that("invalid input of 1 column throws informative error", {
  d <- make_autocor_data()

  expect_error(
    handle_autocorrelations(
      formula = y ~ x1,
      data = d,
      model_type = "glm",
      family = gaussian,
      cols = c("does_not_exist")
    ),
    regexp = "at least two columns"
  )
})

test_that("invalid input of threshold > 1 throws informative error", {
  d <- make_autocor_data()

  expect_error(
    handle_autocorrelations(formula = y ~ x1 + x2 + x3, data = d,
                            cols = c("x1", "x2", "x3"), threshold = 1.5),
    regexp = "outside of range \\(0,1\\)"
  )
})

test_that("input of two invalid columns throw informative error", {
  d <- make_autocor_data()

  expect_error(
    handle_autocorrelations(
      formula = y ~ x1,
      data = d,
      model_type = "glm",
      family = gaussian,
      cols = c("does_not_exist", "also_does_not_exist")
    ),
    regexp = "unable to detect enough valid columns"
  )
})

test_that("error if autocorrelations but no automatic removal is requested", {
  d <- make_tiny_data()

  expect_warning(
    handle_autocorrelations(
      formula = y ~ x1 + x2 + x3,
      data = d,
      model_type = "glm",
      family = gaussian,
      cols = c("x1", "x2", "x3"),
      remove = FALSE,
      threshold = 0.8,
      cor_args = list(method = "spearman",
                      use = "complete.obs")
    ),
    regexp = "Some of your variables are autocorrelated"
  )
})

test_that("relevant autocorrelations columns for glm detected", {
  d <- make_tiny_data()
  formula <- y ~ x1 + x2 + f1 + I(x1^2) + x3:x2

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "glm",
    family = "quasibinomial",
    model_args = list(
      weights = quote(trials),
      subset = quote(x1 > 0),
      contrasts = list(
        f1 = "contr.sum"
      )
    )
  )

  expect_equal(colnames(m_matrix),
               c("x1", "x2", "f11", "f12", "I(x1^2)", "x2:x3"))
  expect_equal(term_map$term,
               c("x1", "x2", "f1", "f1", "I(x1^2)", "x2:x3"))
})

test_that("relevant autocorrelations columns for lm detected", {
  d <- make_significant_factors_data()
  formula <- y ~ x1 + x2 + f1 + I(x1^2) + x3:x2

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "lm",
    family = gaussian
  )

  expect_equal(colnames(m_matrix),
               c("x1", "x2", "f1B", "f1C", "I(x1^2)", "x2:x3"))
})

test_that("relevant autocorrelations columns for glmer detected", {
  d <- d <- make_lmer_data()
  formula <- y ~ x1 + x2 + I(x1^2) + x3:x2 + (1 | grp)

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "glmer",
    family = poisson,
    model_args = list(subset = quote(x1 > -0.5)
    )
  )

  expect_equal(colnames(m_matrix),
               c("x1", "x2", "I(x1^2)", "x2:x3"))
})

test_that("relevant autocorrelations columns for lmer detected", {
  d <- d <- make_lmer_data()
  formula <- y ~ x1 + x2 + I(x1^2) + x3:x2 + (1 | grp)

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "lmer",
    family = gaussian
  )

  expect_equal(colnames(m_matrix),
               c("x1", "x2", "I(x1^2)", "x2:x3"))
})

test_that("relevant autocorrelations columns for gam detected", {
  d <- make_gam_data()
  formula <- y ~ s(x1) + ti(x2, x3) + te(x1, x2) + t2(x2, x3) + x2 * x3

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "gam",
    family = gaussian
  )

  expect_equal(colnames(m_matrix),
               c("x2", "x3", "x2:x3"))
})


test_that("0 variance is rejected from autocorrelation data", {
  d <- make_tiny_data()
  d$x1 <- 1
  formula <- y ~ x1 + x2 + f1 + I(x1^2) + x3:x2

  c(m_matrix, term_map) %<-% format_cor_data(
    d,
    formula,
    model_type = "glm",
    family = "quasibinomial",
    model_args = list(
      weights = quote(trials),
      subset = quote(x1 > 0),
      contrasts = list(
        f1 = "contr.sum"
      )
    )
  )

  (c(term_map, has_no_variance) %<-% check_variance(m_matrix, term_map)) |>
    expect_warning(regexp = "no variance after NA omission")
  expect_equal(term_map$column, c("x2", "f11", "f12", "x2:x3"))
})

test_that("0 variance is rejected from autocorrelation data", {
  d <- make_tiny_data()
  d$x1 <- 1
  formula <- y ~ x1 + x2 + f1 + I(x1^2) + x3:x2

  (res <- handle_autocorrelations(
    formula,
    d,
    model_type = "glm",
    family = "quasibinomial",
    model_args = list(
      weights = quote(trials),
      subset = quote(x1 > 0),
      contrasts = list(
        f1 = "contr.sum"
      )
    )
  )) |>
    expect_warning(regex = "no variance after NA omission")

  expect_equal(res$removed_predictors,
               c("x1", "I(x1^2)"))
})

test_that("formula related columns are recognized from matrix", {
  d <- make_lmer_data()
  formula <- y ~ x1 + x2 + x1:x2 + (1 | grp)

  (res <- handle_autocorrelations(
    formula,
    d,
    model_type = "glmer",
    family = poisson,
    model_args = list(subset = quote(x1 > 0))
  )) |>
    expect_no_error()

  expect_equal(res$removed_predictors,
               c("x1:x2"))
})
