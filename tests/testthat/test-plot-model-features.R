test_that("plot_model_features returns a named list of plots", {
  d <- make_tiny_data()
  term = quote(y ~ x1 + x3 + f1)
  m <- stats::lm(term, data = d)
  m_summary <- coef(summary(m))
  response_frm <- term[[2]]
  categorical_check = setup_categorical_check(d, c('f1'), quote(predictor))
  round_p <- 3
  m_overview <- slightly_expanded_model_summary(m_summary, response_frm, term, categorical_check, round_p, d)
  
  res <- plot_model_features(
    regression_model = m,
    models_overview = m_overview,
    model_type = 'lm',
    model_family = 'gaussian',
    jitter_plots = TRUE,
    plot_type = "violinplot",
    round_p = round_p
  )
  
  expect_type(res, "list")
  expect_true(length(names(res)) > 0)
  
  # plot object existing?
  expect_true(any(vapply(res, function(x) inherits(x, c("ggplot","recordedplot","grob","gtable")) || is.function(x), logical(1))))
})

test_that("plot_model_features does not error on simple model", {
  d <- make_tiny_categorical_data()
  term = quote(y ~ x1 + x3 + f1 + f2 + f1:x1 + x1:x3)
  m <- stats::glm(term, data = d)
  m_summary <- coef(summary(m))
  response_frm <- term[[2]]
  categorical_check = setup_categorical_check(d, c('f1', 'f2'), quote(predictor))
  round_p <- 3
  m_overview <- slightly_expanded_model_summary(m_summary, response_frm, term, categorical_check, round_p, d)

  expect_no_error(
    plot_model_features(regression_model = m,
                        models_overview = m_overview,
                        model_type = 'glm',
                        test = 't.test',
                        model_family = 'gaussian',
                        jitter_plots = TRUE,
                        plot_type = "boxplot",
                        round_p = round_p,
                        remove_insignificant = TRUE)
  )
})
