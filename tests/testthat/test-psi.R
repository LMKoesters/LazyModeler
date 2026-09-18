test_that("Using PSI returns PSI model info", {
  d <- make_significant_factors_data()
  formula <- y ~ f2 + f1 + x1 + x2 + x1:f1
  model_type <- "glm"

  final_model <- create_model(formula, d, model_type = "glm")
  p_values <- get_model_p_values(final_model,
                                 model_type)

  res <- run_psi(
    final_model,
    d,
    final_p_values = p_values,
    model_type = model_type
  )

  expect_type(res, "list")
  expect_named(res, c("stat_result", "psi_model", "psi_plot"))
})
