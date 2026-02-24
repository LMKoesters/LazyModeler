# Extracted from test-psi.R:24

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "LazyModeler", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- make_tiny_data()
res <- optimize_model(
    df = d,
    term = quote(y ~ x1 + x2 + x3 + f1 + f1:x1),
    autocorrelation_cols = c("x1", "x2", "x3"),
    automatic_removal = TRUE,
    autocorrelation_threshold = 0.8,
    correlation_method = "spearman",
    model_type = "lm",
    evaluation_methods = c("aic"),
    simplification_direction = "backward",
    omit_na = "overall",
    scale_predictor = FALSE,
    plot_quality_assessment = "baseR",
    plot_relationships = FALSE,
    use_psi = TRUE,
    trace = FALSE
  )
expect_type(res, "list")
expect_named(res$models_with_info$backward,
               c("overview", "final_model", "plots", "psi"))
