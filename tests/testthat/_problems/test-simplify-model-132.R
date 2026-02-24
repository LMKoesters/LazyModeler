# Extracted from test-simplify-model.R:132

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "LazyModeler", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- make_tiny_glm_data()
res <- simplify_model(
    df = d,
    model_type = "glm",
    term = y ~ x1 + x2 + f1,
    evaluation_methods = c("aic"),
    direction = "backward",
    categorical_vars = "f1",
    backward_simplify_model = TRUE,
    trace = FALSE
  )$backward
m <- res$final_model %||% res$model %||% res$fit
expect_s3_class(m, "glm")
expect_identical(stats::family(m)$family, "binomial")
