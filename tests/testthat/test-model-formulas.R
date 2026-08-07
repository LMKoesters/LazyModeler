test_that("glm with cbind on lhs computes", {
  d <- make_grouped_data()

  expect_no_error({
    res <- create_model(
      formula = cbind(success, failure) ~ x,
      data = d,
      model_type = "glm",
      family = binomial,
    )
  })

  expect_named(stats::coef(res), c("(Intercept)",
                                   "x"))
  res_frame <- model.frame(res)
  expect_equal(colnames(model.response(res_frame)),
               c("success", "failure"))
})
