test_that("expanded_model_summary correctly identifies numeric/categorical vars and interactions", {
  d <- make_tiny_categorical_data()
  term = quote(y ~ f1 + f2 + f1:f2 + x1*x2 + x1:f1 + I(x1^2))
  m <- stats::glm(term, data = d)
  m_summary <- coef(summary(m))
  response_frm <- term[[2]]
  round_p <- 3

  m_overview <- expand_model_summary(m_summary,
                                     term,
                                     c("f1", "f2"),
                                     d,
                                     round_p)

  interactions <- dplyr::filter(m_overview, grepl(":", predictor))
  claimed_interactions <- dplyr::filter(m_overview[m_overview$is_interaction,], !grepl(":", predictor))
  non_interactions <- dplyr::filter(m_overview, !grepl(":", predictor))

  expect_all_true(interactions$is_interaction)
  expect_all_false(non_interactions$is_interaction)
  expect_true(nrow(claimed_interactions) == 0)
  expect_true(m_overview[(m_overview$main_effect1 == "x1") &
                           (!is.na(m_overview$main_effect1)),
                         "pred_type"][[1]] == "numeric")
  expect_true(m_overview[(m_overview$main_effect1 == "f2") &
                           (!is.na(m_overview$main_effect1)),
                         "pred_type"][[1]] == "categorical")
  expect_true(m_overview[((m_overview$main_effect2 == "x1") &
                            (!is.na(m_overview$main_effect2))) &
                         (m_overview$main_effect1 == "f1"),
                         "pred_type"][[1]] == "mixed")
})

test_that("expanded_model_summary correctly names main effects", {
  d <- make_tiny_categorical_data()
  term = quote(y ~ f1 + f2 + f1:f2 + x1*x2 + x1:f1 + I(x1^2))
  m <- stats::glm(term, data = d)
  m_summary <- coef(summary(m))
  response_frm <- term[[2]]
  round_p <- 3

  m_overview <- expand_model_summary(m_summary,
                              term,
                              c("f1", "f2"),
                              d,
                              round_p)

  expect_true(m_overview[m_overview$predictor == "f1B",
                         "main_effect1"][[1]] == "f1")
})