test_that("optimize_model does not throw error on
          gaussian backward simplification with psi", {
  data("plants")

  # during psi, some of the categorical traits are not significant
  # this causes NAs that are then replaced by 1s -> throws warning
  # this is completely correct and the warning is expected
  expect_warning(
    res <- optimize_model(
      df = plants,
      term = quote(sexual_seed_prop ~ altitude +
              latitude_gps_n +
              longitude_gps_e +
              (solar_radiation +
                 annual_mean_temperature +
                 isothermality)^2 +
              I(isothermality^2) +
              habitat +
              ploidy),
      autocorrelation_cols = c(
        "solar_radiation",
        "annual_mean_temperature",
        "isothermality",
        "altitude",
        "latitude_gps_n",
        "longitude_gps_e"),
      automatic_removal = TRUE,
      autocorrelation_threshold = 0.8,
      correlation_method = "spearman",
      model_type = "glm",
      model_family = "gaussian",
      evaluation_methods = c('anova'),
      simplification_direction= "backward",
      omit_na='overall',
      scale_predictor=TRUE,
      plot_quality_assessment='performance',
      round_p = 3,
      cor_use='complete.obs',
      plot_relationships = TRUE,
      jitter_plots = TRUE,
      plot_type="violinplot",
      stat_test='t.test',
      trace=TRUE,
      use_psi = TRUE)
  )

  expect_equal(names(res),
               c("autocorrelations", "models_with_info"))
  expect_equal(names(res$models_with_info),
               c("backward"))
  expect_equal(names(res$models_with_info$backward),
               c("overview",
                 "final_model",
                 "plots",
                 "history",
                 "psi",
                 "model_before_psi"))
})

test_that("optimize_model does not throw error on gaussian forward simplification", {
  data("plants")

  res <- optimize_model(
    df = plants,
    term = quote(sexual_seed_prop ~ altitude +
                   latitude_gps_n +
                   longitude_gps_e +
                   (solar_radiation +
                      annual_mean_temperature +
                      isothermality)^2 +
                   I(isothermality^2) +
                   habitat +
                   ploidy),
    autocorrelation_cols = c(
      "solar_radiation",
      "annual_mean_temperature",
      "isothermality",
      "altitude",
      "latitude_gps_n",
      "longitude_gps_e"),
    automatic_removal = TRUE,
    autocorrelation_threshold = 0.8,
    correlation_method = "spearman",
    model_type = "glm",
    model_family = "gaussian",
    evaluation_methods = c('anova'),
    simplification_direction= "forward",
    omit_na='overall',
    scale_predictor=TRUE,
    plot_quality_assessment='performance',
    round_p = 3,
    cor_use='complete.obs',
    plot_relationships = TRUE,
    jitter_plots = TRUE,
    plot_type="violinplot",
    stat_test='t.test',
    trace=TRUE)

  expect_equal(names(res),
               c("autocorrelations", "models_with_info"))
  expect_equal(names(res$models_with_info),
               c("forward"))
  expect_equal(names(res$models_with_info$forward),
               c("overview",
                 "final_model",
                 "plots",
                 "history"))
})

test_that("optimize_model does not throw error on binomial backward simplification", {
  data("plants")

  res <- optimize_model(
    df = plants,
    term = quote(sexual_seed_prop ~ altitude +
                   latitude_gps_n +
                   longitude_gps_e +
                   (solar_radiation +
                      annual_mean_temperature +
                      isothermality)^2 +
                   I(isothermality^2) +
                   habitat +
                   ploidy),
    autocorrelation_cols = c(
      "solar_radiation",
      "annual_mean_temperature",
      "isothermality",
      "altitude",
      "latitude_gps_n",
      "longitude_gps_e"),
    automatic_removal = TRUE,
    autocorrelation_threshold = 0.8,
    correlation_method = "spearman",
    model_type = "glm",
    model_family = "quasibinomial",
    evaluation_methods = c('anova'),
    simplification_direction= "backward",
    omit_na='overall',
    scale_predictor=TRUE,
    plot_quality_assessment='performance',
    round_p = 3,
    cor_use='complete.obs',
    plot_relationships = TRUE,
    jitter_plots = TRUE,
    plot_type="violinplot",
    stat_test='wilcox',
    trace=TRUE)

  expect_equal(names(res),
               c("autocorrelations", "models_with_info"))
  expect_equal(names(res$models_with_info),
               c("backward"))
  expect_equal(names(res$models_with_info$backward),
               c("overview",
                 "final_model",
                 "plots",
                 "history"))
})

