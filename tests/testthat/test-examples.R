data(plants)

test_that("simplify model example", {
  expect_no_error(simplify_model(plants, "glm", quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + 
                                                        isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + 
                                                        solar_radiation:isothermality + annual_mean_temperature:isothermality), 
                                 c("anova"), model_family='quasibinomial', direction='backward',
                                 categorical_vars=c('habitat', 'ploidy'), backward_simplify_model=TRUE, trace=TRUE, omit.na='overall'))
})

test_that("categorical check", {
  expect_no_error(setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor)))
})

test_that("slightly expanded model summary", {
  model_family = 'quasibinomial'
  categorical_check = setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor))
  final_model = glm(sexual_seed_prop ~ altitude + solar_radiation + 
                      annual_mean_temperature + isothermality + habitat + ploidy + 
                      solar_radiation:isothermality, family=model_family, data=plants)
  models_overview = coef(summary(final_model))
  expect_no_error(slightly_expanded_model_summary(models_overview, 
                                                  quote(sexual_seed_prop), 
                                                  quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + 
                                                          isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + 
                                                          solar_radiation:isothermality + annual_mean_temperature:isothermality), 
                                                  categorical_check, 
                                                  3, 
                                                  plants))
})

test_that("plot model features", {
  model_type = 'glm'
  model_family = 'quasibinomial'
  categorical_check = setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor))
  final_model = glm(sexual_seed_prop ~ altitude + solar_radiation +
                      annual_mean_temperature + isothermality + habitat + ploidy +
                      solar_radiation:isothermality, family=model_family, data=plants)

  models_overview = coef(summary(final_model))
  models_overview = slightly_expanded_model_summary(models_overview, 
                                                  quote(sexual_seed_prop), 
                                                  quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + 
                                                          isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + 
                                                          solar_radiation:isothermality + annual_mean_temperature:isothermality), 
                                                  categorical_check, 
                                                  3,
                                                  plants)
  models_overview = models_overview |>
    mutate(effect_direction = case_when(
      Estimate < 0 ~ 'negative',
      TRUE ~ 'positive'
    ))

  expect_no_error(plot_model_features(final_model,
                                      models_overview,
                                      model_type,
                                      model_family,
                                      plot_type='boxplot',
                                      jitter_plots=TRUE,
                                      test='wilcox',
                                      round_p=3,
                                      remove_insignificant=FALSE))
})

