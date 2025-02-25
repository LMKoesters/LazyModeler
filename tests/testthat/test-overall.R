data(plants)
library(mgcv)

test_that("(g)lm gaussian backward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n	+ longitude_gps_e	+ (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "gaussian", 
                                 evaluation_methods = c('anova'), simplification_direction = "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='t.test', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
})

test_that("(g)lm gaussian forward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "gaussian", 
                                 evaluation_methods=c('anova'), simplification_direction= "forward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='t.test', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm quasibinomial backward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e	+ (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "quasibinomial", 
                                 evaluation_methods=c('anova'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm binomial backward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e	+ (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "binomial", 
                                 evaluation_methods=c('anova'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm binomial backward aic", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n	+ longitude_gps_e	+ (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "binomial", 
                                 evaluation_methods=c("aic"), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
})

test_that("glm binomial backward aicc", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "binomial", 
                                 evaluation_methods=c("aicc"), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm binomial backward bic", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "binomial", 
                                 evaluation_methods=c('bic'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm quasibinomial forward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glm", model_family = "quasibinomial", 
                                 evaluation_methods=c('anova'), simplification_direction= "forward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm quasibinomial backward (autocorrelation_threshold = 0.5, correlation_method = 'pearson', plot_quality_assessment='baseR', plot_type='boxplot')", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.5, correlation_method = "pearson", model_type = "glm", model_family = "quasibinomial", 
                                 evaluation_methods=c('anova'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='baseR', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="boxplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glm quasibinomial forward (autocorrelation_threshold = 0.5, correlation_method = 'pearson', plot_quality_assessment='baseR', plot_type='boxplot'", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	ploidy), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.5, correlation_method = "pearson", model_type = "glm", model_family = "quasibinomial", 
                                 evaluation_methods=c('anova'), simplification_direction= "forward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='baseR', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="boxplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("glmer quasibinomial backward", {
  test_model <- expect_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	(1|ploidy)),
                              autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"),
                              automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "glmer", model_family = "quasibinomial", 
                              evaluation_methods=c('anova'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                              round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='wilcox', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("lmer gaussian backward", {
  test_model <- expect_no_error(optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation	+ annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat	+	(1|ploidy)), 
                                 autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
                                 automatic_removal = TRUE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "lmer", model_family = "gaussian", 
                                 evaluation_methods = c('anova', 'aicc'), simplification_direction = "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                                 round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='t.test', backward_simplify_model = TRUE, trace=TRUE)
  )
  
  plots <- test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["regression_plots"]]
  for (p in names(plots)) {
    expect_no_error(print(plots[[p]]))
  }
  
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["estimate_plot"]]))
  expect_no_error(print(test_model[["models_with_info"]][["backward"]][["plots"]][["regression_plots"]][["effect_size_plot"]]))
  
})

test_that("gam backward", {
  gam_dat <- gamSim(1,n=400,dist="normal",scale=2)
  expect_no_error(
    optimize_model(gam_dat, quote(y ~ s(x0) + s(x1) + s(x2) + s(x3)),
                   autocorrelation_cols = NA, automatic_removal = FALSE, autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "gam", model_family = "gaussian", 
                   evaluation_methods=c('anova'), simplification_direction= "backward", omit.na='overall', scale_predictor=TRUE, plot_quality_assessment='performance', 
                   round_p = 3, cor_use='complete.obs', plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='t.test', backward_simplify_model = TRUE, trace=TRUE)
  )
})

test_that("nlmer backward", {
  data("Orange")
  startvec <- c(Asym = 200, xmid = 725, scal = 350)
  expect_no_error(optimize_model(Orange, quote(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree), 
                                 start = startvec, non_linear = quote(circumference ~ SSlogis(age, Asym, xmid, scal)), random = quote(Asym ~ 1 |Tree), fixed = quote(Asym + xmid + scal ~ 1),
                                 autocorrelation_cols = NA, automatic_removal = FALSE, 
                                 autocorrelation_threshold = 0.8, correlation_method = "spearman", model_type = "nlmer", 
                                 evaluation_methods=c('anova'), simplification_direction = "backward", omit.na='overall', 
                                 scale_predictor=TRUE, plot_quality_assessment='performance', round_p = 3, cor_use='complete.obs', 
                                 plot_relationships = TRUE, jitter_plots = TRUE, plot_type="violinplot", stat_test='t.test', backward_simplify_model = TRUE, trace=TRUE)
  )
})
