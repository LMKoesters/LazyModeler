#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables
#'  that do not significantly predict response variable.
#' @param formula
#'  The formula to be used with the model. Can be either quote() or formula().
#' @param data
#'  Dataframe with response and predictors as columns.
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "lmer", "glmer",
#'  "nlme", "gam", and "nls"
#' @param family
#'  A character string or call describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param directions
#'  Vector of directions of model selection. Default: c("backward")
#' @param simplify_model
#'  Whether or not to simplify the full model. Default: TRUE
#' @param scale_predictors
#'  Whether to apply scaling to predictor variables. Default: FALSE
#' @param remove_autocors
#'  Whether to remove autocorrelated dataframe variables from formula.
#'    Default: TRUE
#' @param trace
#'  Whether to return the full model selection history. Default: TRUE
#' @param base_formula
#'  Lower formula used for forward model selection. Only required if forward
#'    selection is chosen, otherwise NA
#' @param use_psi
#'  Whether to apply post model selection inference to correct p-values.
#'    Default: FALSE
#' @param quality_assessment
#'  The mode of model quality assessment. Either "baseR" or "performance".
#'    Default: "baseR"
#' @param ac_threshold
#'  Threshold at which two variables are to be considered autocorrelated.
#'    Default: 0.7
#' @param ac_columns
#'  Columns to check for autocorrelations. The order of columns dictates
#'  priority basis for removal of predictors. Columns further down the list
#'  are removed first. If empty, columns are determined based on the given
#'  dataframe
#' @param cor_args
#'  Further arguments for [stats::cor()].
#'    Default: method = "pearson" and use = "complete.obs"
#' @param p_threshold
#'  P-value threshold for significance evaluation
#' @param psi_boot_repl
#'  A number or list of psi bootstrap replicates.
#'  Default: 100
#' @param psi_k
#'  The multiple of the number of degrees of freedom used as
#'  penalty in the model selection. The default k = 2 corresponds to the AIC.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param stat_type
#'  Type of Anova test
#' @param psi_label_size
#'  Size of labels within post-selection inference plot.
#'  Default: 2.5
#' @param plot_point_position
#'  Position adjustment for model feature plots. See the position paramter of
#'    [ggplot2::geom_point()] for more information. Default: "jitter"
#' @param plot_relationships
#'  Whether to plot regression, effect size, and estimates.
#'  Default: TRUE
#' @param categorical_stat_test
#'  Either "t.test" or "wilcox".
#'  Used to calculate statistics for regression plots of categorical variables.
#'  Default: "wilcox"
#' @param plot_type
#'  Either "boxplot" or "violin".
#'  Used to plot regression plots for categorical variables.
#'  Default: "boxplot"
#' @param plot_curve
#'  Whether to plot [ggplot2::geom_smooth()] in regression plots. Default: TRUE
#' @return
#'  List with a) information on autocorrelated variables and b)
#'  final simplified/expanded models with further information and plots
#' @examples
#' # setup
#' data("plants")
#'
#' # generate a glm model with the provided term and simplify it
#' # by applying backward simplification
#' simplified_model_info <- optimize_model(
#'   sexual_seed_prop ~ solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     I(isothermality^2) +
#'     habitat +
#'     ploidy +
#'     solar_radiation:annual_mean_temperature +
#'     solar_radiation:isothermality +
#'     annual_mean_temperature:isothermality,
#'   data = plants,
#'   model_type = "glm",
#'   ac_threshold = 0.8,
#'   ac_columns = c("solar_radiation",
#'                  "annual_mean_temperature",
#'                  "isothermality",
#'                  "altitude",
#'                  "latitude_gps_n",
#'                  "longitude_gps_e"),
#'   cor_args = list(method = c("spearman"),
#'                   use = "complete.obs"),
#'   family = "quasibinomial",
#'   directions = c("backward"),
#'   scale_predictors = TRUE,
#'   evaluation_methods = c("anova"),
#'   quality_assessment = "performance",
#'   categorical_stat_test = "wilcox",
#'   plot_type = "violin",
#'   round_p = 3,
#'   trace = TRUE
#' )
#' @export
optimize_model <- function(
    formula,
    data,
    model_type,
    family = stats::gaussian,
    model_args = list(),
    evaluation_methods = c("anova"),
    directions = c("backward"),
    simplify_model = TRUE,
    scale_predictors = TRUE,
    remove_autocors = TRUE,
    trace = TRUE,
    base_formula = NA,
    use_psi = FALSE,
    quality_assessment = "baseR",
    ac_threshold = 0.7,
    ac_columns = c(),
    cor_args = list(method = c("pearson"),
                    use = "complete.obs"),
    p_threshold = 0.05,
    psi_boot_repl = 100,
    psi_k = 2,
    round_p = 5,
    stat_type = 2,
    psi_label_size = 2.5,
    plot_point_position = "jitter",
    plot_relationships = TRUE,
    categorical_stat_test = "wilcox",
    plot_type = "boxplot",
    plot_curve = TRUE) {

  check_model_type(model_type)
  if (!model_type %in% c("glm", "glmer", "gam")) {
    family <- check_model_family(family)
  }
  formula <- check_formula(formula, data)

  out <- list()

  # AUTOCORRELATIONS
  if ((length(ac_columns) == 0) || scale_predictors) {
    ac_columns <- formula_related_cols(formula, data)
  }
  autocorrelations_result <- handle_autocorrelations(formula,
                                                     data,
                                                     ac_columns,
                                                     remove = remove_autocors,
                                                     threshold = ac_threshold,
                                                     cor_args = cor_args)
  out$autocorrelation_result <- autocorrelations_result
  formula <- autocorrelations_result$formula

  # SCALING
  if (scale_predictors) {
    original_data <- data
    for (numerical_var in ac_columns) {
      data[numerical_var] <- as.vector(scale(data[, numerical_var]))
    }
  }

  # MODEL SIMPLIFICATION
  model_out <- list()
  for (direction in directions) {
    model_out[[direction]] <- list()
    if (simplify_model) {
      res <- simplify_model(formula,
                            data,
                            model_type,
                            model_args,
                            evaluation_methods,
                            directions,
                            family,
                            trace,
                            base_formula)
      model_out[[direction]]$model_selection_result <- res
      
      # PSI
      if ((model_type %in% c("glm", "lm")) && use_psi) {
        final_formula <- res$final_model$formula
        final_data <- stats::model.frame(res$final_model)
        final_p_values <- res$p_values
        
        psi_result <- run_psi(final_formula,
                              final_data,
                              final_p_values,
                              model_type,
                              family,
                              model_args,
                              psi_boot_repl,
                              psi_k,
                              round_p,
                              stat_type,
                              p_threshold,
                              psi_label_size)
        
        model_out[[direction]]$psi_result <- psi_result
        model_to_plot <- psi_result$psi_model
      } else {
        model_to_plot <- res$final_model
      }
    } else {
      model_to_plot <- create_model(
        formula,
        data,
        model_type,
        family,
        model_args
      )
      model_out[[direction]]$final_model <- model_to_plot
    }

    if (model_type %in% c("glm", "lm", "gam")) {
      if (scale_predictors) data <- original_data
      plots <- plot_model(model_to_plot,
                          model_type,
                          quality_assessment,
                          plot_relationships,
                          categorical_stat_test,
                          plot_type,
                          plot_curve,
                          round_p,
                          plot_point_position)
      model_out[[direction]]$plots <- plots
    }
  }
  
  out$models_with_info <- model_out
  out
}
