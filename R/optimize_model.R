#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables
#'  that do not significantly predict the response variable.
#' @param formula
#'  The formula to be used with the model. Can be either quote() or formula().
#' @param data
#'  Dataframe with response and predictors as columns. Note that we require
#'    regular sequential rownames. Custom rownames will be lost.
#' @param model_type
#'  Model type to be used as character string.
#'    Options: "lm", "glm", "lmer", "glmer",
#'    "nlme", "gam", and "nls".
#' @param ...
#'  Must be empty. Unsupported arguments will result in an error.
#' @param family
#'  A character string or call describing the family used for model calculation.
#'    See [stats::family] for options. Can also be "automatic".
#'    Default: gaussian
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'    Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'    Default: c("anova")
#' @param directions
#'  Vector of directions of model selection. Default: c("backward")
#' @param simplify_model
#'  Whether or not to simplify the full model. Default: TRUE
#' @param scale_predictors
#'  Whether to apply scaling to predictor variables. Default: TRUE
#' @param detect_autocors
#'  Whether to detect autocorrelated dataframe variables. Note that we
#'    support autocorrelation detection only for model types "lm", "glm",
#'    "lmer", "glmer", and "gam".
#'    Default: TRUE
#' @param remove_autocors
#'  Whether to remove autocorrelated dataframe variables from formula.
#'    Default: TRUE
#' @param trace
#'  Whether to return the full model selection history. Default: TRUE
#' @param base_formula
#'  Lower formula used for forward model selection. Only required if forward
#'    selection is chosen, otherwise NA.
#' @param use_psi
#'  Whether to apply post model selection inference to correct p-values.
#'    Default: FALSE
#' @param quality_assessment
#'  The mode of model quality assessment. Either "baseR" or "performance".
#'    Other values do not produce quality assessment plots.
#'    Default: "baseR"
#' @param ac_threshold
#'  Threshold at which two variables are to be considered autocorrelated.
#'    Default: 0.7
#' @param ac_columns
#'  Columns to check for autocorrelations. The order of columns dictates
#'    priority basis for removal of predictors. Columns further down the list
#'    are removed first. If empty, columns are determined based on the given
#'    dataframe. Provided columns need to be part of the formula.
#'    For instance, when providing the formula y ~ x1 + x2 and ac_columns
#'    c("x1", "x2", "x3"), x3 will not be considered for autocorrelation
#'    testing.
#'    Note that LazyModeler still distinguishes between interactions
#'    and main effects. Removal of predictors will start with the evaluation of
#'    correlations between main effects and will remove blocking interactions
#'    when removing a main effect.
#' @param cor_args
#'  Further arguments for [stats::cor()] and [stats::cor.test()].
#'    Default: method = "pearson" and use = "complete.obs"
#' @param p_threshold
#'  p-value threshold for significance evaluation during autocorrelation
#'    testing and model simplification.
#'    Default: 0.05
#' @param delta
#'  Used as a minimal distance between model performances that needs to be
#'    present for a candidate model to be considered an improvement over the
#'    last computed model within the model selection process. The logic goes
#'    like this: During both backward and forward selection, the bigger
#'    (more complex) model needs to show substantially better metrics than the
#'    smaller (less complex) model, otherwise the smaller model is selected.
#'    A bigger delta will require an even more substantial improvement over the
#'    smaller model for the bigger model to be chosen.
#'    Default: 2
#' @param psi_boot_repl
#'  A number or list of psi bootstrap replicates.
#'    Default: 100
#' @param psi_k
#'  The multiple of the number of degrees of freedom used as
#'    penalty in the model selection. The default k = 2 corresponds to the AIC.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param stat_type
#'  Type of Anova test
#' @param psi_label_size
#'  Size of labels within post-selection inference plot.
#'    Default: 2.5
#' @param plot_point_position
#'  Position adjustment for model feature plots. See the position paramter of
#'    [ggplot2::geom_point()] for more information. Default: "jitter"
#' @param plot_relationships
#'  Whether to plot regression, effect size, and estimates.
#'  Default: TRUE
#' @param categorical_stat_test
#'  Either "t.test" or "wilcox". Used to calculate statistics for regression
#'    plots of categorical variables.
#'    Default: "wilcox"
#' @param plot_type
#'  Either "boxplot" or "violin".
#'    Used to plot regression plots for categorical variables.
#'    Default: "boxplot"
#' @param plot_curve
#'  Whether to plot [ggplot2::geom_smooth()] in regression plots. Default: TRUE
#' @return
#'  List with a) information on autocorrelated variables and b)
#'    final simplified/expanded models with further information and plots
#' @examples
#' # Example 1: Backward glm model simplification using example plant data
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
#'
#' # Example 2: glmer model optimization
#' set.seed(42)
#'
#' n_groups <- 20
#' n_per_group <- 20
#' n <- n_groups * n_per_group
#'
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' grp <- factor(rep(seq_len(n_groups),
#'                   each = n_per_group))
#' group_effect <- rnorm(n_groups,
#'                       mean = 0,
#'                       sd = 0.6)
#' eta <- 0.5 +
#'   0.4 * x1 -
#'   0.3 * x2 +
#'   0.2 * x3 +
#'   group_effect[grp]
#' y <- rpois(n,
#'            lambda = exp(eta))
#' p <- plogis(eta)
#' y_binom <- rbinom(n, size = 1, prob = p)
#'
#' d <- data.frame(
#'   y = y,
#'   y_binom = y_binom,
#'   x1 = x1,
#'   x2 = x2,
#'   x3 = x3,
#'   grp = grp
#' )
#'
#' res <- optimize_model(
#'   formula = y ~ x1 + x2 + x3 + x1:x3 + (1 | grp),
#'   data = d,
#'   model_type = "glmer",
#'   model_args = list(),
#'   evaluation_methods = c("anova"),
#'   directions = c("backward"),
#'   family = poisson
#' )
#'
#' # Example 3: nls model optimization
#' set.seed(42)
#' x <- seq(0, 10, length.out = n)
#' asym <- 5
#' k <- 0.6
#' y <- asym * (1 - exp(-k * x)) + rnorm(n, sd = 0.15)
#'
#' d <- data.frame(
#'   y = y,
#'   x = x
#' )
#' start <- c(Asym = 5, k = 0.6, offset = .1, slope = .5)
#'
#' res <- optimize_model(
#'   formula = y ~ offset + Asym * (1 - exp(-k * x)) + slope * x,
#'   data = d,
#'   model_type = "nls",
#'   model_args = list(
#'     start = start
#'   ),
#'   evaluation_methods = c("anova"),
#'   directions = c("backward"),
#'   scale_predictors = FALSE
#' )
#' @export
optimize_model <- function(
    formula,
    data,
    model_type,
    ...,
    family = stats::gaussian,
    model_args = list(),
    evaluation_methods = c("anova"),
    directions = c("backward"),
    simplify_model = TRUE,
    scale_predictors = TRUE,
    detect_autocors = FALSE,
    remove_autocors = FALSE,
    trace = TRUE,
    base_formula = NA,
    use_psi = FALSE,
    quality_assessment = "baseR",
    ac_threshold = 0.7,
    ac_columns = c(),
    cor_args = list(method = c("pearson"),
                    use = "complete.obs"),
    p_threshold = 0.05,
    delta = 2,
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
  data <- as.data.frame(data)
  if (ncol(data) < 2) {
    stop(
      paste(
        "Your dataframe contains less than 2 columns. You need at least 2",
        "columns (representing the response and the predictor(s))."
      )
    )
  }

  check_user_args(sys.call(), LazyModeler::optimize_model)

  check_model_type(model_type, model_args)
  formula <- check_formula(formula, data)
  if (model_type %in% c("glm", "glmer", "gam")) {
    family <- check_model_family(model_type,
                                 family,
                                 automatic = identical(family, "automatic"),
                                 data = data,
                                 lhs = formula.tools::lhs(formula))
  }
  out <- list()
  rownames(data) <- seq_len(nrow(data))

  # AUTOCORRELATIONS
  autocor_supported <- model_type %in% c("lm", "glm", "lmer", "glmer", "gam")
  if (detect_autocors && autocor_supported) {
    autocorrelations_result <- handle_autocorrelations(
      formula,
      data,
      model_type,
      family,
      ac_columns,
      remove = remove_autocors,
      threshold = ac_threshold,
      p_threshold = p_threshold,
      cor_args = cor_args,
      model_args = model_args
    )
    out$autocorrelation_result <- autocorrelations_result
    formula <- autocorrelations_result$formula
  } else if (detect_autocors && !autocor_supported) {
    warning(
      paste(
        "Unfortunately, we do not support automatic autocorrelation detection",
        "for your chosen model type. Consider checking for autocorrelations",
        "yourself."
      )
    )
  }

  # SCALING
  if (scale_predictors) {
    original_data <- data
    formula_numeric_columns <- formula_numeric_cols(formula, data)
    for (numerical_var in formula_numeric_columns) {
      data[numerical_var] <- as.vector(scale(data[, numerical_var]))
    }
  }

  # MODEL SIMPLIFICATION
  for (direction in directions) {
    model_out <- list()
    if (simplify_model) {
      res <- simplify_model(formula,
                            data,
                            model_type,
                            model_args,
                            evaluation_methods,
                            direction,
                            family,
                            p_threshold,
                            trace,
                            base_formula,
                            delta)
      model_out$model_selection_result <- res

      # PSI
      if ((model_type %in% c("glm", "lm")) && use_psi) {
        final_p_values <- res$p_values

        psi_result <- run_psi(res$final_model,
                              data,
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

        model_out$psi_result <- psi_result
        model_to_plot <- psi_result$psi_model
      } else if (use_psi) {
        warning(
          paste(
            "Unfortunately, we currently only allow for PSI for glm and lm",
            "models. We will therefore proceed without PSI."
          )
        )
        model_to_plot <- res$final_model
      } else {
        model_to_plot <- res$final_model
      }
      model_out$final_model <- model_to_plot
    } else {
      model_to_plot <- create_model(
        formula,
        data,
        model_type,
        family,
        model_args
      )
      model_out$final_model <- model_to_plot
    }

    plotting_allowed <- c("glm", "lm", "glmer", "lmer")
    if (model_type %in% plotting_allowed && plot_relationships) {
      plots <- plot_model(model_to_plot,
                          if (scale_predictors) original_data else data,
                          model_type,
                          family,
                          model_args,
                          quality_assessment,
                          categorical_stat_test,
                          plot_type,
                          plot_curve,
                          round_p,
                          plot_point_position,
                          p_threshold)
      model_out$plots <- plots
    } else if (plot_relationships) {
      warning(
        paste(
          "We currently do not allow for plotting of your model type.",
          "However, plotting for gam models will be added soon."
        )
      )
    }

    out[[direction]] <- model_out
  }

  out
}
