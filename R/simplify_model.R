#' Simplify model
#'
#' Wrapper function for backward simplification and forward selection of a model
#' @param formula
#'  Upper formula to be used for model creation/selection
#' @param data
#'  Underlying model data
#' @param model_type
#'  The type of model to be created. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param direction
#'  Direction of model selection. Default: "backward"
#' @param family
#'  The model family. Default: gaussian
#' @param trace
#'  Whether to return the selection history. Default. TRUE
#' @param base_formula
#'  The lower formula used for forward model selection. Only required if
#'    direction = "forward", otherwise this is NA
#' @returns
#'  The selected model with evaluation metrics and (if trace = TRUE)
#'    the selection history
#' @examples
#' # setup
#' data("plants")
#'
#' # generate a glm model with the provided term and simplify it
#' # by applying backward simplification
#' simplified_model_info <- simplify_model(
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
#'   evaluation_methods = c("anova"),
#'   family = "quasibinomial",
#'   direction = "backward",
#'   trace = TRUE
#' )
#' @export
simplify_model <- function(
    formula,
    data,
    model_type,
    model_args = list(),
    evaluation_methods = c("anova"),
    direction = "backward",
    family = stats::gaussian,
    trace = TRUE,
    base_formula = NA) {
  evaluation_methods <- tolower(evaluation_methods)
  formula <- stats::formula(stats::terms(stats::as.formula(formula),
                                         data = data))

  if (model_type %in% c("nls", "nlme")) {
    return(nls_nlme(formula,
                    data,
                    model_type,
                    model_args,
                    evaluation_methods))
  }

  if (model_type == "gam") {
    model <- create_model(formula, data, model_type, family, model_args,
                          fit = FALSE)

    data <- model$mf
  } else if (model_type %in% c("glm", "lm")) {
    data <- create_model(formula, data, model_type, family, model_args,
                         fit = FALSE)
  } else {
    data <- stats::model.frame(
      formula = reformulas::subbars(stats::as.formula(formula)),
      data = data,
      na.action = stats::na.omit
    )
  }

  if (direction == "backward") {
    out <- optimize_backward(formula, data, model_type, family, model_args,
                             evaluation_methods, trace)
  } else {
    if (typeof(base_formula) != "language") {
      warning(paste("You did not provide a base formula for forward model",
                    "selection. We are continuing with y ~ 1. If you wish",
                    "to use a different lower formula, please adjust your",
                    "method call.",
                    collapse = " "))
      base_formula <- y ~ 1
    }

    out <- optimize_forward(formula, base_formula, data, model_type, family,
                            model_args, evaluation_methods, trace)
  }
  out
}

#' nls/nlme model assessment
#'
#' nls and nlme models do not undergo LazyModeler's model selection procedure
#'  but are created and assessed before being returned
#' @param formula
#'  Formula to be used for model creation
#' @param data
#'  Underlying model data
#' @param model_type
#'  The type of model to be created. Can be either "nlme", or "nls"
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @returns A model with evaluation metrics
nls_nlme <- function(formula,
                     data,
                     model_type,
                     model_args = list(),
                     evaluation_methods = c("anova")) {

  model <- create_model(formula,
                        data,
                        model_type = model_type,
                        model_args = model_args)

  evaluation_methods <- evaluation_methods[evaluation_methods != "anova"]
  old_model_assessments <- as.list(rep(Inf, length(evaluation_methods)))
  old_model_assessments <- stats::setNames(old_model_assessments,
                                           evaluation_methods)
  assessed_models <- compare_models(evaluation_methods,
                                    list(model, model),
                                    "backward",
                                    old_model_assessments)
  p_values <- get_model_p_values(model, model_type)
  out <- list(assessments = assessed_models$assessments,
              final_model = model,
              p_values = p_values)
  out
}

#' Forward model selection
#'
#' Selects a model by sequentially adding terms to a lower formula
#' @param formula
#'  The upper formula containing all terms to be added
#' @param base_formula
#'  The lower formula used as a starting point
#' @param data
#'  Underlying data for model creation
#' @param model_type
#'  The type of model to create. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param family
#'  The model family. Default: gaussian
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param trace
#'  Whether to return the selection history. Default. TRUE
#' @returns
#'  The selected model with evaluation metrics and (if trace = TRUE)
#'    the selection history
optimize_forward <- function(
    formula,
    base_formula,
    data,
    model_type,
    family = stats::gaussian,
    model_args = list(),
    evaluation_methods = c("anova"),
    trace = TRUE) {
  history <- list()
  current_formula <- base_formula
  optimize <- "proceed"
  last_model_info <- NA

  model <- create_model(base_formula,
                        data,
                        model_type,
                        family,
                        model_args)

  c(predictors_to_add, p_col) %<-% get_addable_terms(formula, model, model_type)

  evaluation_wo_anova <- evaluation_methods[evaluation_methods != "anova"]
  old_model_assessments <- as.list(rep(Inf, length(evaluation_methods)))
  old_model_assessments <- stats::setNames(old_model_assessments,
                                           evaluation_methods)
  assessed_models <- compare_models(evaluation_wo_anova,
                                    list(model, model),
                                    "forward",
                                    old_model_assessments)

  last_model_info <- list(
    model = model,
    assessments = assessed_models$assessments
  )
  if (trace) history[[deparse1(current_formula)]] <- last_model_info

  while (optimize == "proceed") {
    c(best_pred_info,
      best_pred,
      best_step_res,
      optimize) %<-% determine_best_pred(predictors_to_add,
                                         current_formula,
                                         formula,
                                         data,
                                         model_type,
                                         last_model_info,
                                         direction = "forward",
                                         family,
                                         model_args,
                                         evaluation_methods)

    if (typeof(best_pred_info) != "list") {
      if (optimize == "revert") {
        final_model = last_model_info$model
      } else {
        final_model = best_step_res$model
      }

      p_values <- get_model_p_values(final_model,
                                     model_type)

      out <- list(
        "final_model" = final_model,
        "p_values" = p_values
      )
      if (trace) out$history <- history
      return(out)
    }

    predictors_to_add <- best_step_res$adjustable_terms

    d <- paste(". ~ . +", best_pred)
    current_formula <- stats::update(stats::as.formula(current_formula), d)

    if (trace) history[[deparse1(current_formula)]] <- list(
      model = best_step_res$model,
      assessments = best_step_res$assessments
    )

    last_model_info <- best_step_res
  }
}

#' Determine best predictor to add/remove
#'
#' Determines the best predictor to add/remove given a list of possible
#'  predictors and methods on which to evaluate the model
#' @param predictors_to_adjust
#'  Predictors that can be removed or added at the current selection step
#' @param current_formula
#'  The formula of the previous model selection step
#' @param full_formula
#'  The full formula containing all terms
#' @param data
#'  Underlying data for model creation
#' @param model_type
#'  The type of model to create. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param last_model_info
#'  Assessment information of the previous model selection step
#' @param direction
#'  Direction of model selection. Default: "backward"
#' @param family
#'  The model family. Default: gaussian
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @returns
#'  The best predictor to add/remove, alongside the updated model and
#'    information on whether to continue the selection process or revert to
#'    the previous selection step
determine_best_pred <- function(predictors_to_adjust,
                                current_formula,
                                full_formula,
                                data,
                                model_type,
                                last_model_info,
                                direction,
                                family = stats::gaussian,
                                model_args = list(),
                                evaluation_methods = c("anova")) {
  best_pred_info <- NA
  best_pred <- NA
  best_step_res <- NA
  simplify_info <- c()
  adjuster <- if (direction == "forward") "+" else "-"

  for (predictor in predictors_to_adjust$predictor) {
    d <- paste(". ~ .", adjuster, predictor)
    formula_to_test <- stats::update(stats::as.formula(current_formula), d)

    step_res <- optimizer_step(formula_to_test,
                               data,
                               model_type,
                               last_model_info,
                               family,
                               model_args,
                               evaluation_methods,
                               direction = direction,
                               full_formula = full_formula)

    if (step_res$optimize == "proceed") {
      if (typeof(best_pred_info) == "list") {
        new_is_better <- new_model_is_better(step_res$assessments,
                                             best_pred_info$assessments)
        if (new_is_better) {
          best_pred_info <- list(assessments = step_res$assessments)
          best_pred <- predictor
          best_step_res <- step_res
        }
      } else {
        best_pred_info <- list(assessments = step_res$assessments)
        best_pred <- predictor
        best_step_res <- step_res
      }
    }

    simplify_info <- append(simplify_info, step_res$optimize)
  }

  if (typeof(best_step_res) == "list") {
    list(best_pred_info, best_pred, best_step_res, "proceed")
  } else {
    list(best_pred_info, best_pred, best_step_res, "revert")
  }
}

#' Backward model selection
#'
#' Selects a model by sequentially removing terms from a formula
#' @param formula
#'  The upper formula containing all terms to be added
#' @param data
#'  Underlying data for model creation
#' @param model_type
#'  The type of model to create. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param family
#'  The model family. Default: gaussian
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param trace
#'  Whether to return the selection history. Default. TRUE
#' @returns
#'  The selected model with evaluation metrics and (if trace = TRUE)
#'    the selection history
optimize_backward <- function(
    formula,
    data,
    model_type,
    family = stats::gaussian,
    model_args = list(),
    evaluation_methods = c("anova"),
    trace = TRUE) {
  history <- list()
  last_model_info <- NA
  optimize <- "proceed"
  current_formula <- formula

  model <- create_model(formula,
                        data,
                        model_type,
                        family,
                        model_args)

  c(predictors_to_remove, p_col) %<-% get_removable_terms(model,
                                                          model_type)

  evaluation_wo_anova <- evaluation_methods[evaluation_methods != "anova"]
  old_model_assessments <- as.list(rep(Inf, length(evaluation_methods)))
  old_model_assessments <- stats::setNames(old_model_assessments,
                                           evaluation_methods)
  assessed_models <- compare_models(evaluation_wo_anova,
                                    list(model, model),
                                    "backward",
                                    old_model_assessments)

  last_model_info <- list(
    model = model,
    assessments = assessed_models$assessments
  )
  if (trace) {
    history[[deparse1(current_formula)]] <- last_model_info
  }

  while (optimize == "proceed") {
    c(best_pred_info,
      best_pred,
      best_step_res,
      optimize) %<-% determine_best_pred(predictors_to_remove,
                                         current_formula,
                                         formula,
                                         data,
                                         model_type,
                                         last_model_info,
                                         direction = "backward",
                                         family,
                                         model_args,
                                         evaluation_methods)

    if (typeof(best_pred_info) != "list") {
      if (optimize == "revert") {
        final_model = last_model_info$model
      } else {
        final_model = best_step_res$model
      }

      p_values <- get_model_p_values(final_model,
                                     model_type)

      out <- list(
        "final_model" = final_model,
        "p_values" = p_values
      )
      if (trace) out$history <- history
      return(out)
    }

    predictors_to_remove <- best_step_res$adjustable_terms

    d <- paste(". ~ . -", best_pred)
    current_formula <- stats::update(stats::as.formula(current_formula), d)

    if (trace) history[[deparse1(current_formula)]] <- list(
      model = best_step_res$model,
      assessments = best_step_res$assessments
    )

    last_model_info <- best_step_res
  }
}

#' Optimization step
#'
#' Takes a step in the model selection process; either forward or backward
#' @param formula
#'  Formula to be used for model creation/selection
#' @param data
#'  Underlying data for model creation
#' @param model_type
#'  The type of model to create. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param last_model_info
#'  Assessment information of the previous model selection step
#' @param family
#'  The model family. Default: gaussian
#' @param model_args
#'  Optional additional arguments to be given directly to the model call
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param direction
#'  Direction of model selection. Default: "backward"
#' @param full_formula
#'  Upper formula to be used with forward model selection.
#'    Should be NA if direction = "backward"
optimizer_step <- function(formula,
                           data,
                           model_type,
                           last_model_info,
                           family = stats::gaussian,
                           model_args,
                           evaluation_methods,
                           direction = "backward",
                           full_formula = NA) {
  has_last_model <- typeof(last_model_info) == "list"
  if (!has_last_model) {
    model <- create_model(formula,
                          data,
                          model_type,
                          family,
                          model_args)
  } else {
    model <- stats::update(last_model_info$model,
                           formula = formula,
                           na.action = stats::na.fail)
  }

  if (direction == "backward") {
    c(adjustable_terms, p_col) %<-% get_removable_terms(model,
                                                        model_type)
  } else {
    c(adjustable_terms, p_col) %<-% get_addable_terms(full_formula,
                                                      model,
                                                      model_type)
  }

  out <- list(
    adjustable_terms = adjustable_terms,
    model = model,
    optimize = "proceed"
  )

  if (!has_last_model) {
    evaluation_methods <- evaluation_methods[evaluation_methods != "anova"]
    last_model_info <- out
    last_model_info$assessments <- stats::setNames(
      as.list(rep(Inf, length(evaluation_methods))),
      evaluation_methods
    )
  }

  assessed_models <- compare_models(evaluation_methods,
                                    list(model, last_model_info$model),
                                    direction,
                                    last_model_info$assessments)
  out$assessments <- assessed_models$assessments

  if (has_last_model &&
        !assessed_models$has_improved) out$optimize <- "revert"

  out
}
