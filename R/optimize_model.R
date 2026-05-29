#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables
#'  that do not significantly predict response variable.
#' @param formula
#'  The formula to be used with the model. Can be either quote() or formula().
#' @param data
#'  Dataframe with response and predictors as columns.
#' @param model_type
#'  Model type to be used.
#'  Options: "lm", "glm", "lmer", "glmer",
#'  "nlme", and "gam". Default: "glm".
#' @param family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options. Default: "gaussian".
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @return
#'  List with a) information on autocorrelated variables and b)
#'  final simplified/expanded models with further information
#'  (see function [LazyModeler::simplify_model()] for further details)
#' @export
optimize_model <- function(
    formula,
    data,
    model_type,
    family = gaussian,
    model_args = list(),
    direction = "backward") {
  
  # CHECKS
  check_model_family(family)
  check_model_type(model_type)
  
  # TODO AUTOCORRELATIONS
  
  # TODO MODEL SIMPLIFICATION
  simplify_model(formula,
                 data,
                 model_type,
                 family,
                 model_args,
                 direction)
  
  # TODO PLOTTING
}
