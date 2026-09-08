#' Create model
#'
#' Creates a regression model
#' @param formula
#'  The underlying formula for model creation as a language object
#' @param data
#'  The underlying data for model creation
#' @param model_type
#'  The type of model to create. Can be either "glm", "lm", "glmer", "lmer",
#'    "gam", "nlme", or "nls"
#' @param family
#'  The model family (distribution)
#' @param model_args
#'  Additional arguments that should be given to the model call directly
#' @param fit
#'  Whether to actually fit the model. Default: TRUE
#' @returns A regression model
#' @examples
#'  # setup
#'  data("plants")
#'  # generate an example glm model for plotting
#'  final_model <- create_model(
#'   sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     habitat +
#'     ploidy +
#'     solar_radiation:isothermality,
#'   data = plants,
#'   model_type = "glm",
#'   family = "quasibinomial"
#'  )
#' @export
create_model <- function(formula,
                         data,
                         model_type,
                         family = stats::gaussian,
                         model_args = list(),
                         fit = TRUE) {
  if (!"na.action" %in% names(model_args)) {
    model_args$na.action <- stats::na.omit
  }

  if (model_type %in% c("glm", "glmer", "gam")) {
    model_args$formula <- formula
    model_args$family <- family
    model_args$data <- data
  } else if (model_type == "nlme") {
    model_args$model <- formula
    model_args$data <- data
  } else {
    model_args$formula <- formula
    model_args$data <- data
  }

  model_fun <- switch(
    as.character(model_type),
    "glm" = {
      if (!fit) model_args$method <- "model.frame"
      stats::glm
    },
    "lm" = {
      if (!fit) model_args$method <- "model.frame"
      stats::lm
    },
    "glmer" = {
      lme4::glmer
    },
    "lmer" = {
      model_args$REML <- FALSE
      lme4::lmer
    },
    "gam" = {
      if (!fit) model_args$fit <- fit
      model_args$method <- "REML"
      mgcv::gam
    },
    "nlme" = {
      nlme::nlme
    },
    "nls" = {
      model_args$model <- TRUE
      stats::nls
    }
  )

  do.call(model_fun, model_args)
}
