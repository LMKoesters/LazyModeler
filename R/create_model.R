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
  model_args$na.action <- stats::na.omit

  model_fun <- switch(
    as.character(model_type),
    "glm" = {
      glm_method <- if (fit) "glm.fit" else "model.frame"
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data,
                                       "method" = glm_method))
      stats::glm
    },
    "lm" = {
      lm_method <- if (fit) "qr" else "model.frame"
      model_args <- c(model_args, list("formula" = formula,
                                       "data" = data,
                                       "method" = lm_method))
      stats::lm
    },
    "glmer" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data))
      lme4::glmer
    },
    "lmer" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "data" = data))
      model_args$REML <- FALSE
      lme4::lmer
    },
    "gam" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data,
                                       "fit" = fit))
      model_args$method <- "REML"
      mgcv::gam
    },
    "nlme" = {
      model_args <- c(model_args, list("model" = formula,
                                       "data" = data))
      nlme::nlme
    },
    "nls" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "data" = data))
      model_args$model <- TRUE
      stats::nls
    }
  )

  do.call(model_fun, model_args)
}
