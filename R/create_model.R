create_model <- function(formula,
                         data,
                         model_type,
                         family = gaussian,
                         model_args = list(),
                         logger = NULL) {
  model_fun <- switch(
    as.character(model_type),
    "glm" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data))
      stats::glm
      },
    "lm" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "data" = data))
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
      lmerTest::lmer
    },
    "gam" = {
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data))
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
      stats::nls
    }
  )
  
  regression_model <- withCallingHandlers(
    do.call(model_fun, model_args),
    warning = function(w) {
      if (grepl("non-integer #successes in a binomial glm", w$message)) {
        logger <- append(logger, w$message)
        tryInvokeRestart("muffleWarning")
      } else {
        warning(w$message)
      }
    }
  )
  
  regression_model
}