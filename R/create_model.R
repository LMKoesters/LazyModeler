create_model <- function(formula,
                         family,
                         data,
                         model_type,
                         model_args,
                         logger) {
  model_fun %<-% switch(
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
      # TODO glmer -> c("formula", "family", "data")
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
      # TODO gam -> c("formula", "family", "data")
      model_args <- c(model_args, list("formula" = formula,
                                       "family" = family,
                                       "data" = data))
      mgcv::gam
    },
    "nlme" = {
      # TODO nlme -> c("model", "data")
      # also important to keep in mind: fixed, random
    },
    "nls" = {
      # TODO nls -> c("formula", "data")
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