x <- function() {
  if (model_type == "lmer") {
    # TODO REML kann weg
    model_args <- list(formula = term, data = df, REML = FALSE)
  } else if (
    (model_type == "lm") ||
    (model_type == "lmer") ||
    (model_type == "nls")
  ) {
    model_args <- list(formula = term, data = df)
  } else if (model_type == "nlme") {
    model_args <- list(data = df)
  } else {
    model_args <- list(formula = term, data = df, family = model_family)
  }
  model_args <- c(model_args, list(...))
  
  
  if (model_type == "glm") {
    regression_model <- withCallingHandlers(
      do.call(stats::glm, model_args),
      warning = function(w) {
        if (grepl("non-integer #successes in a binomial glm", w$message)) {
          errors <- append(errors, w$message)
          tryInvokeRestart("muffleWarning")
        } else {
          warning(w$message)
        }
      }
    )
  } else if (model_type == "lm") {
    regression_model <- do.call(stats::lm, model_args)
  } else if (model_type == "glmer") {
    regression_model <- try(do.call(lme4::glmer, model_args), silent = TRUE)
  } else if (model_type == "lmer") {
    regression_model <- do.call(lmerTest::lmer, model_args)
  } else if (model_type == "gam") {
    regression_model <- do.call(mgcv::gam, model_args)
  } else if (model_type == "nlme") {
    parameter_names <- names(model_args)
    parameter_names[parameter_names == "non_linear"] <- "model"
    names(model_args) <- parameter_names
    regression_model <- do.call(nlme::nlme, model_args)
  } else if (model_type == "nls") {
    regression_model <- do.call(stats::nls, model_args)
  } else {
    stop(
      sprintf(
        "Unknown model type %s. Please choose either (g/n)(l/a)m(er).",
        model_type
      )
    )
  }
}