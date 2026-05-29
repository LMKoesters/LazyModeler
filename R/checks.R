check_model_type <- function(model_type) {
  if (!model_type %in% c("lm", "glm", "lmer", "glmer", "gam", "nlme", "nls")) {
    stop(
      sprintf(
        paste("Unknown model type %s. Please choose either:",
              "`lm`, `glm`, `lmer`, `glmer`, `gam`, `nlme`, `nls`",
              collapse = " "),
        model_type
      ),
      .call = FALSE
    )
  }
}

check_model_family <- function(family) {
  # TODO
  # determine family etc
  # keep in mind that family can be string or function
}