check_model_type <- function(model_type) {
  # TODO this won't work if model_type is something like gaussian()
  # Need to convert to string
  
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
}