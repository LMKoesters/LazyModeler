#' Check that model type is supported by package
#' 
#' Function to check whether model type is among supported model types.
#' @param model_type
#'  The model type to be checked.
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

#' Check that threshold is between 0-1
#'
#' Function to check whether threshold is within appropriate range of (0,1)
#' @param threshold
#'  Threshold value to check
check_correlation_threshold <- function(threshold) {
  if (threshold < 0 || threshold > 1) {
    stop(
      sprintf(
        paste("Your threshold %f is outside of range (0,1).",
          "Please choose an appropriate threshold and run again."
        ),
        threshold
      )
    )
  }
}
