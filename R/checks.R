#' Check that model type is supported by package
#'
#' Function to check whether model type is among supported model types.
#' @param type
#'  The model type to be checked
check_model_type <- function(type) {
  all_options <- c("lm", "glm", "lmer", "glmer", "nlme", "gam", "nls")
  if (!type %in% all_options) {
    stop(
      sprintf(
        paste("Your chosen model type %s is not supported by LazyModeler.",
              "Valid options are: %s.",
              collapse = " "),
        type,
        paste(all_options, collapse = ", ")
      )
    )
  }
}

#' Check validity of model family
#'
#' Checks whether model family is supported and matches response data of
#'  given formula
#' @param family
#'  The chosen model family. Can be NULL if automatic == TRUE. Default: NULL
#' @param automatic
#'  Whether to automatically choose an appropriate family. Default: NULL
#' @param data
#'  The underlying data to use when determining an appropriate model family.
#'    Default: NULL
#' @param lhs
#'  The left-hand-side of the formula used for downstream model creation.
#'    Default: NULL
#' @returns
#'  An appropriate model family for downstream model creation
check_model_family <- function(family = NULL,
                               automatic = TRUE,
                               data = NULL,
                               lhs = NULL) {
  valid_families_info <- determine_model_family(data, lhs)
  valid_families <- valid_families_info$valid_families

  if (automatic) {
    if (length(valid_families) > 1) {
      stop(
        sprintf(
          paste("Your response variable allows for the following",
                "distributions: %s. %s Carefully consider all options, then",
                "rerun optimize_model() with your preferred distribution.",
                collapse = " "),
          paste(valid_families, collapse = " and "),
          valid_families_info$notice
        ),
        call. = FALSE
      )
    } else {
      family <- valid_families[1]

      message(sprintf(
        paste("%s Continuing with %s distribution. If you want to use",
              "a different distribution, make sure your data is formatted",
              "correctly, then rerun optimize_model() with your",
              "preferred distribution.", collapse = " "),
        valid_families_info$notice,
        family
      ))
    }
  } else if (!get_family_character(family) %in% valid_families) {
    warning(
      sprintf(
        paste("Chosen distribution '%s' does not match response values.",
              "Would recommend %s. %s Please check.", collapse = " "),
        family,
        paste(valid_families, collapse = ", or "),
        valid_families_info$notice
      )
    )
  }

  family
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

#' Check for data type of response
#'
#' Checks whether formatting of response column matches data
#' @param response_col
#'  The response column name
#' @param response_data
#'  The response data points
check_response_data_format <- function(response_col, response_data) {
  if (all(is.finite(response_data)) && all(response_data %% 1 == 0)) {
    warning(
      sprintf(
        paste(
          "Your response column %s looks to be of type integer",
          "but is formatted as numeric.",
          "We are going to treat it as numeric. If your response is of",
          "type integer, please correct the column's formatting and restart.",
          collapse = " "
        ),
        response_col
      )
    )
  }
}

#' Checks model family in case of [cbind()]
#'
#' Checks which model family matches data when [cbind()]
#'  was specified in formula
#' @param data
#'  The underlying data to use when determining an appropriate model family
#' @param pasted_response
#'  Response as character
#' @returns A list containing all valid families and a notice explaining
#'  the family choice
check_cbind_model_family <- function(data, pasted_response) {
  if (length(pasted_response) > 3) {
    warning(paste(
      "It seems like you have specified more than two variables",
      "within cbind(). This is not supported. We will only consider the first",
      "two variables. If you wish to correct your formula, please stop and",
      "rerun LazyModeler.",
      collapse = " "
    ))
  }

  col1 <- pasted_response[[2]]
  col2 <- pasted_response[[3]]

  family1 <- determine_model_family(data, str2lang(col1))$family
  family2 <- determine_model_family(data, str2lang(col2))$family

  if ("quasibinomial" %in% family1 || "quasibinomial" %in% family2) {
    valid_families <- c("quasibinomial")
    notice <- paste("Your left-hand-side appears to be a grouped binomial",
                    "outcome and your columns fit a quasibinomial model",
                    "family. If you intended to fit a multivariate linear",
                    "model, unfortunately, we do not support that right now.",
                    collape = " ")
  } else {
    valid_families <- c("binomial")
    notice <- paste("Your left-hand-side appears to be a grouped binomial",
                    "outcome. Please also check for overdispersion;",
                    "if present, quasibinomial may be more appropriate.",
                    "If you intended to fit a multivariate linear model,",
                    "unfortunately, we do not support that right now.",
                    collape = " ")
  }

  list(valid_families = valid_families,
       notice = notice)
}

#' Formula formatting and integrity check
#'
#' Formats formula for downstream process and checks for main effects
#' @param formula
#'  The formula to check
#' @param data
#'  The underlying data for model creation
#' @param add_main
#'  Whether to add main effect. Default: TRUE
#' @returns The updated formula with added main effects if add_main == TRUE
check_formula <- function(formula, data, add_main = TRUE) {
  formula <- stats::formula(stats::terms(stats::as.formula(formula),
                                         data = data))

  terms <- attr(stats::terms.formula(formula), "term.labels")
  interactions <- extract_interactions(terms)
  main_effects <- interactions$main_effect
  missing_main_effects <- setdiff(main_effects, terms)

  if (length(missing_main_effects) > 0) {
    if (!add_main) {
      warning(sprintf(
        paste("We have noticed that certain main effects are missing from",
              "your formula. Main effects are required for hierarchical",
              "integrity and correct coefficient interpretation.",
              "Please add the following main effects to your formula: %s",
              collapse = " "),
        paste(missing_main_effects, collapse = ", ")
      ))
    } else {
      warning(sprintf(
        paste("We have noticed that certain main effects are missing from",
              "your formula. Main effects are required for hierarchical",
              "integrity and correct coefficient interpretation.",
              "We will add the following main effects to your formula: %s",
              collapse = " "),
        paste(missing_main_effects, collapse = ", ")
      ))
      for (main_effect in missing_main_effects) {
        d <- paste(". ~ . +", main_effect)
        formula <- stats::update(stats::as.formula(formula), d)
      }
    }
  }

  formula
}

#' Check whether interactions can be plotted
#'
#' Checks whether given interactions are supported for plotting. Not supported
#'  are interactions of more than two variables and interactions between
#'  factors.
#' @param interactions
#'  The interactions to check
#' @returns A filtered interactions dataframe supported for plotting
check_plot_interactions <- function(interactions) {
  interactions <- interactions |>
    dplyr::group_by(.data$interaction) |>
    dplyr::mutate(not_twoway = dplyr::n() > 2,
                  only_cat = all(!.data$is_numeric))
  if (any(interactions$not_twoway)) {
    warning(paste("We have recognized interactions with more than two",
                  "variables. These interactions will not be plotted.",
                  collapse = " "))
  }
  if (any(interactions$only_cat)) {
    warning(paste("We have recognized interactions between two factors.",
                  "We currently do not support plotting of these types of",
                  "interactions.",
                  collapse = " "))
  }

  interactions[(!interactions$not_twoway) & (!interactions$only_cat), ]
}
