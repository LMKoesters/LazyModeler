#' This is data to be included in my package
#'
#' @name dataset_info
#' @docType data
#' @author Lara Kösters \email{lkoesters@@bgc-jena.mpg.de}
#' @author Kevin Karbstein \email{kkarb@@bgc-jena.mpg.de}
#' @keywords data
NULL

#' Autocorrelation identification and removal
#'
#' Identify autocorrelated variables.
#' Automatically remove autocorrelated variables if specified
#' @param df
#'  Dataframe with response and predictors as columns
#' @param coefficients
#'  List of predictors sorted by their relevance (most to least)
#' @param ...
#'  Arguments given directly to 'cor' method
#' @param automatic_removal
#'  Determines whether to automatically remove autocorrelations
#' @param autocorrelation_threshold
#'  Threshold at which two variables are considered autocorrelated
#' @param correlation_method
#'  Method used for correlation calculation
#' @return
#'  Named list with a) a vector containing all removed predictors
#'  (empty if none were removed), and
#'  b) a dataframe containing autocorrelations and
#'  information on removed variables
remove_autocorrelations <- function(
  df,
  coefficients,
  ...,
  automatic_removal = TRUE,
  autocorrelation_threshold = 0.7,
  correlation_method = "pearson"
) {
  if (autocorrelation_threshold >= 1.0 || autocorrelation_threshold <= 0.0) {
    stop(
      stringr::str_interp(
        "Values for the autocorrelation threshold should be between 0 and 1."
      )
    )
  }

  removed_coefficients <- vector()
  cor_args <- list(x = df[, coefficients], method = correlation_method)
  cor_args <- c(cor_args, list(...))
  correlations <- as.data.frame(do.call(stats::cor, cor_args))

  # compute p-values
  correlation_p_values <- as.data.frame(
    corrplot::cor.mtest(df[, coefficients])$p
  )

  correlations_l <- correlations %>%
    tibble::rownames_to_column(var = "coefficientA") %>%
    tidyr::pivot_longer(!"coefficientA",
      names_to = "coefficientB",
      values_to = "correlation"
    )

  correlations_p_values_l <- correlation_p_values %>%
    tibble::rownames_to_column(var = "coefficientA") %>%
    tidyr::pivot_longer(!"coefficientA",
      names_to = "coefficientB",
      values_to = "p_value"
    )

  correlations_complete <- merge(correlations_l,
    correlations_p_values_l,
    by = c("coefficientA", "coefficientB")
  )

  correlations_complete <- correlations_complete %>%
    dplyr::filter(.data$coefficientA != .data$coefficientB) %>%
    dplyr::mutate(
      col1 = pmin(.data$coefficientA, .data$coefficientB),
      col2 = pmax(.data$coefficientA, .data$coefficientB),
      comparison = paste(.data$col1, .data$col2)
    ) %>%
    dplyr::distinct(.data$comparison, .keep_all = TRUE) %>%
    dplyr::select(!tidyr::any_of(c("col1", "col2", "comparison")))

  autocorrelations <- correlations_complete %>%
    dplyr::filter(
      (
        (.data$correlation >= autocorrelation_threshold) |
          (.data$correlation <= -autocorrelation_threshold)
      ) &
        (.data$p_value < 0.05)
    ) %>%
    tibble::add_column(note = NA)

  if (nrow(autocorrelations) > 0) {
    if (automatic_removal) {
      coefficients_df <- data.frame(
        idx = seq_along(coefficients),
        coefficient = coefficients
      )

      autocorrelations <- merge(autocorrelations,
        coefficients_df[c("coefficient", "idx")],
        by.x = "coefficientA", by.y = "coefficient"
      )
      autocorrelations <- merge(autocorrelations,
        coefficients_df[c("coefficient", "idx")],
        by.x = "coefficientB", by.y = "coefficient", suffixes = c("1", "2")
      )
      autocorrelations <- autocorrelations %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          idx_smaller = sort(c(.data$idx1, .data$idx2))[[1]],
          idx_bigger = sort(c(.data$idx1, .data$idx2))[[2]],
          comb = paste(.data$idx_smaller, .data$idx_bigger, sep = "|")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(.data$idx_bigger))

      for (i in seq_len(nrow(autocorrelations))) {
        # C is the least important and part of comparison
        # B is the more important of comparison
        # A is more important than B

        # check for following constellation: A==B, B==C, A!=C
        autocor_row <- autocorrelations[i, ]
        c <- autocor_row$idx_bigger
        b <- autocor_row$idx_smaller

        b_to_a <- autocorrelations[autocorrelations$idx_bigger == b, ]
        if (nrow(b_to_a) > 0) { # check if there's at least 1 A
          for (a in b_to_a$idx_smaller) {
            # for each variable that is correlated to and more important than B
            coefficient_b <- coefficients_df[b, "coefficient"]
            coefficient_c <- coefficients_df[c, "coefficient"]

            a_b_c <- autocorrelations[
              (autocorrelations$idx_bigger == autocor_row$idx_bigger) &
                (autocorrelations$idx_smaller == a),
            ]
            if (nrow(a_b_c) == 0) {
              # A!=C but A==B and B==C: remove B
              coefficient_a <- coefficients_df[a, "coefficient"]
              coefficients <- coefficients[!coefficients == coefficient_b]
              removed_coefficients <- append(
                removed_coefficients,
                coefficient_b
              )
              autocorrelations[i, "note"] <- sprintf(
                "%s!=%s, but %s==%s and %s==%s; removed %s",
                coefficient_a,
                coefficient_c,
                coefficient_a,
                coefficient_b,
                coefficient_b,
                coefficient_c,
                coefficient_b
              )
            } else {
              coefficients <- coefficients[!coefficients == coefficient_c]
              removed_coefficients <- append(
                removed_coefficients,
                coefficient_c
              )
              autocorrelations[i, "note"] <- stringr::str_interp(
                "removed ${coefficient_c}"
              )
            }
          }
        } else {
          coefficient_c <- coefficients_df[c, "coefficient"]
          coefficients <- coefficients[!coefficients == coefficient_c]
          removed_coefficients <- append(removed_coefficients, coefficient_c)
          autocorrelations[i, "note"] <- stringr::str_interp(
            "removed ${coefficient_c}"
          )
        }
      }
    } else {
      f <- sprintf(
        "autocorrelations_%s.tsv",
        format(Sys.time(), "%m%d%Y_%H%M")
      )
      utils::write.table(autocorrelations, f,
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      stop(
        stringr::str_interp(
          "Some of your variables are autocorrelated. Check ${f} for more info"
        )
      )
    }
  }

  list(
    "removed_predictors" = removed_coefficients,
    "autocorrelations" = autocorrelations[
      ,
      c("coefficientA", "coefficientB", "correlation", "p_value", "note")
    ]
  )
}

#' Model variables categorical check
#'
#' Helper method used for expanding the model summary
#' by adding information on whether or not model variables are categorical
#' @param df
#'  Dataframe with response and predictors as columns
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param col
#'  Column to be used for check. Default is "predictor"
#' @param val_else
#'  Value to be returned by default. Can be a dataframe column.
#' @return
#'  Language object with if/else statements to be applied to
#'  model summary with information on whether model variables are categorical
setup_categorical_check <- function(
  df,
  categorical_vars,
  col = quote(predictor),
  val_else = quote(NA_character_)
) {
  categorical_check <- c()
  for (cat_var in categorical_vars) {
    cat_vals <- as.character(unique(df[[cat_var]]))
    categorical_check <- append(categorical_check, purrr::map(
      cat_vals,
      ~ quo(stringr::str_detect(
        !!col,
        stringr::regex(paste0(!!cat_var, !!.x, "(?![a-z0-9_\\-#])"))
      ) ~ !!cat_var)
    ))
  }

  categorical_check <- c(
    categorical_check,
    rlang::quo(TRUE ~ as.character(!!val_else))
  )

  categorical_check
}

#' Add model assessments
#'
#' Add evaluations of model given model and
#' evaluation methods (anova/aic/aicc/bic)
#' @param regression_model
#'  Regression model that is assessed
#' @param evaluation_methods
#'  List of methods to use for model evaluation.
#'  Options are: anova, aic, aicc, and bic
#' @param is_lmer
#'  Boolean stating whether model is of type lmer
#' @return
#'  Regression model with added assessments based on assessment methods provided
add_assessments <- function(
  regression_model,
  evaluation_methods,
  is_lmer = FALSE
) {
  if ("aic" %in% evaluation_methods) {
    if (is_lmer) {
      attr(regression_model, "aic") <- stats::AIC(regression_model)
    } else {
      regression_model$aic <- stats::AIC(regression_model)
    }
  }
  if ("aicc" %in% evaluation_methods) {
    if (is_lmer) {
      attr(regression_model, "aicc") <- MuMIn::AICc(regression_model)
    } else {
      regression_model$aicc <- MuMIn::AICc(regression_model)
    }
  }
  if ("bic" %in% evaluation_methods) {
    if (is_lmer) {
      attr(regression_model, "bic") <- stats::BIC(regression_model)
    } else {
      regression_model$bic <- stats::BIC(regression_model)
    }
  }

  regression_model
}

#' Regression model was improved
#'
#' Determine if model 1 better fits the data than model 2
#' given two models and evaluation methods (anova/aic/aicc/bic)
#' @param regression_model
#'  Regression model 1
#' @param old_regression_model
#'  Regression model 2
#' @param evaluation_methods
#'  List of methods to use for model evaluation.
#'  Options are: anova, aic, aicc, and bic
#' @param direction
#'  Mode of stepwise model improvement.
#'  Either 'forward' (i.e., forward selection),
#'  or 'backward' (i.e., backward simplification).
#' @param model_type
#'  The model to be used. Options: (g)lm, (g)lmer, nlmer, and gam.
#'  Default: glm.
#' @return
#'  Boolean; TRUE if model 1 better fits the data than model 2
model_improved <- function(
  regression_model,
  old_regression_model,
  evaluation_methods,
  direction,
  model_type
) {
  if ("anova" %in% evaluation_methods) {
    if (model_type == "gam") {
      anova_res <- mgcv::anova.gam(regression_model, old_regression_model)
    } else {
      anova_res <- stats::anova(regression_model, old_regression_model)
    }

    p_col <- colnames(anova_res)[[grep("Pr\\(", colnames(anova_res))]]

    if ((direction == "backward") && (anova_res[2, p_col] < .05)) {
      return(FALSE)
    } else if ((direction == "forward") && (anova_res[2, p_col] > .05)) {
      return(FALSE)
    }
  }

  is_lmer <- model_type == "lmer" || model_type == "glmer"

  if ("aic" %in% evaluation_methods) {
    if (is_lmer) {
      aic_new <- attr(regression_model, "aic")
      aic_old <- attr(old_regression_model, "aic")
    } else {
      aic_new <- regression_model$aic
      aic_old <- old_regression_model$aic
    }

    if (aic_new >= aic_old) {
      return(FALSE)
    }
  }
  if ("aicc" %in% evaluation_methods) {
    if (is_lmer) {
      aicc_new <- attr(regression_model, "aicc")
      aicc_old <- attr(old_regression_model, "aicc")
    } else {
      aicc_new <- regression_model$aicc
      aicc_old <- old_regression_model$aicc
    }

    if (aicc_new >= aicc_old) {
      return(FALSE)
    }
  }
  if ("bic" %in% evaluation_methods) {
    if (is_lmer) {
      bic_new <- attr(regression_model, "bic")
      bic_old <- attr(old_regression_model, "bic")
    } else {
      bic_new <- regression_model$bic
      bic_old <- old_regression_model$bic
    }

    if (bic_new >= bic_old) {
      return(FALSE)
    }
  }

  TRUE
}


#' Generate regression model
#'
#' Helper method to generate a regression model
#' based on a given model type, family and formula
#' @param df
#'  Dataframe with response and predictors as columns
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam. Default: glm.
#' @param term
#'  The formula to be used with the model.
#'  Can be either quote() or formula()
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options.
#' @param ...
#'  Arguments given to model call
#' @return
#'  List with regression model and used model arguments
generate_regression_model <- function(df, model_type, term, model_family, ...) {
  errors <- c()

  if (model_type == "lmer") {
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
    regression_model <- withCallingHandlers(do.call(stats::glm, model_args),
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

  if ("try-error" %in% class(regression_model)) {
    e <- attr(regression_model, "condition")
    if (grepl("converge", e$message)) {
      warning(
        stringr::str_interp(
          "${e$message}. Please change the random effect or x-variable set."
        )
      )
      errors <- append(
        errors,
        stringr::str_interp(
          "${e$message}. Please change the random effect or x-variable set."
        )
      )
      return(list(
        regression_model = regression_model,
        model_args = model_args,
        errors = errors
      ))
    }
  }

  list(
    regression_model = regression_model,
    model_args = model_args,
    errors = errors
  )
}

#' Forward model selection
#'
#' Forward model selection
#' @param df
#'  Dataframe with response and predictors as columns
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam. Default: glm.
#' @param term
#'  The formula to be used with the model.
#'  Can be either quote() or formula().
#' @param evaluation_methods
#'  Methods to be used for model evaluation.
#'  Options: anova, aic, aicc, bic. Default: anova.
#' @param ...
#'  Parameters to be directly used with model call
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param trace
#'  Store and return model selection history. Default: FALSE.
#' @param omit_na
#'  Either 'overall' or 'stepwise'.
#'  If 'overall', NAs are removed before modeling.
#'  If 'stepwise', NAs are removed per step based
#'  on the variables in the current formula.
#' @return
#'  List containing the final regression model,
#'  the significant and marginally significant model
#'  variables and the selection history if trace is TRUE
forward_selection <- function(
  df,
  model_type,
  term,
  evaluation_methods,
  categorical_vars,
  ...,
  model_family = "gaussian",
  trace = FALSE,
  omit_na = "overall"
) {
  errors <- c()
  single_vars <- c()
  interactions <- c()

  split_frm <- find_call(term, return = "fixed")

  if (omit_na == "overall") df <- remove_nas(df, term, model_type)

  # remove random effect from single_vars because it included in start term
  random_effects <- lapply(reformulas::findbars(term), function(y) deparse(y))
  if (length(random_effects) > 0) {
    d <- paste(". ~ (", paste(random_effects, sep = " + "), ")")
  } else {
    d <- ". ~ 1"
  }

  for (expr in split_frm) {
    if (stringr::str_detect(deparse(expr), ":")) {
      interactions <- append(interactions, deparse(expr))
    } else if (!(deparse(expr) %in% random_effects)) {
      single_vars <- append(single_vars, deparse(expr))
    }
  }

  frm <- stats::update(stats::as.formula(term), d)
  regression_model <- generate_regression_model(
    df,
    model_type,
    frm,
    model_family,
    ...
  )
  model_args <- regression_model$model_args
  regression_model <- regression_model$regression_model

  if (trace) history <- list()

  blank_start <- TRUE

  for (vars in list(single_vars, interactions)) {
    simplify <- TRUE
    while (simplify) {
      res <- mo_step(
        df,
        regression_model,
        vars,
        "forward",
        evaluation_methods,
        categorical_vars,
        model_type,
        blank_start = blank_start,
        stepwise_omit = (omit_na == "stepwise"),
        model_args = model_args
      )
      errors <- append(errors, res$errors)
      if (blank_start) blank_start <- FALSE
      simplify <- res$simplify
      regression_model <- res$regression_model
      vars <- res$vars
      if (trace && simplify) {
        history[[
          stringr::str_interp("${length(history)+1}_${res$var}")
        ]] <- regression_model
      }
    }
  }

  if (model_type == "gam") {
    model_sum <- summary(regression_model)$p.table
  } else {
    model_sum <- stats::coef(summary(regression_model))
  }

  model_summary <- expand_model_summary(
    model_summary = model_sum,
    term = stats::formula(regression_model),
    categorical_vars = categorical_vars,
    df,
  )

  main_effects <- c(model_summary$main_effect1, model_summary$main_effect2)
  main_effects <- main_effects[!is.na(main_effects)]

  for (main_effect in main_effects) {
    if (!(main_effect %in% model_summary$predictor)) {
      if (model_type == "lmer") {
        # update formula
        model_args$formula <- stats::update(
          stats::formula(regression_model),
          paste(". ~ . +", main_effect)
        )

        # update data
        if (omit_na == "stepwise") {
          model_args$data <- remove_nas(
            model_args$data,
            model_args$formula,
            model_type
          )
        }

        regression_model <- do.call(lmerTest::lmer, model_args)
      } else {
        regression_model <- withCallingHandlers(
          stats::update(regression_model, paste("~ . +", main_effect)),
          warning = function(w) {
            if (grepl("non-integer #successes in a binomial glm", w$message)) {
              errors <- append(errors, w$message)
              tryInvokeRestart("muffleWarning")
            } else {
              warning(w$message)
            }
          }
        )
      }
      if (trace) {
        history[[
          stringr::str_interp("${length(history)+1}_${main_effect}")
        ]] <- regression_model
      }
    }
  }

  out <- list(
    final_model = regression_model,
    significant_variables = model_summary[
      model_summary$p_value < 0.05,
      "predictor"
    ],
    marginally_significant_variables = model_summary[
      (model_summary$p_value < 0.1) & (model_summary$p_value >= 0.05),
      "predictor"
    ]
  )

  if (trace) out$history <- history
  out$errors <- errors

  out
}

#' Take step in model simplification/selection
#'
#' Used to either remove or add a variable to
#'  an already existing regression model
#' @param df
#'  Dataframe with response and predictors as columns
#' @param regression_model
#'  The current regression model
#' @param vars
#'  List of model variables that are yet to be added/removed
#' @param direction
#'  Mode of stepwise model improvement.
#'  Either 'forward' (i.e., forward selection),
#'  or 'backward' (i.e., backward simplification)
#' @param evaluation_methods
#'  Methods to be used for model evaluation.
#'  Options: anova, aic, aicc, bic, default: anova.
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam. Default: glm.
#' @param blank_start
#'  Used for first run in forward model selection.
#'  Creates empty model.
#' @param stepwise_omit
#'  If TRUE, NAs are removed per step based on
#'  the variables in the current formula
#' @param model_args
#'  Arguments given directly to model call
#' @return
#'  List containing a) information on whether model was simplified/expanded,
#'  b) the (new) current regression model,
#'  c) updated list of variables to be added/removed,
#'  d) added/removed variable, and
#'  e) summary of the model
mo_step <- function(
  df,
  regression_model,
  vars,
  direction,
  evaluation_methods,
  categorical_vars,
  model_type,
  blank_start = FALSE,
  stepwise_omit = FALSE,
  model_args = NA
) {
  errors <- c()
  is_lmer <- model_type == "lmer" || model_type == "glmer"
  final_model <- NA
  adjusted_var <- NA
  model_summary <- NA
  min_p_value <- Inf
  reason <- NA
  simplified <- FALSE

  interactions <- unique(vars[stringr::str_detect(vars, ":")])

  for (var in vars) { # loop important for forward simplification
    # add or remove var from model
    if (direction == "backward") d <- ". ~ . -"

    # check for interactions when removing single vars and remove them as well
    var_m <- var
    if (direction == "backward" && !stringr::str_detect(var, ":")) {
      for (interaction in interactions) {
        if (
          stringr::str_detect(
            interaction,
            stringr::fixed(stringr::str_interp("${var}:"))
          ) ||
            stringr::str_detect(
              interaction,
              stringr::fixed(stringr::str_interp(":${var}"))
            )
        ) {
          var_m <- paste(var_m, "-", interaction)
        }
      }
    }

    # update model
    if (direction == "forward") {
      # stepwise: more restrictive model is regression_model_updated
      if (is_lmer) {
        # update formula
        model_args$formula <- stats::update(
          stats::formula(regression_model),
          paste(". ~ ", var_m, " + .")
        )

        # update data
        if (stepwise_omit) {
          model_args$data <- remove_nas(
            model_args$data,
            model_args$formula,
            model_type
          )
        }

        regression_model_updated <- do.call(lmerTest::lmer, model_args)

        if (stepwise_omit) {
          model_args$formula <- stats::formula(regression_model)
          regression_model <- do.call(lmerTest::lmer, model_args)
        }
      } else {
        regression_model_updated <- withCallingHandlers(
          stats::update(
            regression_model,
            paste(". ~ ", var_m, " + ."),
            data = df
          ),
          warning = function(w) {
            if (grepl("non-integer #successes in a binomial glm", w$message)) {
              errors <- append(errors, w$message)
              tryInvokeRestart("muffleWarning")
            } else {
              warning(w$message)
            }
          }
        )

        # if interaction, make sure order is same as before
        if (stringr::str_detect(var, ":")) {
          if (!stringr::str_detect(
            Reduce(
              paste,
              deparse(regression_model_updated$formula)
            ),
            var
          )) {
            vars <- vars[vars != var]
            var_spl <- as.vector(unlist(stringr::str_split(var, ":")))
            var <- paste0(var_spl[[2]], ":", var_spl[[1]])
            vars <- c(vars, var)
          }
        }

        if (stepwise_omit) {
          tryCatch(
            expr = {
              regression_model <- stats::update(
                regression_model,
                ". ~ .",
                data = regression_model_updated$model
              )
            },
            error = function(e) {
              regression_model <- stats::update(
                regression_model,
                ". ~ .",
                data = regression_model_updated$data
              )
            }
          )
        }
      }

      if (model_type == "gam") {
        model_sum <- summary(regression_model_updated)$p.table
      } else {
        model_sum <- stats::coef(summary(regression_model_updated))
      }
    } else {
      # stepwise: more restrictive model is regression_model
      if (is_lmer) {
        # update formula
        model_args$formula <- stats::update(
          stats::formula(regression_model),
          paste(d, var_m)
        )

        # update data
        if (stepwise_omit) {
          model_args$data <- remove_nas(
            model_args$data,
            stats::formula(regression_model),
            model_type
          )
        }

        if (model_type == "glmer") {
          regression_model_updated <- try(
            do.call(lme4::glmer, model_args),
            silent = TRUE
          )
        } else {
          regression_model_updated <- do.call(lmerTest::lmer, model_args)
        }

        if ("try-error" %in% class(regression_model_updated)) {
          e <- attr(regression_model_updated, "condition")
          if (grepl("converge", e$message)) {
            warning(
              sprintf(
                "%s. Please change the random effect or x-variable set.",
                e$message
              )
            )
            errors <- append(
              errors,
              sprintf(
                "%s. Please change the random effect or x-variable set.",
                e$message
              )
            )
            return(list(
              simplify = simplified,
              regression_model = regression_model,
              vars = vars,
              model_summary = model_summary,
              var = adjusted_var,
              reason = reason,
              errors = errors
            ))
          }
        }
      } else {
        regression_model_updated <- withCallingHandlers(
          stats::update(
            regression_model,
            paste(d, var_m),
            data = regression_model$model
          ),
          warning = function(w) {
            if (grepl("non-integer #successes in a binomial glm", w$message)) {
              errors <- append(errors, w$message)
              tryInvokeRestart("muffleWarning")
            } else {
              warning(w$message)
            }
          }
        )
      }

      if (model_type == "gam") {
        model_sum <- summary(regression_model_updated)$p.table
      } else {
        model_sum <- stats::coef(summary(regression_model_updated))
      }
    }

    model_summary <- expand_model_summary(
      model_summary = model_sum,
      term = stats::formula(regression_model_updated),
      categorical_vars = categorical_vars,
      df = df,
    )

    regression_model_updated <- add_assessments(
      regression_model_updated,
      evaluation_methods,
      is_lmer = model_type == "lmer"
    )

    # first check if variables are significant
    if (!blank_start) {
      model_has_improved <- model_improved(
        regression_model_updated,
        regression_model,
        evaluation_methods,
        direction,
        model_type
      )
    }

    if (((direction == "forward") &&
      (model_summary[1, "smallest_p_value"] > 0.1)) ||
      (!blank_start && !model_has_improved)) {
      reason <- "the p-value was not significant"
      next
    }

    this_p_value <- NA
    if (direction == "forward") {
      # get smallest p-value of added variable
      #' if categorical: f1low, f1high, f1shallow,... -> f1 smallest p-value
      #' if part of interaction: x1 -> smallest p-value across x1;x1:x2;...
      this_p_value <- model_summary[
        (model_summary$predictor == var) |
          ((model_summary$main_effect1 == var) &
            (is.na(model_summary$main_effect2))),
        "smallest_p_value"
      ]
      if (length(this_p_value) > 1) this_p_value <- this_p_value[[1]]
    }

    # check if variable should be added/removed
    # this is either the case with backward simplification
    # or when the p-value of the added variable is the smallest encountered
    # across all add-able variables
    if ((direction == "backward") || (this_p_value < min_p_value)) {
      final_model <- regression_model_updated
      adjusted_var <- var
      min_p_value <- this_p_value
      simplified <- TRUE
    } else {
      reason <- "the p-value was not significant"
    }

    if (direction == "backward") break
  }

  if (simplified) {
    list(
      simplify = simplified,
      regression_model = final_model,
      vars = vars[vars != adjusted_var],
      var = adjusted_var,
      model_summary = model_summary,
      errors = errors
    )
  } else {
    list(
      simplify = simplified,
      regression_model = regression_model,
      vars = vars,
      model_summary = model_summary,
      var = adjusted_var,
      reason = reason,
      errors = errors
    )
  }
}


#' Model simplification
#'
#' Simplify or expand a regression model given the data and model parameters
#' @param df
#'  Dataframe with response and predictors as columns
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam.
#' @param term
#'  The formula to be used with the model.
#'  Can be either quote() or formula()
#' @param evaluation_methods
#'  Methods to be used for model evaluation.
#'  Options: anova, aic, aicc, bic
#' @param ...
#'  Parameters to be directly used with model call
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param direction
#'  Mode of stepwise model improvement.
#'  Either 'forward' (i.e., forward selection),
#'  or 'backward' (i.e., backward simplification),
#'  or 'both' (default: both).
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param backward_simplify_model
#'  If FALSE, the model and information on significant and
#'  marginally significant variables are returned without any improvements
#' @param trace
#'  Store and return model selection history.
#'  Default: FALSE
#' @param omit_na
#'  Either 'overall' or 'stepwise'.
#'  If 'overall', NAs are removed before modeling.
#'  If 'stepwise', NAs are removed per step based on
#'  the variables in the current formula
#' @return
#'  List of results from forward and/or backward model selection.
#'  Backward/forward result will have the following structure:
#'  List containing the final regression model, the significant and
#'  marginally significant model variables and, if trace is TRUE,
#'  the selection history
#' @examples
#' # setup
#' data("plants")
#'
#' # generate a glm model with the provided term and simplify it
#' # by applying backward simplification
#' simplified_model_info <- simplify_model(plants,
#'   "glm",
#'   quote(sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     I(isothermality^2) +
#'     habitat +
#'     ploidy +
#'     solar_radiation:annual_mean_temperature +
#'     solar_radiation:isothermality +
#'     annual_mean_temperature:isothermality),
#'   c("anova"),
#'   model_family = "quasibinomial",
#'   direction = "backward",
#'   categorical_vars = c("habitat", "ploidy"),
#'   backward_simplify_model = TRUE,
#'   trace = TRUE,
#'   omit_na = "overall"
#' )
#' @export
simplify_model <- function(
  df,
  model_type,
  term,
  evaluation_methods,
  ...,
  model_family = "gaussian",
  direction = "both",
  categorical_vars = c(),
  backward_simplify_model = TRUE,
  trace = FALSE,
  omit_na = "overall"
) {
  if (!(omit_na %in% c("overall", "stepwise"))) {
    stop(
      sprintf(
        "Unknown option %s for option omit_na. Please choose a valid option.",
        omit_na
      )
    )
  }

  models <- list()

  if (model_type == "nlme") {
    generated_model <- generate_regression_model(
      df,
      model_type,
      term,
      model_family,
      ...
    )
    regression_model <- generated_model$regression_model
    regression_model <- add_assessments(
      regression_model,
      evaluation_methods,
      FALSE
    )
    anova_res <- nlme::anova.lme(regression_model)
    regression_model$anova <- anova_res
    regression_model$aic <- stats::AIC(regression_model)
    regression_model$bic <- stats::BIC(regression_model)
    regression_model$aicc <- MuMIn::AICc(regression_model)

    models$nlme <- list(final_model = regression_model)
    if (trace) models$nlme$history <- NA
    return(models)
  }

  if (direction == "both" || direction == "forward") {
    models$forward <- forward_selection(
      df,
      model_type,
      term,
      evaluation_methods,
      categorical_vars,
      ...,
      model_family = model_family,
      trace = trace,
      omit_na = omit_na
    )
  }
  if (direction == "both" || direction == "backward") {
    models$backward <- backward_simplification(
      df,
      model_type,
      term,
      evaluation_methods,
      categorical_vars,
      ...,
      model_family = model_family,
      simplify_model = backward_simplify_model,
      trace = trace,
      omit_na = omit_na
    )
  }

  models
}

#' Remove NAs
#'
#' Small helper method to remove NAs from dataframe given model formula
#' @param df
#'  Dataframe with response and predictors as columns
#' @param term
#'  The formula to be used with the model. Can be either quote() or formula()
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam.
#' @return
#'  List containing the final regression model,
#'  the significant and marginally significant model variables
#'  and the selection history if trace is TRUE
remove_nas <- function(df, term, model_type) {
  columns_in_formula <- as.character(
    find_call(term, return = "atomic", df_cols = colnames(df))
  )
  response <- stats::terms(stats::as.formula(term))[[2]]
  cols_w_y <- append(as.character(response), columns_in_formula)

  if (model_type == "nlme") cols_w_y <- cols_w_y[cols_w_y %in% names(df)]

  df <- df[, cols_w_y]
  df <- stats::na.omit(df)
}

#' Backward regression model simplification
#'
#' Simplify the regression model by eliminating model variables one by one
#' @param df
#'  Dataframe with response and predictors as columns
#' @param model_type
#'  The model to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam. Default: glm.
#' @param term
#'  The formula to be used with the model.
#'  Can be either quote() or formula()
#' @param evaluation_methods
#'  Methods to be used for model evaluation.
#'  Options: anova, aic, aicc, bic, default: anova.
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param ...
#'  Parameters to be directly used with model call
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param simplify_model
#'  If FALSE, the model and information on significant and
#'  marginally significant variables are returned without any improvements.
#' @param trace
#'  Store and return model selection history (default: FALSE)
#' @param omit_na
#'  Either 'overall' or 'stepwise'.
#'  If 'overall', NAs are removed before modeling.
#'  If 'stepwise', NAs are removed per step based on
#'  the variables in the current formula
#' @return
#'  List containing the final regression model,
#'  the significant and marginally significant model variables
#'  and, if trace is TRUE, the selection history
backward_simplification <- function(
  df,
  model_type,
  term,
  evaluation_methods,
  categorical_vars,
  ...,
  model_family = "gaussian",
  simplify_model = TRUE,
  trace = FALSE,
  omit_na = "overall"
) {
  if (omit_na == "overall") df <- remove_nas(df, term, model_type)

  errors <- c()
  generated_model <- generate_regression_model(
    df,
    model_type,
    term,
    model_family,
    ...
  )
  model_args <- generated_model$model_args
  regression_model <- generated_model$regression_model

  if (trace) {
    history <- list()
    history[[
      stringr::str_interp(
        "${length(history)+1}_starting_model"
      )
    ]] <- regression_model
  }
  out <- list("final_model" = NA)
  if (model_type == "gam") {
    model_sum <- summary(regression_model)$p.table
  } else {
    model_sum <- stats::coef(summary(regression_model))
  }
  model_summary <- expand_model_summary(
    model_summary = model_sum,
    term = stats::formula(regression_model),
    categorical_vars = categorical_vars,
    df = df
  )

  if ((model_summary[1, "smallest_p_value"] < 0.1) || (!simplify_model)) {
    # All p-values are below p<0.1 --> accept the original model
    out$significant_variables <- model_summary[
      model_summary$p_value < 0.05,
      "predictor"
    ]
    out$marginally_significant_variables <- model_summary[
      (model_summary$p_value < 0.1) & (model_summary$p_value >= 0.05),
      "predictor"
    ]
    out$final_model <- regression_model
  } else {
    regression_model <- add_assessments(
      regression_model,
      evaluation_methods,
      model_type == "lmer"
    )
    simplify <- TRUE
    while (simplify) {
      res <- mo_step(
        df,
        regression_model,
        unique(
          model_summary[
            model_summary["smallest_p_value"] > 0.1,
            "predictor_term"
          ]
        ),
        "backward",
        evaluation_methods,
        categorical_vars,
        model_type,
        model_args = model_args
      )
      simplify <- res$simplify
      errors <- append(errors, res$errors)
      if (simplify) {
        regression_model <- res$regression_model
        if (omit_na == "stepwise") {
          current_term <- stats::formula(regression_model)
          regression_model <- generate_regression_model(
            df,
            model_type,
            current_term,
            model_family,
            ...
          )
          regression_model <- add_assessments(
            regression_model,
            evaluation_methods,
            model_type == "lmer"
          )
          if (model_type == "gam") {
            model_sum <- summary(regression_model)$p.table
          } else {
            model_sum <- stats::coef(summary(regression_model))
          }
          model_summary <- expand_model_summary(
            model_summary = model_sum,
            term = stats::formula(regression_model),
            categorical_vars = categorical_vars,
            df = df
          )
        } else {
          model_summary <- res$model_summary
        }
        if (trace) {
          history[[stringr::str_interp(
            "${length(history)+1}_${res$var}"
          )]] <- regression_model
        }
      }
    }

    var_cnt <- length(find_call(term, return = "fixed"))
    if ((var_cnt > 2) && (stats::formula(regression_model) == term)) {
      warning(
        sprintf(
          "We initialised the model successfully,
          but weren't able to simplify it, because %s.
          Please check overfitting and multicollinearity among
          predictors to exclude potential modeling issues.",
          res$reason
        )
      )
    }

    out$final_model <- regression_model
    out$significant_variables <- model_summary[
      model_summary$p_value < 0.05,
      "predictor"
    ]
    out$marginally_significant_variables <- model_summary[
      (model_summary$p_value < 0.1) & (model_summary$p_value >= 0.05),
      "predictor"
    ]
    if (trace) out$history <- history
  }

  out$errors <- errors
  out
}


#' Split and search model formula
#'
#' Extract part of model formula
#' @param term
#'  The formula to be used with the model. Can be either quote() or formula()
#' @param return
#'  Can be 'all', 'atomic', 'interactions', 'main_effects',
#'  'negative', or 'match'. Default: 'all'
#' @param pred_full
#'  Only needed for 'match' and 'negative'.
#'  The term to search for within the formula
#' @param df_cols
#'  Only needed if return == 'atomic'
#' @return
#'  The extracted part(s) of the model formula matching a given pattern.
#'  Can be either a call or a vector
find_call <- function(term, pred_full = NA, return = "all", df_cols = NA) {
  expr_terms <- stats::terms(stats::as.formula(term))
  term_labels <- colnames(attr(expr_terms, "factors"))

  if (return == "all") { # return as is
    out <- term_labels
  } else if (return == "fixed") {
    randoms <- reformulas::findbars(stats::as.formula(term))
    vars <- term_labels
    response <- formula.tools::lhs(stats::as.formula(term))
    if (length(response) > 1) {
      response <- deparse(formula.tools::lhs(response))
      random_idx <- sapply(vars, function(y) grepl(y, randoms, fixed = TRUE))
      out <- vars[(vars != response) & !random_idx]
    } else {
      response <- deparse(response)
      out <- vars[vars != response]
    }
  } else if (return == "random") {
    return(as.character(reformulas::findbars(stats::as.formula(term))))
  } else if (return == "atomic") {
    # return column names of fixed effects, i.e., remove transforms
    randoms <- reformulas::findbars(stats::as.formula(term))
    vars <- all.vars(stats::as.formula(term))
    response <- formula.tools::lhs(stats::as.formula(term))
    if (length(response) > 1) {
      response <- deparse(formula.tools::lhs(response))
    } else {
      response <- deparse(response)
    }

    out <- vars[vars != response]
  } else if (return == "interactions") {
    # subset for interactions
    out <- term_labels[stringr::str_detect(term_labels, ":")]
  } else if (return == "main_effects") {
    # subset for main effects
    out <- term_labels[!stringr::str_detect(term_labels, ":")]
  } else if (return == "negative") {
    # get term labels that don't match search term
    out <- term_labels[!stringr::str_detect(
      term_labels,
      stringr::fixed(pred_full)
    )]
  } else if (return == "match") {
    # get term labels matching search term
    out <- term_labels[stringr::str_detect(
      term_labels,
      stringr::fixed(pred_full)
    )]
  }

  out_v <- list()
  for (x in out) out_v <- append(out_v, parse(text = x)[[1]])
  unique(out_v)
}

#' Determine regression model family
#'
#' Determine model family automatically based on response column in dataframe
#' @param df
#'  Dataframe with response and predictors as columns
#' @param response_frm
#'  Response formula
#' @return
#'  Model family to be used with regression model
determine_model_family <- function(df, response_frm) {
  response_col <- paste(response_frm)
  df[response_col] <- with(df, eval(parse(text = response_frm)))

  response <- df[!is.na(df[[response_col]]), ][[response_col]]
  is_num <- is.numeric(response)
  is_int <- is.integer(response)
  gaus_or_pois <- is_num && ((min(response) < 0) || (max(response) > 1))

  if (is_int || is_num) {
    if ((min(response) < 0) || (max(response) > 100)) {
      model_family <- c("gaussian")
    } else if (is_int || (gaus_or_pois)) {
      model_family <- c("gaussian", "poisson")
    } else if (is_num && (min(response) >= 0) && (max(response) <= 1)) {
      model_family <- c("gaussian", "poisson", "quasibinomial")
    } else {
      model_family <- c()
    }
  } else if (is.logical(response)) {
    model_family <- c("binomial")
  } else {
    model_family <- c()
  }

  model_family
}

#' A minor version of expanded model summary
#'
#' Expand model summary by adding column names and
#'  information on interaction groups.
#' @param model_summary
#'  Initial model summary.
#' @param term
#'  The formula to be used with the model.
#'  Can be either quote() or formula().
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns).
#' @param df
#'  Dataframe with response and predictors as columns.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#'  Default: 3.
#' @return
#'  Expanded summary of a regression model with added information
#'  on base variables, variable types and significance levels.
#' @examples
#' # setup
#' data("plants")
#'
#' # generate an example glm model
#' final_model <- glm(
#'   sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     habitat +
#'     ploidy +
#'     solar_radiation:isothermality,
#'   family = "quasibinomial",
#'   data = plants
#' )
#'
#' # generate an extended model summary with information on predictor columns
#' # and variable types
#' models_overview <- stats::coef(summary(final_model))
#' expand_model_summary(
#'   models_overview,
#'   quote(sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     I(isothermality^2) +
#'     habitat +
#'     ploidy +
#'     solar_radiation:annual_mean_temperature +
#'     solar_radiation:isothermality +
#'     annual_mean_temperature:isothermality),
#'   c("habitat", "ploidy"),
#'   plants,
#'   3
#' )
#' @export
expand_model_summary <- function(
  model_summary,
  term,
  categorical_vars,
  df,
  round_p = 3
) {
  df_cols <- colnames(df)
  response_frm <- term[[2]]

  # define function to find main effects
  find_main_effect <- function(this_coeff) {
    if (grepl("^I\\(.+\\^[0-9]+\\)", this_coeff)) {
      this_coeff_small <- stringr::str_extract(
        this_coeff, "(?<=I\\().+(?=\\^[0-9])"
      )

      if (is.call(str2lang(this_coeff_small))) {
        # transformed var -> no main effect
        return(this_coeff_small)
      } else {
        for (el in df_cols) {
          ptrn <- stringr::str_interp("(?<!:)${el}(?![a-z0-9_\\-#:])")
          if (stringr::str_detect(
            this_coeff,
            stringr::regex(ptrn, ignore_case = TRUE)
          )) {
            return(el)
          }
        }
      }
    }

    this_coeff
  }

  get_pred_type <- function(coef1, coef2) {
    is_numeric <- c()

    for (coef in c(coef1, coef2)) {
      if (!is.na(coef)) {
        if (coef %in% df_cols) {
          is_numeric <- c(is_numeric, is.numeric(df[[coef]]))
        } else {
          for (el in df_cols) {
            ptrn <- stringr::str_interp("(?<!:)${el}(?![a-z0-9_\\-#:])")
            if (stringr::str_detect(
              coef,
              stringr::regex(ptrn, ignore_case = TRUE)
            )) {
              is_numeric <- c(is_numeric, is.numeric(df[[el]]))
            }
          }
        }
      }
    }

    if (length(unique(is_numeric)) == 1) {
      if (is_numeric[[1]] == TRUE) {
        "numeric"
      } else {
        "categorical"
      }
    } else {
      "mixed"
    }
  }

  get_min_p <- function(coef, coef_type) {
    if (coef_type == "numeric") {
      min_p <- min(
        expanded_model_summary[
          (expanded_model_summary$predictor == coef) |
            (((!is.na(expanded_model_summary$main_effect1)) &
              (expanded_model_summary$main_effect1 == coef)) |
              ((!is.na(expanded_model_summary$main_effect2)) &
                (expanded_model_summary$main_effect2 == coef))),
          "p_value"
        ]
      )
    } else if (coef_type == "categorical") {
      min_p <- min(
        expanded_model_summary[
          (expanded_model_summary$predictor == coef) |
            (expanded_model_summary$main_effect1 == coef),
          "p_value"
        ]
      )
    } else {
      min_p <- 0.0
    }

    min_p
  }

  elim_cat_trait <- function(predictor, main_effect1, main_effect2) {
    pred_wo_trait <- stringr::str_replace(
      predictor,
      sprintf(
        "(?<=\\:|^)[^\\:]*%s[^\\:]*(?=\\:|$)",
        main_effect2
      ),
      main_effect2
    )

    pred_wo_trait <- stringr::str_replace(
      pred_wo_trait,
      sprintf(
        "(?<=\\:|^)[^\\:]*%s[^\\:]*(?=\\:|$)",
        main_effect2
      ),
      main_effect2
    )
    pred_wo_trait
  }

  # define p-value column
  model_summary <- as.data.frame(model_summary)
  p_col <- intersect(
    c("Pr(>|t|)", "Pr(>|z|)", "Corrected_p"),
    colnames(model_summary)
  )
  if (length(p_col) != 1) {
    stop(stringr::str_interp("Not able to calculate p-values."))
  }

  # connect categorical trait values to traits
  cat_traits <- c()
  for (cat_var in categorical_vars) {
    for (cat_trait in paste(cat_var, unique(df[, cat_var]), sep = "")) {
      cat_traits[[cat_trait]] <- cat_var
    }
  }

  expanded_model_summary <- model_summary %>%
    tibble::rownames_to_column(var = "predictor") %>%
    dplyr::filter(.data$predictor != "(Intercept)") %>%
    dplyr::mutate(
      response = paste(response_frm),
      .before = "predictor"
    ) %>%
    dplyr::rename(p_value = tidyr::all_of(p_col)) %>%
    dplyr::mutate(
      is_interaction = dplyr::case_when(
        stringr::str_detect(.data$predictor, ":") ~ TRUE,
        TRUE ~ FALSE
      ),
      main_effect1 = {
        # remove transforms
        pred1 <- dplyr::case_when(
          .data$is_interaction ~ purrr::map_chr(
            stringr::str_extract(predictor, "^[^:]+"),
            find_main_effect
          ),
          TRUE ~ purrr::map_chr(
            predictor,
            find_main_effect
          )
        )

        # reduce categories to column names
        if (length(categorical_vars) == 0) {
          as.character(pred1)
        } else {
          dplyr::case_when(
            !!!setup_categorical_check(
              df, categorical_vars, quote(pred1),
              quote(pred1)
            )
          )
        }
      },
      main_effect2 = {
        pred2 <- dplyr::case_when(
          is_interaction ~ purrr::map_chr(
            stringr::str_extract(predictor, "[^:]+$"),
            find_main_effect
          )
        )

        # reduce categories to column names
        if (length(categorical_vars) == 0) {
          as.character(pred2)
        } else {
          dplyr::case_when(
            !!!setup_categorical_check(
              df, categorical_vars, quote(pred2),
              quote(pred2)
            )
          )
        }
      },
      pred_type = purrr::map2_chr(
        .data$main_effect1,
        .data$main_effect2,
        get_pred_type
      ),
      significance = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::mutate_if(is.numeric, round, digits = round_p) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      predictor_term = {
        if (.data$pred_type == "numeric") {
          .data$predictor
        } else if (.data$pred_type == "categorical") {
          .data$main_effect1
        } else {
          elim_cat_trait(
            .data$predictor,
            .data$main_effect1,
            .data$main_effect2
          )
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  expanded_model_summary <- expanded_model_summary %>%
    dplyr::mutate(
      smallest_p_value = dplyr::case_when(
        .data$pred_type == "categorical" ~ purrr::map2_dbl(
          .data$main_effect1,
          .data$pred_type,
          get_min_p
        ),
        .data$pred_type == "numeric" & .data$predictor == .data$main_effect1 ~
          purrr::map2_dbl(
            .data$predictor,
            .data$pred_type,
            get_min_p
          ),
        TRUE ~ .data$p_value
      )
    )


  non_sign_interactions <- expanded_model_summary[
    ((expanded_model_summary$p_value > 0.1) |
      (is.na(expanded_model_summary$p_value))) &
      (expanded_model_summary$is_interaction),
  ]

  if (nrow(non_sign_interactions) > 0) {
    expanded_model_summary <- expanded_model_summary %>%
      dplyr::arrange(
        dplyr::desc(.data$is_interaction),
        dplyr::desc(.data$smallest_p_value)
      )
  } else {
    expanded_model_summary <- expanded_model_summary %>%
      dplyr::arrange(dplyr::desc(.data$smallest_p_value))
  }

  expanded_model_summary
}


#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables
#'  that do not significantly predict response variable.
#' @param df
#'  Dataframe with response and predictors as columns.
#' @param term
#'  The formula to be used with the model. Can be either quote() or formula().
#' @param ...
#'  Arguments given directly to model call
#' @param autocorrelation_cols
#'  Sorted list of columns in dataframe to be considered when
#'  removing autocorrelations. Should be sorted by priority.
#'  Last element gets eliminated first.
#' @param automatic_removal
#'  Whether to automatically remove autocorrelations. Default: FALSE.
#' @param autocorrelation_threshold
#'  Threshold at which two variables are considered autocorrelated.
#'  Default: 0.7.
#' @param correlation_method
#'  The method used for correlation calculation.
#'  Default: "pearson".
#' @param cor_use
#'  Parameter 'use' for 'cor' method.
#'  Describes handling of missing values.
#'  Default: "complete.obs".
#' @param model_type
#'  Model type to be used.
#'  Options: (g)lm, (g)lmer, nlmer, and gam. Default: glm.
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param evaluation_methods
#'  Methods to be used for model evaluation.
#'  Options: anova, aic, aicc, bic. Default: anova.
#' @param simplification_direction
#'  Mode of stepwise model improvement.
#'  Either "forward" (i.e., forward selection),
#'  or "backward" (i.e., backward simplification),
#'  or "both". Default: "both".
#' @param backward_simplify_model
#'  If FALSE, the model and information on significant and
#'  marginally significant variables are returned without
#'  backward simplification. Default: TRUE.
#' @param omit_na
#'  Either "overall" or "stepwise".
#'  If "overall", NAs are removed before modeling.
#'  If "stepwise", NAs are removed per step based on
#'  the variables in the current formula. Default: "overall".
#' @param scale_predictor
#'  Whether to apply scaling to predictor variables. Default: FALSE.
#' @param plot_quality_assessment
#'  Module to use for plots for model quality assessment.
#'  Options: "performance", "baseR". Default: "baseR".
#' @param plot_relationships
#'  Whether to plot regression, effect size, and estimates.
#'  Default: FALSE.
#' @param jitter_plots
#'  Whether geom_point plots should use to jitter.
#'  Default: FALSE.
#' @param plot_type
#'  Either "boxplot" or "violin".
#'  Used to plot regression plots for categorical variables.
#'  Default: "boxplot".
#' @param stat_test
#'  Either "t.test" or "wilcox".
#'  Used to calculate statistics for regression plots of categorical variables.
#'  Default: "wilcox".
#' @param use_psi
#'  Whether to apply post-selection inference to optimized model.
#'  Only applicable to lm and glm models.
#'  Note that running post-selection inference may take some time.
#'  Default: FALSE.
#' @param psi_k
#'  The multiple of the number of degrees of freedom used as
#'  penalty in the post-selection inference model selection.
#'  The default k = 2 corresponds to the AIC. Default: 2.
#' @param psi_boot_repl
#'  A number or list of bootstrap replicates.
#'  The default is no bootstrapping. Default: 100.
#' @param psi_p_threshold
#'  P-value threshold to be used for significance filtering.
#'  Default: 0.05.
#' @param psi_label_size
#'  Size of test labels within post-selection inference plot.
#'  Default: 2.5.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#'  Default: 5.
#' @param trace
#'  Store and return model selection history. Default: FALSE.
#' @return
#'  List with a) information on autocorrelated variables and b)
#'  final simplified/expanded models with further information
#'  (see function [LazyModeler::simplify_model()] for further details)
#' @examples
#' # setup
#' data("plants")
#'
#' # generate and optimize a glm model using the provided
#' # dataframe and term. Check for correlations between the values
#' # of the provided dataframe columns and remove autocorrelated variables.
#' # Apply backward simplification to the model and plot the final model.
#' optimized_model_result <- optimize_model(
#'   plants,
#'   quote(sexual_seed_prop ~ altitude +
#'     latitude_gps_n +
#'     longitude_gps_e +
#'     (solar_radiation +
#'       annual_mean_temperature +
#'       isothermality)^2 +
#'     I(isothermality^2) +
#'     habitat +
#'     ploidy),
#'   autocorrelation_cols = c(
#'     "solar_radiation",
#'     "annual_mean_temperature",
#'     "isothermality",
#'     "altitude",
#'     "latitude_gps_n",
#'     "longitude_gps_e"
#'   ),
#'   automatic_removal = TRUE,
#'   autocorrelation_threshold = 0.8,
#'   correlation_method = "spearman",
#'   model_type = "glm",
#'   model_family = "quasibinomial",
#'   evaluation_methods = c("anova"),
#'   simplification_direction = "backward",
#'   backward_simplify_model = TRUE,
#'   omit_na = "overall",
#'   scale_predictor = TRUE,
#'   plot_quality_assessment = "performance",
#'   round_p = 3,
#'   cor_use = "complete.obs",
#'   plot_relationships = TRUE,
#'   jitter_plots = TRUE,
#'   plot_type = "boxplot",
#'   stat_test = "wilcox",
#'   trace = TRUE
#' )
#' @export
optimize_model <- function(
  df,
  term,
  ...,
  autocorrelation_cols = NA,
  automatic_removal = FALSE,
  autocorrelation_threshold = 0.7,
  correlation_method = "pearson",
  cor_use = "complete.obs",
  model_type = "glm",
  model_family = "gaussian",
  evaluation_methods = c("anova"),
  simplification_direction = "both",
  backward_simplify_model = TRUE,
  omit_na = "overall",
  scale_predictor = FALSE,
  plot_quality_assessment = "baseR",
  plot_relationships = FALSE,
  jitter_plots = FALSE,
  plot_type = "boxplot",
  stat_test = "wilcox",
  use_psi = FALSE,
  psi_k = 2,
  psi_boot_repl = 100,
  psi_p_threshold = 0.05,
  psi_label_size = 2.5,
  round_p = 5,
  trace = FALSE
) {
  errors <- c()
  # check for autocorrelations
  if (!all(is.na(autocorrelation_cols))) {
    are_cols_numeric <- unlist(
      lapply(autocorrelation_cols, function(pr) is.numeric(df[, pr]))
    )
    if (!all(are_cols_numeric)) {
      warning(
        "You asked us to check non-numeric columns for autocorrelations.
          We're gonna skip them for you, but please don't do it again. ;)"
      )
      autocorrelation_cols <- autocorrelation_cols[are_cols_numeric]
    }

    autocorrelations_and_preds <- remove_autocorrelations(
      df,
      autocorrelation_cols,
      use = cor_use,
      automatic_removal = automatic_removal,
      autocorrelation_threshold = autocorrelation_threshold,
      correlation_method = correlation_method
    )
    autocorrelations <- autocorrelations_and_preds$autocorrelations
  } else {
    autocorrelations <- NA
  }

  # some basic checks
  plot_perf <- plot_quality_assessment == "performance"
  if ((plot_perf) && (!rlang::is_installed("see"))) {
    warning(
      "Install 'see' for performance-based diagnostic plots.
    Falling back to base R plots."
    )
    plot_quality_assessment <- "baseR"
  }
  if (model_type == "nlme") {
    warning(
      "We're not using the lme4 package so we'd rather
      refer to this model type as nlme. ;)"
    )
    model_type <- "nlme"
  }
  if ((model_type %in% c("lm", "lmer")) && (model_family != "gaussian")) {
    warning(
      "When opting for lm and lmer models,
      you should always specify gaussian as the model family.
      But don't ya worry, I'll change the family for you. ;)"
    )
    model_family <- "gaussian"
  }
  if (!("anova" %in% evaluation_methods) && (grepl("quasi", model_family))) {
    warning(
      "Anova is the only method (that we support) that works with
      quasibinomial distributions. Imma change it for ya real quick,
      but please remember this for your future endeavors. ;)"
    )
    evaluation_methods <- c("anova")
  }

  # check that interactions + main effects are in formula (i.e., x:y + x + y)
  main_effects <- as.character(find_call(term, return = "main_effects"))
  interactions <- find_call(term, return = "interactions")
  need_to_add <- FALSE

  extract_main_effects <- function(interaction) {
    interaction_w_tr <- interaction[[1]] != ":"

    if (!interaction_w_tr) {
      stringr::str_split(deparse(interaction), ":")[[1]]
    } else {
      extract_main_effects(interaction[[2]])
    }
  }

  term <- stats::as.formula(term)
  for (interaction in interactions) {
    effects <- extract_main_effects(interaction)

    for (effect in effects) {
      new_term <- paste("~ . +", effect)
      if (!(effect %in% main_effects)) {
        term <- stats::update(term, new_term)
        need_to_add <- TRUE
      }
    }
  }

  if (need_to_add) {
    warning(
      "Please remember to add all main effects as separate items
      in your formula when including interactions.
      We've added all main effects for you this time.
      Might not next time, though, so better include 'em yourself. ;)"
    )
  }

  response_frm <- term[[2]]
  predictor_frm <- term[[3]]

  # remove autocorrelations from formula and gather categorical variables
  split_frm <- find_call(term, return = "all")

  categorical_vars <- vector()
  numerical_vars <- vector()
  remaining_pred_call <- vector()
  for (frm_part in split_frm) {
    frm_str <- deparse(frm_part)

    if ((frm_str %in% colnames(df)) && (is.factor(df[, frm_str]))) {
      categorical_vars <- append(categorical_vars, frm_str)
    } else if ((frm_str %in% colnames(df)) && (is.numeric(df[, frm_str]))) {
      numerical_vars <- append(numerical_vars, frm_str)
    }

    if (!all(is.na(autocorrelation_cols)) && (nrow(autocorrelations) > 0)) {
      ptrn <- sprintf(
        "(%s)(?![a-z0-9_\\-#])",
        paste(autocorrelations_and_preds$removed_predictors, collapse = "|")
      )
      if (!stringr::str_detect(
        frm_str,
        stringr::regex(ptrn, ignore_case = TRUE)
      )) {
        remaining_pred_call <- append(remaining_pred_call, frm_part)
      }
    }
  }

  if (!all(is.na(autocorrelation_cols)) && (nrow(autocorrelations) > 0)) {
    predictor_frm <- remaining_pred_call[[1]]
    for (i in 2:length(remaining_pred_call)) {
      predictor_frm <- call("+", predictor_frm, remaining_pred_call[[i]])
    }
  }

  if (scale_predictor) {
    for (numerical_var in numerical_vars) {
      df[
        stringr::str_interp("${numerical_var}_aB3cD5eF6G")
      ] <- df[, numerical_var]
      df[numerical_var] <- as.vector(scale(df[, numerical_var]))
    }
  }

  var_cnt <- length(find_call(term, return = "atomic"))
  if (nrow(df) / var_cnt <= 4) {
    warning(
      stringr::str_interp(
        "There are probably too many variables (${var_cnt})
        in comparison to datapoints (${nrow(df)}).
        Please add data, or remove x variables to reduce overfitting,
        multicollinearity, convergence issues,
        statistical significance and interpretability.",
      )
    )
  }

  term <- substitute(
    response ~ predictor,
    list(response = response_frm, predictor = predictor_frm)
  )

  if (model_type != "nlme") {
    possible_model_family <- determine_model_family(df, response_frm)
    if (model_family == "automatic") {
      if (length(possible_model_family) > 1) {
        print(
          sprintf(
            "Your response variable allows for %s distribution.
            Which one would you prefer?",
            paste(possible_model_family, collapse = " and ")
          )
        )
        model_family <- readline()
      } else if (length(possible_model_family) == 0) {
        stop("No appropriate family found. Please check data types of columns.")
      } else {
        model_family <- possible_model_family[1]
        print(stringr::str_interp(
          "Continuing with ${model_family} distribution."
        ))
      }
    } else if (!(model_family %in% possible_model_family)) {
      warning(
        sprintf(
          "Chosen distribution '%s' does not match response values.
          Would recommend %s. Please check.",
          model_family,
          paste(possible_model_family, collapse = ", or ")
        )
      )
    }
  }

  simplified_models <- simplify_model(
    df,
    model_type,
    term,
    evaluation_methods,
    ...,
    model_family = model_family,
    categorical_vars = categorical_vars,
    direction = simplification_direction,
    backward_simplify_model = backward_simplify_model,
    trace = trace,
    omit_na = omit_na
  )

  for (direction in names(simplified_models)) {
    simplified_model_info <- simplified_models[[direction]]
    for (error in unique(simplified_model_info$errors)) warning(error)
    final_model <- simplified_model_info$final_model

    if ((model_type %in% c("glm", "lm")) && use_psi) {
      final_model_before_psi <- final_model
      models_overview_for_psi <- expand_model_summary(
        stats::coef(summary(final_model)),
        term,
        categorical_vars,
        df,
        round_p
      )

      psi_result <- apply_psi(
        df,
        final_model,
        model_type,
        model_family,
        term,
        categorical_vars,
        models_overview_for_psi,
        boot_repl = psi_boot_repl,
        k = psi_k,
        round_p = round_p
      )
      psi_plot <- plot_psi(
        psi_result$result,
        p_threshold = psi_p_threshold,
        label_size = psi_label_size
      )

      psi_result$plot <- psi_plot
      final_model <- psi_result$psi_model
    }

    if (scale_predictor) {
      for (numerical_var in numerical_vars) {
        df[numerical_var] <- df[
          ,
          stringr::str_interp("${numerical_var}_aB3cD5eF6G")
        ]
        df[numerical_var] <- df[numerical_var] %>%
          dplyr::select(!c(!!numerical_var))
      }
    }

    if ((!is.list(final_model)) && (!typeof(final_model) == "S4")) {
      simplified_models[[direction]] <- simplified_model_info
      next
    }

    all_sign_vars <- c(
      simplified_model_info$significant_variables,
      simplified_model_info$marginally_significant_variables
    )
    all_sign_vars <- all_sign_vars[!is.na(all_sign_vars)]

    if (length(all_sign_vars) != 0) {
      # CHECK MODEL
      model_plots <- list()
      if (plot_quality_assessment == "performance") {
        if (model_type == "gam") {
          plot(performance::check_model(final_model, residual_type = "normal"))
        } else if (model_family == "binomial") {
          withCallingHandlers(
            plot(
              performance::check_model(final_model, type = "discrete_both")
            ),
            warning = function(w) {
              if (grepl("stat_density", w$message)) {
                warning(
                  sprintf(
                    "%s. Please choose another model family.",
                    w$message
                  )
                )
                tryInvokeRestart("muffleWarning")
              } else {
                warning(w$message)
              }
            }
          )
        } else {
          withCallingHandlers(
            plot(performance::check_model(final_model)),
            message = function(w) {
              if (grepl("unknown labels", w$message)) {
                tryInvokeRestart("muffleMessage")
              }
            }
          )
        }
        checked_model <- grDevices::recordPlot()
        grDevices::dev.off()
      } else if (plot_quality_assessment == "baseR") {
        graphics::par(mfrow = c(2, 2), mar = rep(2, 4))
        plot(final_model)
        checked_model <- grDevices::recordPlot()
        grDevices::dev.off()
      } else {
        checked_model <- NA
      }
      model_plots$model_check <- checked_model

      if (model_type == "gam") {
        models_overview <- summary(final_model)$p.table
        models_overview_s <- summary(final_model)$s.table
        models_overview_s <- expand_model_summary(
          models_overview_s,
          term,
          categorical_vars,
          df,
          round_p
        )
      } else {
        models_overview <- stats::coef(summary(final_model))
      }

      models_overview <- expand_model_summary(
        models_overview,
        term,
        categorical_vars,
        df,
        round_p
      )
      models_overview <- models_overview %>%
        dplyr::mutate(effect_direction = dplyr::case_when(
          .data$Estimate < 0 ~ "negative",
          TRUE ~ "positive"
        ))

      if (plot_relationships && (nrow(models_overview) > 0)) {
        plot_out <- plot_model_features(
          final_model,
          models_overview,
          model_type,
          model_family,
          plot_type,
          jitter_plots,
          stat_test,
          round_p,
          !backward_simplify_model
        )
        model_plots <- append(model_plots, plot_out)
      }
      if (plot_relationships && (model_type == "gam")) {
        plot_out <- plot_model_features(
          final_model,
          models_overview_s,
          model_type,
          model_family,
          plot_type,
          jitter_plots,
          stat_test,
          round_p,
          !backward_simplify_model
        )
        model_plots <- append(model_plots, plot_out)
      }
    } else {
      simplified_models[[direction]] <- simplified_model_info
      next
    }

    models_overview <- models_overview %>%
      dplyr::rename(
        Response = "response",
        Predictor = "predictor",
        `Variable Type` = "pred_type",
        Significance = "significance",
        `Effect Direction` = "effect_direction"
      ) %>%
      dplyr::select(-c(
        "main_effect1",
        "main_effect2"
      ))

    if (model_type == "gam") {
      models_overview <- list(
        p.table = models_overview,
        s.table = models_overview_s
      )
    }

    simplified_models[[direction]] <- list(
      overview = models_overview,
      final_model = simplified_model_info$final_model,
      plots = model_plots
    )

    if (trace) {
      simplified_models[[direction]]$history <- simplified_model_info$history
    }

    if ((model_type %in% c("glm", "lm")) && use_psi) {
      simplified_models[[direction]]$psi <- psi_result
      simplified_models[[direction]]$model_before_psi <- final_model_before_psi
    }
  }

  list(
    autocorrelations = autocorrelations,
    models_with_info = simplified_models
  )
}

#' Plot model features
#'
#' Plot estimate, regression and effect size
#' @param regression_model
#'  The (simplified) regression model.
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()]).
#' @param model_type
#'  Type of regression model (e.g., 'glm').
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options.
#' @param plot_type
#'  Whether to plot as boxplot or violin plot
#'  for categorical variables. Default: "boxplot".
#' @param jitter_plots
#'  Whether to use jitter when generating geom_point plots.
#'  Default: FALSE.
#' @param test
#'  Test to be used for pairwise significance testing for
#'  categorical variables. Either wilcox or t.test. Default: wilcox.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#'  Default: 5.
#' @param remove_insignificant
#'  Used to exclude insignificant relationships from both estimate
#'  and effect size plots and omit curves from regression plots.
#'  Default: FALSE.
#' @return
#'  A list of plots (estimate, regression, and effect size)
#'  and statistics for categorical variables.
#' @examples
#' # setup
#' library(dplyr)
#' data("plants")
#' model_type <- "glm"
#' model_family <- "quasibinomial"
#'
#' # generate an example glm model for plotting
#' final_model <- glm(
#'   sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     habitat +
#'     ploidy +
#'     solar_radiation:isothermality,
#'   family = model_family,
#'   data = plants
#' )
#'
#' # generate extended model summary with information on predictor columns
#' # and variable types
#' models_overview <- stats::coef(summary(final_model))
#' models_overview <- LazyModeler::expand_model_summary(
#'   models_overview,
#'   quote(sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     I(isothermality^2) +
#'     habitat +
#'     ploidy +
#'     solar_radiation:annual_mean_temperature +
#'     solar_radiation:isothermality +
#'     annual_mean_temperature:isothermality),
#'   c("habitat", "ploidy"),
#'   plants,
#'   3
#' )
#' models_overview <- models_overview %>%
#'   mutate(effect_direction = case_when(
#'     Estimate < 0 ~ "negative",
#'     TRUE ~ "positive"
#'   ))
#'
#' # plot the model's features including plotting of estimates,
#' # numeric and categorical variables using boxplots where possible
#' plot_model_features(
#'   final_model,
#'   models_overview,
#'   model_type,
#'   model_family,
#'   plot_type = "boxplot",
#'   jitter_plots = TRUE,
#'   test = "wilcox",
#'   round_p = 3,
#'   remove_insignificant = FALSE
#' )
#' @export
plot_model_features <- function(
  regression_model,
  models_overview,
  model_type,
  model_family,
  plot_type = "boxplot",
  jitter_plots = FALSE,
  test = "wilcox",
  round_p = 5,
  remove_insignificant = FALSE
) {
  if (remove_insignificant) {
    models_overview_ns <- models_overview[
      models_overview$significance != "ns",
    ]
  } else {
    models_overview_ns <- models_overview
  }

  model_plots <- list(
    regression_plots = plot_regression(
      regression_model,
      models_overview,
      model_type,
      model_family,
      plot_type,
      jitter_plots,
      test,
      round_p,
      !remove_insignificant
    )
  )

  if ("Estimate" %in% colnames(models_overview)) {
    model_plots[["estimate_plot"]] <- plot_estimate(
      models_overview_ns,
      regression_model,
      model_type
    )
    model_plots[["effect_size_plot"]] <- plot_effect_size(models_overview_ns)
  }

  model_plots
}

#' Plot regression model estimate
#'
#' Plot estimate, regression and effect size
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()])
#' @param regression_model
#'  The (simplified) regression model
#' @param model_type
#'  Type of regression model (e.g., 'glm')
#' @return
#'  A plot of the estimates listed within the model summary
plot_estimate <- function(models_overview, regression_model, model_type) {
  if ((model_type == "lmer") || (model_type == "glmer")) {
    df <- stats::model.frame(regression_model)
  } else {
    df <- regression_model$model
  }

  new_rows <- lapply(
    unique(models_overview[
      models_overview["pred_type"] == "categorical",
      "main_effect1"
    ]),
    function(cat_var) {
      cat_var_states <- unique(df[[cat_var]])
      reference <- setdiff(
        cat_var_states,
        sub(cat_var, "", models_overview$predictor)
      )
      data.frame(
        predictor = paste0(cat_var, reference),
        main_effect1 = cat_var,
        Estimate = 0,
        `Std. Error` = 0,
        significance = "ref",
        stringsAsFactors = FALSE
      )
    }
  )
  models_overview <- dplyr::bind_rows(models_overview, do.call(rbind, new_rows))
  models_overview <- models_overview %>%
    dplyr::arrange(
      dplyr::desc(.data$main_effect1),
      dplyr::desc(.data$predictor)
    ) %>%
    dplyr::mutate(
      formatted_predictor = dplyr::case_when(
        (.data$significance != "ref") & (.data$pred_type == "categorical") ~
          paste0("italic(", .data$predictor, ")"),
        TRUE ~ .data$predictor
      ),
      formatted_predictor = factor(
        .data$formatted_predictor,
        levels = .data$formatted_predictor
      )
    )

  formatted_labels <- sapply(models_overview$formatted_predictor, function(x) {
    if (grepl("^italic", x)) {
      parse(text = as.character(x))
    } else {
      as.character(x)
    }
  })

  p <- ggplot2::ggplot(
    models_overview,
    ggplot2::aes(x = .data$Estimate, y = .data$formatted_predictor)
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 0)
    ) +
    ggplot2::labs(y = "Predictor") +
    ggplot2::scale_y_discrete(labels = formatted_labels) +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = .5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Estimate)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = .data$Estimate - .data$`Std. Error`,
        xmax = .data$Estimate + .data$`Std. Error`,
        color = .data$Estimate
      ),
      width = .0
    ) +
    colorspace::scale_color_continuous_diverging(rev = TRUE)

  p
}

#' Plot model regression (categorical)
#'
#' Plot categorical regression of model variables
#' @param cat_var
#'  Categorical variable to plot
#' @param df
#'  Dataframe object extracted from model
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()])
#' @param test
#'  Test to be used for pairwise significance testing
#'  for categorical variables. Either wilcox or t.test
#' @param response_str
#'  Response as string
#' @param response_col
#'  Response column name
#' @param plot_type
#'  Whether to plot as boxplot or violin plot for categorical variables
#' @param regression_plots
#'  A list storing regression plots
#' @param stat_results
#'  A list storing results of statistics
#' @return
#'  Regression plots of variables listed within
#'  the model summary + statistics
plot_regression_categorical <- function(
  cat_var,
  df,
  models_overview,
  test,
  response_str,
  response_col,
  plot_type,
  regression_plots,
  stat_results
) {
  models_overview_group <- models_overview %>%
    dplyr::filter(.data$main_effect1 == cat_var) %>%
    dplyr::mutate(factor_group = stringr::str_remove(
      .data$predictor,
      .data$main_effect1
    ))

  new_rows <- lapply(
    unique(
      models_overview_group[
        models_overview_group["pred_type"] == "categorical",
        "main_effect1"
      ]
    ), function(cat_var) {
      cat_var_states <- unique(df[[cat_var]])
      reference <- setdiff(
        cat_var_states, sub(cat_var, "", models_overview_group$predictor)
      )
      data.frame(
        response = models_overview_group[1, "response"],
        predictor = paste0(cat_var, reference),
        main_effect1 = cat_var,
        Estimate = 0,
        significance = "ref",
        stringsAsFactors = FALSE
      )
    }
  )
  models_overview_group <- dplyr::bind_rows(
    models_overview_group, do.call(rbind, new_rows)
  )

  df <- merge(
    df, models_overview_group,
    by.x = cat_var, by.y = "factor_group",
    all.x = TRUE, all.y = FALSE
  )

  cat_var <- parse(text = cat_var)[[1]]
  response_sym <- parse(text = response_str)[[1]]

  if (nrow(models_overview_group) > 1) { # multiple vars
    stat_out <- run_stats(df, response_col, cat_var, test = test)
    df <- stat_out$df
    stat_results[[cat_var]] <- stat_out$stat_results

    df <- df %>%
      dplyr::mutate(significance = ifelse(
        is.na(.data$significance), "ref", .data$significance
      ))
    p <- ggplot2::ggplot(
      data = df,
      ggplot2::aes(x = !!cat_var, y = !!response_sym)
    )
    if (plot_type == "boxplot") {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(color = !!cat_var, fill = !!cat_var)
      )
    } else {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(color = !!cat_var, fill = !!cat_var)
      )
    }

    p <- p +
      ggplot2::scale_color_viridis_d(option = "G") +
      ggplot2::scale_fill_viridis_d(option = "G") +
      ggplot2::geom_text(
        ggplot2::aes(x = !!cat_var, label = .data$letter),
        y = max(df[!is.na(df[response_str]), response_str]) * 1.15,
        check_overlap = TRUE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = !!cat_var, label = .data$significance),
        y = max(df[!is.na(df[response_str]), response_str] * 1.05),
        check_overlap = TRUE
      ) +
      ggplot2::ylim(
        min(df[!is.na(df[response_str]), response_str]),
        max(df[!is.na(df[response_str]), response_str]) * 1.15
      )
  } else { # single vars
    p <- ggplot2::ggplot(
      data = df,
      ggplot2::aes(x = !!cat_var, y = !!response_sym)
    )
    if (plot_type == "boxplot") {
      p <- p + ggplot2::geom_boxplot(color = "#21918c")
    } else {
      p <- p + ggplot2::geom_violin(color = "#21918c")
    }
  }

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  regression_plots[[cat_var]] <- p

  list(
    regression_plots = regression_plots,
    stat_results = stat_results
  )
}

#' Plot model regression (numeric & mixed)
#'
#' Plot numeric/mixed regression of model variables.
#' @param predictor_full
#'  Numeric/mixed variable to plot.
#' @param regression_model
#'  The (simplified) regression model.
#' @param df
#'  Dataframe object extracted from model.
#' @param term
#'  The formula to be used with the model.
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()]).
#' @param response_str
#'  Response as string.
#' @param regression_plots
#'  A list storing regression plots.
#' @param jitter
#'  Whether to use jitter when generating dotplots.
#' @param plot_curve
#'  Whether to plot regression plot curve using geom_smooth.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options.
#' @param model_type
#'  Type of regression model (e.g., 'glm').
#' @return
#'  Returns regression plots of variables listed within the model summary.
plot_regression_numeric <- function(
  predictor_full,
  regression_model,
  df,
  term,
  models_overview,
  response_str,
  regression_plots,
  jitter,
  plot_curve,
  round_p,
  model_family,
  model_type
) {
  if ("Estimate" %in% colnames(models_overview)) {
    estimate <- round(
      models_overview[
        models_overview["predictor"] == predictor_full,
        "Estimate"
      ],
      digits = round_p
    )
  } else {
    estimate <- "NA"
  }

  sign <- models_overview[
    models_overview["predictor"] == predictor_full,
    "significance"
  ]

  # check if columns are numeric
  predictor_data_type <- models_overview[
    models_overview$predictor == predictor_full,
    "pred_type"
  ][[1]]

  pred_main_effects <- c(models_overview[
    models_overview$predictor == predictor_full,
    c("main_effect1", "main_effect2")
  ])
  pred_main_effects <- pred_main_effects[!is.na(pred_main_effects)]

  transformed_vars <- transform_variables(
    regression_model,
    models_overview,
    df,
    term,
    predictor_full,
    pred_main_effects,
    predictor_data_type
  )

  df <- transformed_vars$df
  predictor_full <- transformed_vars$predictor_full
  preds_w_tr <- transformed_vars$preds_w_tr
  col_name <- transformed_vars$col_name
  
  response_sym <- parse(text = response_str)[[1]]

  if ((length(pred_main_effects) == 1) || (predictor_data_type == "numeric")) {
    # either one variable or only continuous variables
    if (length(pred_main_effects) == 1) {
      x <- predictor_full
    } else {
      x <- col_name
    }

    x <- parse(text = x)[[1]]

    if (stringr::str_detect(deparse(x), stringr::regex("I\\(.+\\^2\\)"))) {
      x_base <- stringr::str_extract(
        deparse(x),
        stringr::regex("(?<=I\\().+(?=\\^2\\))")
      )
      new_x <- sprintf(
        "%s_squared_%s",
        x_base,
        format(Sys.time(), format = "%Y%m%d%H%M")
      )
      df <- df %>%
        dplyr::rename(!!new_x := deparse(x))
      x <- parse(text = new_x)[[1]]
      df[[new_x]] <- as.numeric(df[[new_x]])
    } else if (stringr::str_detect(deparse(x), ":")) {
      new_x <- gsub(":", "_", deparse(x))
      df <- as.data.frame(stats::model.matrix(regression_model))
      df <- df %>%
        dplyr::rename(!!new_x := deparse(x))
      df[[response_str]] <- as.data.frame(
        stats::model.frame(regression_model)
        )[,response_str]
      x <- parse(text = new_x)[[1]]
    }

    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = !!x, y = !!response_sym)) +
      ggplot2::geom_point(
        color = "white",
        fill = "#21918c",
        pch = 21,
        position = jitter
      )

    if (plot_curve && sign != "ns") {
      p <- p + ggplot2::geom_smooth(
        method = model_type,
        method.args = list(family = model_family),
        formula = y ~ x,
        se = FALSE,
        color = "#21918c",
        fill = "#21918c"
      )
    }
  } else if (predictor_data_type == "mixed") { # >1 variable with one continuous
    x <- sapply(preds_w_tr, function(pr) is.numeric(df[, pr]))
    sign_vars_sorted <- c(preds_w_tr[x], preds_w_tr[!x])

    x <- parse(text = sign_vars_sorted[1])[[1]]
    color <- parse(text = sign_vars_sorted[2])[[1]]

    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = !!x, y = !!response_sym)) +
      ggplot2::geom_point(ggplot2::aes(color = !!color), position = jitter) +
      ggplot2::scale_color_viridis_d(option = "G") +
      ggplot2::scale_fill_viridis_d(option = "G")

    if (plot_curve && sign != "ns") {
      p <- p + ggplot2::geom_smooth(
        method = model_type,
        formula = y ~ x,
        method.args = list(family = model_family),
        se = TRUE,
        ggplot2::aes(color = !!color, fill = !!color)
      )
    }
  }

  x_max <- max(df[[deparse(x)]])
  x_min <- min(df[[deparse(x)]])

  est_len <- length(as.character(estimate))

  if (estimate > 0 || is.na(estimate)) {
    x <- x_min + ((x_max - x_min) * (.1 + est_len * .01))
  } else {
    x <- x_max - ((x_max - x_min) * (.1 + est_len * .01))
  }

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = "none") +
    ggplot2::annotate(
      "label",
      x = x,
      y = Inf,
      label = stringr::str_interp("estimate=${estimate}\n${sign}"),
      vjust = 1,
      hjust = .5
    )
  regression_plots[[predictor_full]] <- p

  regression_plots
}

#' Transform dataframe variable
#'
#' Transform variable name and values within given dataframe
#'  according to full predictor name and term.
#' @param regression_model
#'  The (simplified) regression model.
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()]).
#' @param df
#'  Dataframe object extracted from model.
#' @param term
#'  The complete model term.
#' @param predictor_full
#'  Coefficient as stated in term.
#' @param main_effects
#'  Main effects included in predictor_full.
#' @param pred_data_type
#'  Information on whether main_effects are numeric.
#'  Either: "numeric", "categorical", or "mixed".
#' @return
#'  Returns list with a) dataframe with transformed variables,
#'  b) adjusted predictor_full corresponding to transformed
#'  column name, c) a list with variables with transforms, and
#'  d) the name of a generated column in case of a coefficient interaction.
transform_variables <- function(regression_model,
                                models_overview,
                                df,
                                term,
                                predictor_full,
                                main_effects,
                                pred_data_type) {
  preds_w_tr <- vector()
  col_name <- NA_character_
  # transform df variables according to formula
  if ((pred_data_type != "numeric") &&
    (lengths(stringr::str_extract_all(predictor_full, ":")) > 1)) {
    warning(
      stringr::str_interp(
        "We don't allow for plotting of >2 interactions for
        categorical variables. Skipping plotting for ${predictor_full}."
      )
    )
    return()
  } else if (length(main_effects) == 2) { # interaction
    # retrieve categorical info if any
    cat_traits <- c()
    for (main_effect in main_effects) {
      pred_type <- models_overview[
        (models_overview$main_effect1 == main_effect) |
          (is.na(models_overview$main_effect2)),
        "pred_type"
      ][[1]]
      if (pred_type == "categorical") {
        for (cat_trait in paste(
          main_effect, unique(df[, main_effect]),
          sep = ""
        )) {
          cat_traits[[cat_trait]] <- main_effect
        }
      }
    }

    # === turn predictor into formula call +
    # === retrieve transform interaction info
    complete_call <- str2lang(predictor_full)
    interaction_transformed <- complete_call[[1]] != ":"

    # === transform individual variables
    for (pred_w_tr in main_effects) {
      if (pred_w_tr %in% names(cat_traits)) {
        preds_w_tr[[length(preds_w_tr) + 1]] <- cat_traits[[pred_w_tr]]
      } else {
        if (is.call(str2lang(pred_w_tr))) {
          new_col <- gsub("\\^|\\(|\\)", ".", pred_w_tr,
            perl = TRUE
          )
          df[new_col] <- with(df, eval(str2lang(pred_w_tr)))
        } else {
          new_col <- pred_w_tr
        }
        preds_w_tr[[length(preds_w_tr) + 1]] <- new_col
      }
    }

    # === compute interaction if continuous
    if (pred_data_type == "numeric") { # only continuous variables
      col_name <- if (interaction_transformed) {
        stringr::str_interp("${preds_w_tr[[1]]}_${preds_w_tr[[2]]}")
      } else {
        predictor_full
      }
      df[col_name] <- df[preds_w_tr[[1]]] * df[preds_w_tr[[2]]]
    }

    # === check if transform needs to be applied for interaction
    if ((pred_data_type == "numeric") && interaction_transformed) {
      # transform around interaction
      interaction_call <- deparse(complete_call[[2]])
      tr_inter_str <- gsub(
        interaction_call,
        col_name,
        deparse(complete_call, width.cutoff = 500)
      )
      df[predictor_full] <- with(df, eval(parse(text = tr_inter_str)))
    } else if (length(complete_call) != 3) {
      # 3 == no transform around interaction,
      # different value needs to be checked
      warning(
        stringr::str_interp(
          "Internal error. Check interaction call.
          No plot generated for ${predictor_full}"
        )
      )
    }
  } else if (stringr::str_detect(predictor_full, stringr::regex("s\\(.+\\)"))) {
    # no interaction + s() transform
    col_name <- stringr::str_extract(predictor_full, "(?<=\\().+(?=\\))")
    df[stringr::str_interp("${col_name}_smooth")] <- stats::predict(
      regression_model,
      type = "terms"
    )[, predictor_full]
    predictor_full <- stringr::str_interp("${col_name}_smooth")
  }

  list(
    df = df, predictor_full = predictor_full,
    preds_w_tr = preds_w_tr, col_name = col_name
  )
}

#' Plot model regression
#'
#' Plot regression of model variables.
#' @param regression_model
#'  The (simplified) regression model.
#' @param models_overview
#'  Custom model summary (see [LazyModeler::expand_model_summary()]).
#' @param model_type
#'  Type of regression model (e.g., 'glm').
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options.
#' @param plot_type
#'  Whether to plot as boxplot or violin plot for
#'  categorical variables. Default: boxplot.
#' @param jitter_plots
#'  Whether to use jitter when generating dotplots. Default: FALSE.
#' @param test
#'  Test to be used for pairwise significance testing for
#'  categorical variables. Either wilcox or t.test. Default: wilcox.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#' @param plot_curve
#'  Whether to plot regression plot curve using geom_smooth.
#' @return
#'  Regression plots of variables listed within the model summary.
plot_regression <- function(
  regression_model,
  models_overview,
  model_type,
  model_family,
  plot_type = "boxplot",
  jitter_plots = FALSE,
  test = "wilcox",
  round_p = 5,
  plot_curve = FALSE
) {
  if (jitter_plots) {
    jitter <- "jitter"
  } else {
    jitter <- "identity"
  }
  if ((model_type == "lmer") || (model_type == "glmer")) {
    df <- stats::model.frame(regression_model)
  } else {
    df <- regression_model$model
  }

  term <- stats::formula(regression_model)
  response_col <- parse(text = all.vars(term)[1])[[1]]

  regression_plots <- list()
  stat_results <- list()

  # transform response
  response_frm <- term[[2]]
  response_str <- paste(response_frm)
  if (is.call(response_frm)) {
    df[response_str] <- with(df, eval(parse(text = response_frm)))
  }

  # plot categorical vars (if no interaction with numeric)
  for (cat_var in unique(
    models_overview[
      models_overview$pred_type == "categorical",
      "main_effect1"
    ]
  )) {
    res <- plot_regression_categorical(
      cat_var,
      df,
      models_overview,
      test,
      response_str,
      response_col,
      plot_type,
      regression_plots,
      stat_results
    )
    regression_plots <- res$regression_plots
    stat_results <- res$stat_results
  }

  # for continuous or mixed variables
  relevant_predictors <- models_overview[
    models_overview$pred_type %in% c("numeric", "mixed"),
  ] %>%
    dplyr::filter(
      (.data$pred_type != "mixed") |
        (!duplicated(.data$main_effect1, .data$main_effect2))
    )
  for (predictor_full in relevant_predictors$predictor) {
    regression_plots <- plot_regression_numeric(
      predictor_full,
      regression_model,
      df,
      term,
      models_overview,
      response_str,
      regression_plots,
      jitter,
      plot_curve,
      round_p,
      model_family,
      model_type
    )
  }

  list(regression_plots = regression_plots, stat_results = stat_results)
}

#' Plot regression model effect sizes
#'
#' Plot effect sizes of model variables.
#' @param models_overview
#'  Custom model summary (see function [LazyModeler::expand_model_summary()]).
#' @return
#'  A plot of the effect sizes of the variables listed within the model summary.
plot_effect_size <- function(models_overview) {
  models_overview <- models_overview %>%
    dplyr::mutate(
      Estimate_abs = abs(.data$Estimate),
      Est_sum = sum(.data$Estimate_abs),
      `Effect size` = dplyr::case_when(
        .data$Estimate > 0 ~ (.data$Estimate_abs / .data$Est_sum) * 100,
        .data$Estimate < 0 ~ (.data$Estimate_abs / .data$Est_sum) * -100
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$`Effect size`))

  p <- ggplot2::ggplot(
    models_overview,
    ggplot2::aes(
      y = .data$`Effect size`,
      x = .data$predictor,
      fill = .data$`Effect size`
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    colorspace::scale_fill_continuous_diverging(rev = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::ylab("Effect size [%]") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}

#' Function for calculating statistics for given categorical variable
#'
#' Calculate pairwise statistics for given categorical variable and
#'  add letters indicating significant differences
#' @param df
#'  Dataframe with response and predictors as columns.
#' @param response
#'  Character string representing the response/y variable (column name)
#' @param predictor
#'  Character string representing the predictor/x variable (column name)
#' @param test
#'  Test to be used for pairwise significance testing for categorical variables.
#'  Either wilcox or t.test.
#' @return
#'  Output of the statistical test and the updated dataframe with letters
#'  indicating significant differences between traits of a categorical variable
run_stats <- function(df, response, predictor, test = "wilcox") {
  add_letters <- function(df, stat_results, test) {
    df_sorted <- df %>%
      dplyr::group_by(!!predictor) %>%
      dplyr::mutate(resp_mean = mean(!!response)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(!!predictor, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(!!response)) %>%
      dplyr::mutate(idx = dplyr::row_number())
    if (any(stat_results$p_value < 0.05)) {
      stat_results <- merge(
        stat_results,
        df_sorted[c(deparse(predictor), "idx")],
        by.x = "var1",
        by.y = deparse(predictor)
      )
      stat_results <- merge(
        stat_results,
        df_sorted[c(deparse(predictor), "idx")],
        by.x = "var2", by.y = deparse(predictor),
        suffixes = c("1", "2")
      )
      # for accessibility
      stat_results <- stat_results %>%
        dplyr::rowwise() %>%
        dplyr::mutate(comb = paste(
          sort(c(.data$idx1, .data$idx2)),
          collapse = "|"
        )) %>%
        dplyr::ungroup()

      letter <- "A"
      for (i in seq_len(nrow(df_sorted))) {
        ab_check <- paste(i - 1, i, sep = "|")

        if (i == 1) { # highest accuracy
          df_sorted[1, "letter"] <- letter
        } else if (i != nrow(df_sorted)) {
          # everything between highest and lowest accuracy
          # capital letter in comment is the model we're assigning
          # the letter to at the moment ;)
          aab_check <- stat_results[stat_results$comb == paste(
            i, i + 1,
            sep = "|"
          ), ]
          aabb_check <- stat_results[stat_results$comb == paste(
            i - 1, i + 1,
            sep = "|"
          ), ]

          if (stat_results[stat_results$comb == ab_check, ]$p_value < 0.05) {
            # i and i-1 are significantly different (a-B lettering)
            new_letter <- LETTERS[match(letter, LETTERS) + 1]
            df_sorted[i, "letter"] <- new_letter
            letter <- new_letter
          } else if (aab_check$p_value < 0.05) {
            # a-A-b
            df_sorted[i, "letter"] <- letter
          } else if (aabb_check$p_value < 0.05) {
            # a-AB-b
            new_letter <- LETTERS[match(letter, LETTERS) + 1]
            df_sorted[i, "letter"] <- paste(letter, new_letter, sep = "")
            letter <- new_letter
          } else {
            # a-A-a
            df_sorted[i, "letter"] <- letter
          }
        } else { # lowest accuracy
          if (stat_results[stat_results$comb == ab_check, ]$p_value < 0.05) {
            # i and i-1 are significantly different (a-B lettering)
            df_sorted[i, "letter"] <- LETTERS[match(letter, LETTERS) + 1]
          } else {
            # a-A-end
            df_sorted[i, "letter"] <- letter
          }
        }
      }
      # add letters to results dataframe
      df <- merge(
        df,
        df_sorted[c(deparse(predictor), "letter")],
        by = deparse(predictor)
      )
    } else {
      df["letter"] <- "A"
    }

    df
  }

  if (test == "wilcox") {
    stat_result <- withCallingHandlers(
      as.data.frame(
        stats::pairwise.wilcox.test(
          df[[response]],
          df[[predictor]]
        )$p.value
      ),
      warning = function(w) {
        if (grepl("cannot compute exact p-value with ties", w$message)) {
          tryInvokeRestart("muffleWarning")
          as.data.frame(stats::pairwise.wilcox.test(
            df[[response]],
            df[[predictor]],
            exact = FALSE
          )$p.value)
        } else {
          warning(w$message)
        }
      }
    )
    stat_results <- stat_result %>%
      tibble::rownames_to_column(var = "var1") %>%
      tidyr::pivot_longer(
        cols = -c("var1"),
        names_to = "var2",
        values_to = "p_value"
      ) %>%
      dplyr::filter((.data$var1 != .data$var2) & (!is.na(.data$p_value)))
  } else if (test == "t.test") {
    stat_results <- as.data.frame(stats::pairwise.t.test(
      df[[response]],
      df[[predictor]]
    )$p.value) %>%
      tibble::rownames_to_column(var = "var1") %>%
      tidyr::pivot_longer(
        cols = -c("var1"),
        names_to = "var2",
        values_to = "p_value"
      ) %>%
      dplyr::filter((.data$var1 != .data$var2) & (!is.na(.data$p_value)))
  } else {
    warning("We only allow wilcox and t.test for statistical testing.")
  }

  df <- add_letters(
    df[, c(deparse(response), deparse(predictor), "significance")],
    stat_results, test
  )

  list(df = df, stat_results = stat_results)
}



#' Calculate post-selection inference
#'
#' Calculate post-selection inference for glm and lm models
#' @param df
#'  Dataframe with response and predictors as columns.
#' @param regression_model
#'  A starter glm model with a complete term.
#' @param model_type
#'  The model to be used.
#'  Options: glm or lm.
#' @param model_family
#'  A character string describing the family used for model calculation.
#'  See [stats::family] for options.
#' @param term
#'  The call to be used with the model.
#' @param categorical_vars
#'  List of categorical variables within the model
#'  (base names, i.e., names of columns)
#' @param models_overview
#'  Output of [LazyModeler::expand_model_summary()].
#' @param boot_repl
#'  A number or list of bootstrap replicates.
#'  The default is no bootstrapping.
#' @param k
#'  The multiple of the number of degrees of freedom used as
#'  penalty in the model selection. The default k = 2 corresponds to the AIC.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values.
#' @return
#'  Returns a dataframe including the model's independent variables and
#'  their respective raw and expected p-values
#'  as well as the updated glm model
apply_psi <- function(df,
                      regression_model,
                      model_type,
                      model_family,
                      term,
                      categorical_vars,
                      models_overview,
                      boot_repl = 100,
                      k = 2,
                      round_p = 6) {
  # extend categorical into trait columns
  y <- stats::model.response(stats::model.frame(regression_model))
  x <- stats::model.matrix(regression_model)

  dat_num <- data.frame(y = y, x[, -1, drop = FALSE])
  model_args <- list(
    formula = y ~ .,
    data = dat_num
  )

  if (model_type == "glm") {
    model_args[["family"]] <- model_family
    mod_num <- do.call(stats::glm, model_args)
  } else if (model_type == "lm") {
    mod_num <- do.call(stats::lm, model_args)
  }

  # run selcorr
  postsel <- tryCatch(
    selcorr::selcorr(mod_num, boot.repl = boot_repl, k = k, quiet = TRUE),
    error = function(e) {
      warning(paste("selcorr() failed:", e$message))
      NULL
    }
  )

  # extract correct column name for p-values (t or z)
  coef_table <- summary(mod_num)$coefficients
  p_col <- intersect(
    c("Pr(>|t|)", "Pr(>|z|)"),
    colnames(coef_table)
  )
  if (length(p_col) == 0) {
    stop("No valid p-value column found in summary(mod_num)$coefficients.")
  }
  p_raw <- coef_table[, p_col]

  if (is.null(postsel)) {
    warning("Skipping correction: selcorr() failed.")
    result <- data.frame(
      Term = names(p_raw),
      Raw_p = round(p_raw, 6),
      Corrected_p = NA,
      Difference = NA,
      Raw_sig = ifelse(p_raw <= 0.05, "Yes", "No"),
      Corr_sig = NA
    )
    return(result)
  }

  # Extract corrected p-values
  coef_corr <- tryCatch(summary(postsel)$coefficients, error = function(e) NULL)
  if (is.null(coef_corr)) {
    warning("Could not extract corrected p-values from selcorr() output.")
    p_corr <- rep(NA, length(p_raw))
  } else {
    p_col_corr <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_corr))
    if (length(p_col_corr) == 0) p_col_corr <- 1L
    p_corr <- coef_corr[, p_col_corr]
  }

  result <- merge(p_raw, p_corr, by = "row.names", all = TRUE)
  colnames(result) <- c("Coefficient", "Raw_p", "Corrected_p")
  rownames(result) <- result$Coefficient

  # Summarize
  result <- result %>%
    dplyr::mutate(
      Difference = ifelse(
        !is.na(.data$Raw_p) & !is.na(.data$Corrected_p),
        round(.data$Corrected_p - .data$Raw_p, 6),
        NA
      ),
      Raw_p = ifelse(
        !is.na(.data$Raw_p),
        round(.data$Raw_p, 6),
        NA
      ),
      Corrected_p = ifelse(
        !is.na(.data$Corrected_p),
        round(.data$Corrected_p, 6),
        NA
      ),
      Changed_sig = dplyr::case_when(
        (.data$Raw_p <= 0.05) &
          ((.data$Corrected_p > 0.05) | (is.na(.data$Corrected_p))) ~
          "lost significance after correction",
        TRUE ~ dplyr::case_when(
          (.data$Raw_p > 0.05) & (.data$Corrected_p <= 0.05) ~
            "became significant after correction",
          TRUE ~ "no change"
        )
      )
    ) %>%
    dplyr::filter(.data$Coefficient != "(Intercept)")

  # map original term to coefficient in results df
  models_overview_psi <- models_overview %>%
    dplyr::rowwise() %>%
    dplyr::mutate(selcorr_coeff = gsub("\\^|\\(|\\)|\\:", ".",
      .data$predictor,
      perl = TRUE
    )) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  result <- merge(
    result,
    models_overview_psi,
    all.x = TRUE,
    all.y = FALSE,
    by.x = "Coefficient",
    by.y = "selcorr_coeff"
  )

  # determine which coefficients to remove
  # making sure that categorical is evaluated based on smallest p-value
  # and interactions retain main effects
  interactions_to_keep <- result[
    ((is.na(result$Corrected_p)) |
      (result$Corrected_p <= 0.05)) &
      (result$is_interaction),
  ]
  main_effects_to_keep <- unique(c(
    interactions_to_keep$main_effect1,
    interactions_to_keep$main_effect2
  ))

  result <- result %>%
    dplyr::mutate(to_remove = dplyr::case_when(
      .data$pred_type == "categorical" ~ ifelse(
        ((is.na(.data$Corrected_p)) |
          (.data$Corrected_p > 0.05)) &
          ((is.na(.data$smallest_p_value)) |
            (.data$smallest_p_value > 0.05)),
        TRUE,
        FALSE
      ),
      !.data$is_interaction ~ ifelse(
        !(.data$main_effect1 %in% main_effects_to_keep) &
          ((is.na(.data$Corrected_p)) | (.data$Corrected_p > 0.05)),
        TRUE,
        FALSE
      ),
      (is.na(.data$Corrected_p)) | (.data$Corrected_p > 0.05) ~ TRUE,
      TRUE ~ FALSE
    ))

  # remove non-significant coefficients from term
  to_remove_rows <- result[result$to_remove, ] %>%
    dplyr::filter(
      (.data$pred_type != "mixed") |
        (!duplicated(.data$main_effect1, .data$main_effect2))
    )
  single_var_rows <- result[!result$is_interaction, ]
  if (nrow(to_remove_rows) > 0) {
    for (i in seq_len(nrow(to_remove_rows))) {
      i_row <- to_remove_rows[i, ]
      if (i_row$is_interaction) {
        if (i_row$pred_type == "mixed") {
          predictor <- i_row$predictor
          main_effect1 <- i_row$main_effect1
          main_effect2 <- i_row$main_effect2

          if ((main_effect1 %in% single_var_rows$main_effect1) &&
            (single_var_rows[
              single_var_rows$main_effect1 == main_effect1,
              "pred_type"
            ][[1]] == "categorical")) {
            coef_to_remove <- stringr::str_replace(
              predictor,
              sprintf(
                "(?<=\\:|^)[^\\:]*%s[^\\:]*(?=\\:|$)",
                main_effect1
              ),
              main_effect1
            )
          } else if ((main_effect2 %in% single_var_rows$main_effect1) &&
            (single_var_rows[
              single_var_rows$main_effect1 == main_effect2,
              "pred_type"
            ][[1]] == "categorical")) {
            coef_to_remove <- stringr::str_replace(
              predictor,
              sprintf(
                "(?<=\\:|^)[^\\:]*%s[^\\:]*(?=\\:|$)",
                main_effect2
              ),
              main_effect2
            )
          }
        } else {
          coef_to_remove <- i_row$predictor
        }
      } else if (i_row$pred_type == "categorical") {
        coef_to_remove <- i_row$main_effect1
      } else {
        coef_to_remove <- i_row$predictor
      }

      term <- stats::update(
        stats::as.formula(term),
        paste(". ~ . -", coef_to_remove)
      )
    }
  }

  model_args <- list(
    formula = term,
    data = df
  )
  if (model_type == "glm") {
    model_args[["family"]] <- model_family
    final_model <- do.call(stats::glm, model_args)
  } else if (model_type == "lm") {
    final_model <- do.call(stats::lm, model_args)
  }

  result <- result %>%
    dplyr::select(
      tidyr::all_of(
        c("Coefficient", "Raw_p", "Corrected_p", "Difference", "Changed_sig")
      )
    )

  list(
    result = result,
    psi_model = final_model
  )
}


#' Plot PSI p-value changes
#'
#' Plot p-values of model coefficients before and
#'  after post selection inference.
#' @param psi_df
#'  Dataframe provided by [LazyModeler::apply_psi()].
#' @param p_threshold
#'  P-value threshold to be used for significance filtering.
#' @param label_size
#'  Size of test labels within plot.
#' @return
#'  Plot displaying p-values before and after post selection inference
plot_psi <- function(psi_df, p_threshold = 0.05, label_size = 2.5) {
  psi_df$Term <- trimws(as.character(psi_df$Coefficient))

  if ((any(is.na(psi_df$Raw_p))) || (any(is.na(psi_df$Corrected_p)))) {
    warning("Found NA p-values. Replacing with 1.")
  }

  psi_df %>%
    dplyr::mutate(
      Raw_p = ifelse(
        .data$Raw_p == 0,
        1e-12,
        ifelse(is.na(.data$Raw_p), 1, .data$Raw_p)
      ),
      Corrected_p = ifelse(
        .data$Corrected_p == 0,
        1e-12,
        ifelse(is.na(.data$Corrected_p), 1, .data$Corrected_p)
      ),
      Raw_sig = ifelse(.data$Raw_p <= p_threshold,
        "Yes",
        "No"
      ),
      Corr_sig = ifelse(.data$Corrected_p <= p_threshold,
        "Yes",
        "No"
      )
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = psi_df,
      ggplot2::aes(
        x = 1, xend = 2,
        y = log10(.data$Raw_p), yend = log10(.data$Corrected_p),
        color = abs(log10(.data$Corrected_p) - log10(.data$Raw_p))
      ),
      linewidth = 0.8, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = psi_df,
      ggplot2::aes(x = 1, y = log10(.data$Raw_p), fill = .data$Raw_sig),
      shape = 21, color = "black", size = 3
    ) +
    ggplot2::geom_point(
      data = psi_df,
      ggplot2::aes(x = 2, y = log10(.data$Corrected_p), fill = .data$Corr_sig),
      shape = 21, color = "black", size = 3
    ) +
    ggrepel::geom_text_repel(
      data = psi_df,
      ggplot2::aes(x = 1, y = log10(.data$Raw_p), label = .data$Term),
      size = label_size,
      hjust = 1.2,
      nudge_x = -0.05,
      segment.size = 0.25,
      box.padding = 0.2,
      direction = "y",
      max.overlaps = Inf
    ) +
    ggplot2::geom_hline(
      yintercept = log10(p_threshold),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Raw", "Corrected")
    ) +
    ggplot2::scale_color_gradient(
      low = "gray70",
      high = "purple",
      name = "|\u0394 log10(p)|"
    ) +
    ggplot2::scale_fill_manual(
      values = c("Yes" = "#1b9e77", "No" = "#d95f02"),
      name = "Significant"
    ) +
    ggplot2::labs(
      x = "P-Value",
      y = expression(log[10](italic(p))),
      title = "P-values before and after correction"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  p
}
