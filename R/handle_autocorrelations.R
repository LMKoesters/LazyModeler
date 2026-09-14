#' Handle autocorrelations
#'
#' Processes autocorrelations between columns in input dataframe by calculating
#'  correlation scores, determining which column best to remove if
#'  more than two columns are autocorrelated, and removing related terms from
#'  a given formula.
#' @param formula
#'  A formula used for downstream model creation and simplification
#' @param data
#'  Underlying data for autocorrelation detection and downstream
#'    model creation
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "lmer", "glmer",
#'  "nlme", "gam", and "nls"
#' @param family
#'  A character string or call describing the family used for model calculation.
#'    See [stats::family] for options.
#' @param cols
#'  Columns to check for autocorrelations. The order of columns dictates
#'  priority basis for removal of predictors. Columns further down the list
#'  are removed first
#' @param remove
#'  Boolean switch for automatic removal of autocorrelated variables.
#'    Default: TRUE
#' @param threshold
#'  The threshold at which two variables are to be considered autocorrelated
#' @param p_threshold
#'  p-value threshold for significance evaluation
#' @param cor_args
#'  Further arguments for [stats::cor()] and [stats::cor.test()].
#'    Default: method = "pearson" and use = "complete.obs"
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @return
#'  Named list with
#'  a) a vector with removed predictors, and
#'  b) a dataframe with comprehensive information on autocorrelations;
#'    NULL if no autocorrelations were detected
#'  c) an updated formula without autocorrelated variables
handle_autocorrelations <- function(
    formula,
    data,
    model_type,
    family = NA,
    cols = c(),
    remove = TRUE,
    threshold = 0.7,
    p_threshold = 0.05,
    cor_args = list(method = c("pearson"),
                    use = "complete.obs"),
    model_args = list()) {
  if (length(cols) == 1) {
    stop(
      paste(
        "You have provided only a single column for autocorrelation",
        "testing. You need to specify at least two columns.",
        "Alternatively, you can let us determine relevant columns",
        "based on the model formula."
      )
    )
  }

  # SETUP
  check_correlation_threshold(threshold)
  cor_args <- check_cor_args(cor_args)

  if (length(cols) > 0) {
    terms <- attr(stats::terms.formula(formula), "term.labels")
    cols <- cols[cols %in% terms]
    cor_formula <- stats::reformulate(cols)
  } else {
    cor_formula <- formula
  }

  # FORMAT DATA AND DETERMINE COLUMNS TO TEST
  c(data_ext, term_map_all) %<-% format_cor_data(
    data,
    cor_formula,
    model_type,
    family,
    model_args
  )

  if (length(cols) > 0) {
    term_map_cor <- term_map_all[term_map_all$term %in% cols, ]
    term_map_cor <- term_map_cor[order(match(term_map_cor$term, cols)), ]
  } else {
    term_map_cor <- term_map_all
  }

  if (length(unique(term_map_cor$term)) < 2) {
    stop(
      sprintf(
        paste(
          "We were unable to detect enough valid columns for testing against",
          "autocorrelations. Please check whether the columns included in your",
          "dataframe and formula match. Columns found were: %s"
        ),
        paste(term_map_cor$term, collapse = " ")
      )
    )
  }

  c(term_map_cor, has_no_variance) %<-% check_variance(data_ext, term_map_cor)
  term_map_cor <- sort_term_map(term_map_cor)
  term_map_all <- sort_term_map(term_map_all)
  cols <- unique(term_map_cor$column)

  # COMPUTE CORRELATIONS
  cor_args <- c(cor_args, list(x = data_ext[, cols]))
  if (cor_args$use == "all.obs" && anyNA.data.frame(cor_args$x)) {
    stop(
      paste(
        "You specified use=all.obs for autocorrelation detection",
        "but your data contains NAs. Please check and run again."
      )
    )
  } else if (cor_args$use == "complete.obs") {
    complete <- stats::complete.cases(cor_args$x)
    if (!any(complete)) {
      stop(
        paste(
          "You have specified use=complete.obs, but your data does not",
          "contain complete cases. Please check and run again."
        )
      )
    }
  }
  correlations <- as.data.frame(
    do.call(stats::cor, cor_args[names(cor_args) %in% c("x", "use", "method")])
  )
  correlations_l <- cor_pivot_longer(correlations, "correlation") |>
    dplyr::filter(.data$correlation >= threshold |
                    .data$correlation <= -threshold)
  if (nrow(correlations_l) == 0) {
    out <- list("autocorrelations_info" = NULL,
                "problematic_predictors" = c(has_no_variance),
                "formula" = formula)
    if (remove) {
      removable_terms <- removed_preds_to_terms(out$problematic_predictors,
                                                term_map_cor)
      out$formula <- remove_autocor_predictors(formula,
                                               unique(removable_terms))
    } else {
      out$formula <- formula
    }
    return(out)
  }

  # COMPUTE P-VALUES
  ## remove NAs
  if (cor_args$use == "complete.obs" || cor_args$use == "na.or.complete") {
    complete <- stats::complete.cases(cor_args$x)
    cor_args$x <- cor_args$x[complete, ]
  }
  correlations_w_p <- extract_cor_p_values(
    correlations_l,
    cor_args
  )

  # SORT & FILTER BY P-THRESHOLD
  correlations_w_p <- cor_sort_and_filter(
    correlations_w_p,
    threshold,
    p_threshold
  )

  if (nrow(correlations_w_p) == 0) {
    out <- list("autocorrelations_info" = NULL,
                "problematic_predictors" = c(has_no_variance))
  } else if (remove) {
    c(autocorrelations,
      problematic_predictors) %<-% remove_autocorrelations(
      correlations_w_p,
      term_map_all,
      term_map_cor
    )
    problematic_predictors <- c(problematic_predictors, has_no_variance)
    autocorrelations <- autocorrelations[, c("coefficientA",
                                             "coefficientB",
                                             "correlation",
                                             "p_value",
                                             "note")]
    out <- list("autocorrelations_info" = autocorrelations,
                "problematic_predictors" = problematic_predictors)
  } else {
    warning(
      paste("Some of your variables are autocorrelated.",
            "Check $autocorrelations for more info.",
            collapse = " ")
    )
    out <- list("autocorrelations_info" = correlations_w_p,
                "problematic_predictors" = c(has_no_variance))
  }

  # UPDATE FORMULA
  if (remove) {
    removable_terms <- removed_preds_to_terms(out$problematic_predictors,
                                              term_map_cor)
    out$formula <- remove_autocor_predictors(formula,
                                             unique(removable_terms))
  } else {
    out$formula <- formula
  }

  out
}

#' Handle autocorrelations
#'
#' Removes autocorrelations between columns of an input dataframe.
#' @param correlations_w_p
#'  Dataframe of sorted and indexed autocorrelated variables with p-values.
#' @param term_map_all
#'  Terms with columns represented in the data.
#' @param term_map_cor
#'  Terms with columns to check for autocorrelations. The order of columns
#'  dictates priority basis for removal of predictors. Columns further down
#'  the list are removed first.
#' @return
#'  A list with
#'  a) a vector with removed predictors (NULL if none were removed), and
#'  b) a dataframe with comprehensive information on autocorrelations
remove_autocorrelations <- function(
    correlations_w_p,
    term_map_all,
    term_map_cor) {
  cols <- unique(term_map_cor$column)
  c(coefficients, autocorrelations) %<-% cor_prep_autocor(correlations_w_p,
                                                          cols)
  interactions <- term_map_all[term_map_all$is_interaction, ]
  term_map_distinct <- term_map_cor[, c("column", "term")] |>
    dplyr::distinct(.keep_all = TRUE)
  autocorrelations <- autocorrelations |>
    dplyr::left_join(term_map_distinct,
                     by = c("coefficientA" = "column")) |>
    dplyr::left_join(term_map_distinct,
                     by = c("coefficientB" = "column"),
                     suffix = c("A", "B")) |>
    dplyr::filter(.data$termA != .data$termB)

  autocors_int <- autocorrelations[
    (autocorrelations$coefficientA %in% interactions$column) &
      (autocorrelations$coefficientB) %in% interactions$column,
  ]
  autocors_int_me <- autocorrelations[
    ((autocorrelations$coefficientA %in% interactions$column) |
       (autocorrelations$coefficientB %in% interactions$column)) &
      (!rownames(autocorrelations) %in% rownames(autocors_int)),
  ]
  autocors_me <- autocorrelations[
    (!autocorrelations$coefficientA %in% interactions$column) &
      (!autocorrelations$coefficientB %in% interactions$column),
  ]

  problematic_predictors <- c()
  autocors_notes <- c()
  for (autocors in list(autocors_me, autocors_int_me, autocors_int)) {
    # first: remove main effects and blocking interactions in one step
    c(autocors_notes_i,
      problematic_predictors) %<-% determine_removable_predictors(
      autocors,
      coefficients,
      term_map_all,
      interactions,
      problematic_predictors
    )
    autocors_notes <- rbind(autocors_notes, autocors_notes_i)
  }

  list(autocors_notes, problematic_predictors)
}

#' Determine predictors to remove
#'
#' Decides which autocorrelated predictors should be removed while
#'  considering interactions and main effects
#' @param autocorrelations
#'  A dataframe with indexed variables and a dataframe with pairs of
#'    autocorrelated variables.
#' @param coefficients
#'  A dataframe with indexed variables.
#' @param term_map
#'  Columns to check for autocorrelations. The order of columns dictates
#'    priority basis for removal of predictors. Columns further down the list
#'    are removed first.
#' @param interactions
#'  Map of columns to terms with only interactions.
#' @param problematic_predictors
#'  A character vector with columns so far removed due to autocorrelations.
#' @return
#'  A list with
#'  a) a vector with removed predictors (NULL if none were removed), and
#'  b) a dataframe with comprehensive information on autocorrelations
determine_removable_predictors <- function(
    autocorrelations,
    coefficients,
    term_map,
    interactions,
    problematic_predictors = c()) {
  # NOTE: the smaller the index, the more important the coefficient
  for (i in seq_len(nrow(autocorrelations))) {
    # C is the least important and part of comparison
    # B is the more important of comparison
    # A is more important than B
    autocor_row <- autocorrelations[i, ]
    c <- autocor_row$idx_bigger
    coefficient_c <- coefficients[c, "coefficient"]
    b <- autocor_row$idx_smaller
    coefficient_b <- coefficients[b, "coefficient"]

    b_to_a <- autocorrelations[autocorrelations$idx_bigger == b, ]
    # check if there's at least 1 A
    if (nrow(b_to_a) > 0) {
      # for each variable that is correlated to and more important than B
      for (a in b_to_a$idx_smaller) {
        coefficient_a <- coefficients[a, "coefficient"]
        # check A!=C
        a_b_c <- autocorrelations[
          (autocorrelations$idx_bigger == autocor_row$idx_bigger) &
            (autocorrelations$idx_smaller == a),
        ]
        if (nrow(a_b_c) == 0) {
          # A!=C but A==B and B==C: remove B
          if (coefficient_b %in% problematic_predictors ||
                (coefficient_a %in% problematic_predictors &&
                   coefficient_c %in% problematic_predictors)) {
            coefficient_to_remove <- NA
            already_removed <- paste(coefficient_b, "was already removed")
          } else if (coefficient_a %in% problematic_predictors) {
            coefficient_to_remove <- coefficient_c
            already_removed <- paste(coefficient_a, "was already removed")
          } else if (coefficient_c %in% problematic_predictors) {
            coefficient_to_remove <- coefficient_b
            already_removed <- paste(coefficient_c, "was already removed")
          } else {
            coefficient_to_remove <- coefficient_b
            already_removed <- ""
          }

          if (!is.na(coefficient_to_remove)) {
            problematic_predictors <- append(
              problematic_predictors,
              coefficient_to_remove
            )

            # check for interaction with coefficient_b as main effect
            blocking_interactions <- unique(
              interactions[interactions$column_main == coefficient_to_remove,
                           "column"]
            )
            if (length(blocking_interactions) > 0) {
              problematic_predictors <- append(
                problematic_predictors,
                blocking_interactions
              )

              autocorrelations[i, "note"] <- sprintf(
                paste("%s==%s and %s==%s, but %s!=%s; removed %s",
                      already_removed,
                      "and blocking interaction(s): %s"),
                coefficient_a, coefficient_b,
                coefficient_b, coefficient_c,
                coefficient_a, coefficient_c,
                coefficient_b, paste(blocking_interactions, collapse = " ")
              )
            } else {
              autocorrelations[i, "note"] <- sprintf(
                "%s==%s and %s==%s, but %s!=%s; removed %s",
                coefficient_a, coefficient_b,
                coefficient_b, coefficient_c,
                coefficient_a, coefficient_c,
                coefficient_b
              )
            }
          }
        } else if (!((coefficient_c %in% problematic_predictors) ||
                       (coefficient_b %in% problematic_predictors))) {
          problematic_predictors <- append(
            problematic_predictors,
            coefficient_c
          )

          # check whether there's interaction with coefficient_c as main effect
          blocking_interactions <- unique(
            interactions[interactions$column_main == coefficient_c, "column"]
          )
          if (length(blocking_interactions) > 0) {
            problematic_predictors <- append(
              problematic_predictors,
              blocking_interactions
            )

            autocorrelations[i, "note"] <- sprintf(
              "removed %s and blocking interaction(s): %s",
              coefficient_c, paste(blocking_interactions, collapse = " ")
            )
          } else {
            autocorrelations[i, "note"] <- sprintf("removed %s", coefficient_c)
          }
        } else {
          if (coefficient_c %in% problematic_predictors &&
                coefficient_b %in% problematic_predictors) {
            already_removed <- paste(coefficient_b, coefficient_c, sep = " & ")
          } else if (coefficient_b %in% problematic_predictors) {
            already_removed <- coefficient_b
          } else {
            already_removed <- coefficient_c
          }
          autocorrelations[i, "note"] <- sprintf("already removed %s",
                                                 already_removed)
        }
      }
    } else if (!((coefficient_c %in% problematic_predictors) ||
                   (coefficient_b %in% problematic_predictors))) {
      problematic_predictors <- append(
        problematic_predictors,
        coefficient_c
      )

      # check whether there's interaction with coefficient_c as main effect
      blocking_interactions <- unique(
        interactions[interactions$column_main == coefficient_c, "column"]
      )
      if (length(blocking_interactions) > 0) {
        problematic_predictors <- append(
          problematic_predictors,
          blocking_interactions
        )

        autocorrelations[i, "note"] <- sprintf(
          "removed %s and blocking interaction(s): %s",
          coefficient_c, paste(blocking_interactions, collapse = " ")
        )
      } else {
        autocorrelations[i, "note"] <- sprintf("removed %s", coefficient_c)
      }
    } else {
      if (coefficient_c %in% problematic_predictors &&
            coefficient_b %in% problematic_predictors) {
        already_removed <- paste(coefficient_b, coefficient_c, sep = " & ")
      } else if (coefficient_b %in% problematic_predictors) {
        already_removed <- coefficient_b
      } else {
        already_removed <- coefficient_c
      }
      autocorrelations[i, "note"] <- sprintf("already removed %s",
                                             already_removed)
    }
  }

  list(autocorrelations, problematic_predictors)
}

#' Computes p-values for correlations
#'
#' Uses [stats::cor.test()] to compute p-values of correlations
#'  detected using [stats::cor()].
#' @param correlations_l
#'  A dataframe with autocorrelations, but pivoted longer
#' @param cor_args
#'  Further arguments for [stats::cor()] and [stats::cor.test()].
#'    Default: method = "pearson" and use = "complete.obs"
#' @return
#'  An updated dataframe of autocorrelations with added p-values
extract_cor_p_values <- function(correlations_l, cor_args) {
  for (i in seq_len(nrow(correlations_l))) {
    cor_args_cp <- cor_args[!names(cor_args) %in% c("use")]
    coef_a <- correlations_l[i, "coefficientA"][[1]]
    coef_b <- correlations_l[i, "coefficientB"][[1]]
    cor_args_cp$y <- cor_args_cp$x[, coef_b]
    cor_args_cp$x <- cor_args_cp$x[, coef_a]
    if ((cor_args$use == "everything") &&
          (anyNA(cor_args$x) || anyNA(cor_args$y))) {
      correlations_l[i, "p_value"] <- NA_real_
    } else if (length(cor_args_cp$x) == 0) {
      correlations_l[i, "p_value"] <- NA_real_
    } else {
      cor_stats <- do.call(
        stats::cor.test,
        cor_args_cp
      )
      
      correlations_l[i, "p_value"] <- cor_stats$p.value
    }
  }

  correlations_l
}
