#' Handle autocorrelations
#'
#' Function to process autocorrelations between columns in input dataframe
#' @param data
#'  Dataframe with response and predictors as columns.
#' @param cols
#'  Columns to check for autocorrelations. The order of columns dictates
#'  priority basis for removal of predictors. Columns further down the list
#'  are removed first.
#' @param remove
#'  Boolean switch for automatic removal of autocorrelated variables.
#' @param threshold
#'  The threshold at which two variables are to be considered as autocorrelated.
#' @param cor_args
#'  Further arguments for [stats::cor()].
#' @return
#'  Named list with
#'  a) a vector with removed predictors (NULL if none were removed), and
#'  b) a dataframe with comprehensive information on autocorrelations
handle_autocorrelations <- function(
  data,
  cols,
  remove = TRUE,
  threshold = 0.7,
  cor_args = list(method = c("pearson"),
    use = "complete.obs"
  )
) {
  # SETUP
  check_correlation_threshold(threshold)

  # COMPUTE CORRELATIONS
  cor_args <- c(cor_args, list(x = data[, cols]))
  correlations <- as.data.frame(do.call(stats::cor, cor_args))

  # COMPUTE P-VALUES
  correlations_p_val <- as.data.frame(
    corrplot::cor.mtest(data[, cols])$p
  )

  # PIVOT LONGER & MERGE
  correlations_l <- cor_pivot_longer(correlations, "correlation")
  correlations_p_val_l <- cor_pivot_longer(correlations_p_val, "p_value")
  correlations_w_p <- merge(
    correlations_l,
    correlations_p_val_l,
    by = c("coefficientA", "coefficientB")
  )

  correlations_w_p <- cor_sort_and_filter(correlations_w_p, threshold)
  if (nrow(correlations_w_p) == 0) {
    list("autocorrelations_info" = NULL,
         "removed_predictors" = NULL)
  } else if (remove) {
    c(autocorrelations,
      removed_predictors) %<-% remove_autocorrelations(correlations_w_p, cols)
    autocorrelations <- autocorrelations[, c("coefficientA",
                                             "coefficientB",
                                             "correlation",
                                             "p_value",
                                             "note")]
    list("autocorrelations_info" = autocorrelations,
         "removed_predictors" = removed_predictors)
  } else {
    warning(
      paste("Some of your variables are autocorrelated.",
            "Check $autocorrelations for more info.",
            collapse = " ")
    )
    list("autocorrelations_info" = correlations_w_p,
         "removed_predictors" = NULL)
  }
}

#' Handle autocorrelations
#'
#' Function to remove autocorrelations between columns in input dataframe
#' @param correlations_w_p
#'  Dataframe of sorted and indexed autocorrelated variables with p-values.
#' @param cols
#'  Columns to check for autocorrelations. The order of columns dictates
#'  priority basis for removal of predictors. Columns further down the list
#'  are removed first.
#' @return
#'  A list with
#'  a) a vector with removed predictors (NULL if none were removed), and
#'  b) a dataframe with comprehensive information on autocorrelations
remove_autocorrelations <- function(correlations_w_p, cols) {
  c(coefficients, autocorrelations) %<-% cor_prep_autocor(correlations_w_p,
                                                          cols)
  removed_predictors <- vector()
  for (i in seq_len(nrow(autocorrelations))) {
    # C is the least important and part of comparison
    # B is the more important of comparison
    # A is more important than B

    autocor_row <- autocorrelations[i, ]
    c <- autocor_row$idx_bigger
    b <- autocor_row$idx_smaller

    b_to_a <- autocorrelations[autocorrelations$idx_bigger == b, ]
    if (nrow(b_to_a) > 0) {
      # check if there's at least 1 A
      for (a in b_to_a$idx_smaller) {
        # for each variable that is correlated to and more important than B
        coefficient_b <- coefficients[b, "coefficient"]
        coefficient_c <- coefficients[c, "coefficient"]

        a_b_c <- autocorrelations[
          (autocorrelations$idx_bigger == autocor_row$idx_bigger) &
            (autocorrelations$idx_smaller == a),
        ]
        if (nrow(a_b_c) == 0) {
          # A!=C but A==B and B==C: remove B
          coefficient_a <- coefficients[a, "coefficient"]
          removed_predictors <- append(
            removed_predictors,
            coefficient_b
          )
          autocorrelations[i, "note"] <- sprintf(
            "%s==%s and %s==%s, but %s!=%s; removed %s",
            coefficient_a, coefficient_b,
            coefficient_b, coefficient_c,
            coefficient_a, coefficient_c,
            coefficient_b
          )
        } else {
          removed_predictors <- append(
            removed_predictors,
            coefficient_c
          )
          autocorrelations[i, "note"] <- sprintf("removed %s", coefficient_c)
        }
      }
    } else {
      coefficient_c <- coefficients[c, "coefficient"]
      removed_predictors <- append(removed_predictors, coefficient_c)
      autocorrelations[i, "note"] <- sprintf("removed %s", coefficient_c)
    }
  }

  list(autocorrelations, removed_predictors)
}
