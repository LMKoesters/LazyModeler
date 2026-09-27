#' Add letters to factor levels
#'
#' Adds letters indicating significant differences between factor levels
#' @param response
#'  String representing the response/y variable (column name)
#' @param predictor
#'  String representing the predictor/x variable (column name)
#' @param m_frame
#'  Dataframe with response and predictors as columns
#' @param stat_result
#'  Results of either wilcox or t.test
#' @return
#'  Updated dataframe with letters
add_letters <- function(response,
                        predictor,
                        m_frame,
                        stat_result) {
  df_sorted <- m_frame |>
    dplyr::group_by(.data[[predictor]]) |>
    dplyr::mutate(resp_mean = mean(.data[[response]])) |>
    dplyr::ungroup() |>
    dplyr::distinct(.data[[predictor]], .keep_all = TRUE) |>
    dplyr::arrange(dplyr::desc(.data$resp_mean)) |>
    dplyr::mutate(idx = dplyr::row_number())

  if (any(stat_result$p_value < 0.05)) {
    stat_result <- prepare_stats_for_letters(stat_result,
                                             df_sorted,
                                             predictor)

    letter <- "A"
    for (i in seq_len(nrow(df_sorted))) {
      ab_check <- paste(i - 1, i, sep = "|")

      if (i == 1) {
        # highest mean response
        df_sorted[1, "letter"] <- letter
      } else if (i != nrow(df_sorted)) {
        # everything between highest and lowest mean response
        aab_check <- stat_result[
          stat_result$comb ==
            paste(
              i,
              i + 1,
              sep = "|"
            ),
        ]
        aabb_check <- stat_result[
          stat_result$comb ==
            paste(
              i - 1,
              i + 1,
              sep = "|"
            ),
        ]

        if (stat_result[stat_result$comb == ab_check, ]$p_value < 0.05) {
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
      } else {
        # lowest mean response
        if (stat_result[stat_result$comb == ab_check, ]$p_value < 0.05) {
          # i and i-1 are significantly different (a-B lettering)
          df_sorted[i, "letter"] <- LETTERS[match(letter, LETTERS) + 1]
        } else {
          # a-A-end
          df_sorted[i, "letter"] <- letter
        }
      }
    }
    m_frame <- m_frame |>
      dplyr::left_join(df_sorted[c(predictor, "letter")],
                       by = predictor)
  } else {
    m_frame["letter"] <- "A"
  }

  m_frame
}

#' Format statistical results
#'
#' Formats results of statistical test for addition of factor level letters
#' @param stat_result
#'  Initial result of statistical test (wilcox or t.test)
#' @param df_sorted
#'  Dataframe with sorted mean responses per factor
#' @param predictor
#'  Predicting factor as string
#' @returns
#'  Formatted statistical result for letter addition
prepare_stats_for_letters <- function(stat_result,
                                      df_sorted,
                                      predictor) {
  stat_result <- merge(
    stat_result,
    df_sorted[c(predictor, "idx")],
    by.x = "var1",
    by.y = predictor
  )
  stat_result <- merge(
    stat_result,
    df_sorted[c(predictor, "idx")],
    by.x = "var2",
    by.y = predictor,
    suffixes = c("1", "2")
  )
  stat_result <- stat_result |>
    dplyr::rowwise() |>
    dplyr::mutate(
      comb = paste(
        sort(c(.data$idx1, .data$idx2)),
        collapse = "|"
      )
    )

  stat_result
}

#' Compute statistics
#'
#' Runs either wilcox or t.test
#' @param m_frame
#'  Data for test
#' @param response
#'  Response variable as string
#' @param predictor
#'  Predicting factor as string
#' @param test
#'  Test to run. Either "wilcox" or "t.test"
#' @param p_threshold
#'  p-value threshold for significance evaluation.
#'    Default: 0.05
#' @returns
#'  List with
#'    a) updated m_frame with letters indicating significant differences
#'    between factor levels
#'    b) result of statistical test
run_stats <- function(m_frame,
                      response,
                      predictor,
                      test = "wilcox",
                      p_threshold = 0.05) {
  if (test == "wilcox") {
    global_test <- "kruskal"
  } else if (test == "t.test") {
    global_test <- "anova"
  } else {
    stop(paste0(
      "We only allow wilcox and t.test for statistical testing.",
      "Please adjust your choice accordingly."
    ))
  }

  c(global_result, global_p) %<-% run_global_test(
    global_test,
    m_frame,
    response,
    predictor
  )

  if (global_p < p_threshold) {
    stat_result <- run_posthoc_test(
      test,
      m_frame,
      response,
      predictor
    )

    m_frame <- add_letters(
      response,
      predictor,
      m_frame[, c(response, predictor)],
      stat_result
    )
  } else {
    stat_result <- NA
    m_frame$letter <- NA_character_
  }

  stat_result <- list(
    global = global_result,
    posthoc = stat_result
  )

  list(m_frame = m_frame,
       stat_result = stat_result)
}

#' Run global statistical test
#'
#' Runs either Kruskal-Wallis or Welch-Anova
#' @param test
#'  Statistical test to run (either kruskal or anova)
#' @param m_frame
#'  Data for statistical test
#' @param response
#'  Response variable as string
#' @param predictor
#'  Predicting factor as string
#' @returns Returns boolean for (non-)significant global result
run_global_test <- function(test,
                            m_frame,
                            response,
                            predictor) {
  if (test == "kruskal") {
    global_test <- stats::kruskal.test(
      m_frame[[response]] ~ m_frame[[predictor]]
    )
    global_p <- global_test$p.value
  } else if (test == "anova") {
    global_test <- stats::oneway.test(
      m_frame[[response]] ~ m_frame[[predictor]],
      var.equal = FALSE
    )
    global_p <- summary(global_test)[[1]][["Pr(>F)"]][1]
  }

  list(global_test, global_p)
}

#' Run statistical test
#'
#' Runs either wilcox or t.test
#' @param test
#'  Statistical test to run (either wilcox or t.test)
#' @param m_frame
#'  Data for statistical test
#' @param response
#'  Response variable as string
#' @param predictor
#'  Predicting factor as string
#' @returns Result of statistical test, pivoted from wide to long
run_posthoc_test <- function(test,
                             m_frame,
                             response,
                             predictor) {
  if (test == "wilcox") {
    stat_result_post <- withCallingHandlers(
      as.data.frame(
        stats::pairwise.wilcox.test(
          m_frame[[response]],
          m_frame[[predictor]]
        )$p.value
      ),
      warning = function(w) {
        if (grepl("cannot compute exact p-value with ties", w$message)) {
          tryInvokeRestart("muffleWarning")
          as.data.frame(
            stats::pairwise.wilcox.test(
              m_frame[[response]],
              m_frame[[predictor]],
              exact = FALSE
            )$p.value
          )
        } else {
          warning(w$message)
        }
      }
    )
  } else if (test == "t.test") {
    stat_result_post <- as.data.frame(
      stats::pairwise.t.test(
        m_frame[[response]],
        m_frame[[predictor]],
        var.equal = FALSE
      )$p.value
    )
  }

  stat_result_post <- stat_result_post |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(
      cols = -c("var1"),
      names_to = "var2",
      values_to = "p_value"
    ) |>
    dplyr::filter((.data$var1 != .data$var2) & (!is.na(.data$p_value)))

  stat_result_post
}
