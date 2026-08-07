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
#' @returns
#'  List with
#'    a) updated m_frame with letters indicating significant differences
#'    between factor levels
#'    b) result of statistical test
run_stats <- function(m_frame,
                      response,
                      predictor,
                      test = "wilcox") {
  if (test == "wilcox") {
    stat_func <- stats::pairwise.wilcox.test
  } else if (test == "t.test") {
    stat_func <- stats::pairwise.t.test
  } else {
    warning("We only allow wilcox and t.test for statistical testing.")
    return(list())
  }

  stat_result <- run_test(stat_func,
                          m_frame,
                          response,
                          predictor)

  m_frame <- add_letters(
    response,
    predictor,
    m_frame[, c(response, predictor)],
    stat_result
  )

  list(m_frame = m_frame,
       stat_result = stat_result)
}

#' Run statistical test
#'
#' Runs either wilcox or t.test
#' @param stat_func
#'  Statistical function to run
#' @param m_frame
#'  Data for statistical test
#' @param response
#'  Response variable as string
#' @param predictor
#'  Predicting factor as string
#' @returns Result of statistical test, pivoted from wide to long
run_test <- function(stat_func,
                     m_frame,
                     response,
                     predictor) {
  stat_result <- withCallingHandlers(
    as.data.frame(
      stat_func(
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
  stat_result <- stat_result |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(
      cols = -c("var1"),
      names_to = "var2",
      values_to = "p_value"
    ) |>
    dplyr::filter((.data$var1 != .data$var2) & (!is.na(.data$p_value)))

  stat_result
}
