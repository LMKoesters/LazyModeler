cor_pivot_longer <- function(correlations, values_to) {
  correlations |>
    tibble::rownames_to_column(var = "coefficientA") |>
    tidyr::pivot_longer(
      !"coefficientA",
      names_to = "coefficientB",
      values_to = values_to
    ) |>
    dplyr::filter(.data$coefficientA != .data$coefficientB)
}

cor_sort_and_filter <- function(correlations_w_p, threshold) {
  correlations_w_p |>
    dplyr::mutate(
      sorted_coefA = pmin(.data$coefficientA, .data$coefficientB),
      sorted_coefB = pmax(.data$coefficientA, .data$coefficientB),
      comparison = paste(.data$sorted_coefA, .data$sorted_coefB)
    ) |>
    dplyr::distinct(.data$comparison, .keep_all = TRUE) |>
    dplyr::select(!tidyr::any_of(c("sorted_coefA",
                                   "sorted_coefB",
                                   "comparison"))) |>
    dplyr::filter(
      ((.data$correlation >= threshold) |
         (.data$correlation <= -threshold)) &
        (.data$p_value < 0.05)
    ) |>
    tibble::add_column(note = NA)
}

cor_prep_autocor <- function(autocorrelations, coefficients) {
  coefficients_df <- data.frame(
    idx = seq_along(coefficients),
    coefficient = coefficients
  )

  autocorrelations <- autocorrelations |>
    dplyr::left_join(coefficients_df,
                     by = dplyr::join_by(
                       .data$coefficientA == .data$coefficient
                     )) |>
    dplyr::left_join(coefficients_df,
                     by = dplyr::join_by(
                       .data$coefficientB == .data$coefficient
                     ),
                     suffix = c("1", "2")) |>
    dplyr::mutate(
      idx_smaller = pmin(.data$idx1, .data$idx2),
      idx_bigger = pmax(.data$idx1, .data$idx2)
    ) |>
    dplyr::arrange(dplyr::desc(.data$idx_bigger)) |>
    dplyr::select(-c("idx1", "idx2"))

  list(coefficients_df, autocorrelations)
}
