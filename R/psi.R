#' Run Post-Selection Inference
#'
#' Runs [selcorr::selcorr()] and creates plot for comparison of raw and
#'  corrected p-values
#' @param final_formula
#'  Formula of the final model selected in [LazyModeler::simplify_model()]
#' @param final_data
#'  Data of the final model selected in [LazyModeler::simplify_model()]
#' @param final_p_values
#'  P-values of the final model selected in [LazyModeler::simplify_model()]
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm" or "glm"
#' @param family
#'  A character string or call describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @param boot_repl
#'  A number or list of psi bootstrap replicates.
#'  Default: 100
#' @param k
#'  The multiple of the number of degrees of freedom used as
#'  penalty in the model selection. The default k = 2 corresponds to the AIC.
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param stat_type
#'  Type of Anova test. Default: 3
#' @param p_threshold
#'  P-value threshold for significance evaluation. Default: 0.05
#' @param label_size
#'  Size of labels within post-selection inference plot.
#'  Default: 2.5
#' @returns
#'  A list with a model overview with raw and corrected p-values,
#'    a plot for comparison of raw and corrected p-values,
#'    and the updated model
run_psi <- function(final_formula,
                    final_data,
                    final_p_values,
                    model_type,
                    family = stats::gaussian,
                    model_args = list(),
                    boot_repl = 100,
                    k = 2,
                    round_p = 6,
                    stat_type = 3,
                    p_threshold = 0.05,
                    label_size = 2.5) {
  final_model <- prepare_num_model(final_formula,
                                   final_data,
                                   model_type,
                                   family,
                                   model_args)

  postsel <- tryCatch(
    selcorr::selcorr(final_model, boot.repl = boot_repl, k = k, quiet = TRUE),
    error = function(e) {
      warning(paste("selcorr() failed:", e$message))
      print(e$message)
      NULL
    }
  )

  if (is.null(postsel)) {
    warning("Skipping correction: selcorr() failed.")
    return(list())
  }

  corrected_p_values <- get_corrected_p_values(final_model,
                                               postsel,
                                               model_type,
                                               stat_type)
  c(updated_formula,
    corrected_p_w_terms) %<-% update_final_formula(final_formula,
                                                   final_data,
                                                   corrected_p_values)
  updated_model <- create_model(updated_formula,
                                final_data,
                                model_type,
                                family,
                                model_args)

  psi_plot <- plot_psi(
    psi_info = corrected_p_w_terms,
    p_threshold = p_threshold,
    label_size = label_size
  )

  list(
    stat_result = corrected_p_values,
    psi_model = updated_model,
    psi_plot = psi_plot
  )
}

#' Update final formula
#'
#' @param final_formula
#'  Formula of the final model selected in [LazyModeler::simplify_model()]
#' @param final_data
#'  Data of the final model selected in [LazyModeler::simplify_model()]
#' @param corrected_p_values
#'  A dataframe with model coefficients, and raw and corrected p-values for
#'    filtering of coefficients within the formula
#' @returns
#'  An updated model formula and the adjusted corrected_p_values dataframe
#'    with factor levels
update_final_formula <- function(final_formula,
                                 final_data,
                                 corrected_p_values) {
  terms <- corrected_p_values$predictor
  interactions <- extract_interactions(terms)
  cat_vars <- extract_categorical_vars(stats::model.frame(final_data))

  factor_info <- interactions |>
    dplyr::left_join(cat_vars, by = c("main_effect" = "trait")) |>
    dplyr::mutate(main_effect = ifelse(is.na(.data$var),
                                       .data$main_effect,
                                       .data$var)) |>
    dplyr::group_by(.data$predictor) |>
    dplyr::filter(!all(is.na(.data$var))) |>
    dplyr::mutate(
      numeric_effect = dplyr::first(.data$main_effect[is.na(.data$var)])
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$var)) |>
    dplyr::mutate(
      term = paste(.data$numeric_effect, .data$var, sep = ":")
    ) |>
    dplyr::select(-c("var", "main_effect", "numeric_effect")) |>
    rbind(cat_vars |> dplyr::rename(predictor = "trait",
                                    term = "var"))

  terms_to_remove <- corrected_p_values |>
    dplyr::select(c("predictor", "corrected_p_value")) |>
    dplyr::left_join(factor_info, by = "predictor") |>
    dplyr::mutate(term = ifelse(is.na(.data$term),
                                .data$predictor,
                                .data$term)) |>
    dplyr::group_by(.data$term) |>
    dplyr::filter(all(is.na(.data$corrected_p_value) |
                        .data$corrected_p_value >= 0.1))

  updated_formula <- final_formula
  for (term in terms_to_remove$predictor) {
    d <- paste(". ~ . -", term)
    updated_formula <- stats::update(stats::as.formula(updated_formula), d)
  }

  corrected_p_w_terms <- corrected_p_values |>
    dplyr::left_join(factor_info, by = "predictor") |>
    dplyr::mutate(term = ifelse(is.na(.data$term),
                                .data$predictor,
                                .data$term))

  list(updated_formula, corrected_p_w_terms)
}

#' Get corrected p-values
#'
#' @param final_model
#'  The final model selected in [LazyModeler::simplify_model()]
#' @param postsel
#'  The model selected by [selcorr::selcorr()]
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm" or "glm"
#' @param stat_type
#'  Type of Anova test. Default: 3
#' @returns
#'  A dataframe with model coefficients, and raw and corrected p-values for
#'    filtering of coefficients within the formula
get_corrected_p_values <- function(final_model,
                                   postsel,
                                   model_type,
                                   stat_type = 3) {
  corrected_p_value <- get_model_p_values(postsel, model_type, stat_type)
  p_col <- colnames(corrected_p_value)[[
    grep("Pr\\(", colnames(corrected_p_value))
  ]]
  corrected_p_value <- corrected_p_value |>
    dplyr::rename(corrected_p_value = dplyr::all_of(p_col)) |>
    as.data.frame() |>
    dplyr::select(c("predictor", "corrected_p_value"))

  raw_p_value <- get_model_p_values(final_model, model_type, stat_type)
  p_col <- colnames(raw_p_value)[[
    grep("Pr\\(", colnames(raw_p_value))
  ]]
  p_values <- raw_p_value |>
    dplyr::rename(raw_p_value = dplyr::all_of(p_col)) |>
    as.data.frame() |>
    dplyr::full_join(corrected_p_value,
                     by = "predictor") |>
    dplyr::filter(!is.na(.data$raw_p_value)) |>
    dplyr::mutate(
      diff = ifelse(
        !is.na(.data$raw_p_value) & !is.na(.data$corrected_p_value),
        round(.data$corrected_p_value - .data$raw_p_value, 6),
        NA
      ),
      significance_changed = dplyr::case_when(
        (.data$raw_p_value <= 0.05) &
          ((.data$corrected_p_value > 0.05) |
             (is.na(.data$corrected_p_value))) ~
          "lost significance after correction",
        TRUE ~ dplyr::case_when(
          (.data$raw_p_value > 0.05) & (.data$corrected_p_value <= 0.05) ~
            "became significant after correction",
          TRUE ~ "no change"
        )
      )
    )

  p_values
}

#' Prepare numeric model
#'
#' @param final_formula
#'  Formula of the final model selected in [LazyModeler::simplify_model()]
#' @param final_data
#'  Data of the final model selected in [LazyModeler::simplify_model()]
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm" or "glm"
#' @param family
#'  A character string or call describing the family used for model calculation.
#'  See [stats::family] for options. Default: gaussian.
#' @param model_args
#'  A named list of additional arguments given directly to model call
#' @returns
#'  A numeric-only version of the final model, where categorical variables
#'    are split into separate columns
prepare_num_model <- function(final_formula,
                              final_data,
                              model_type,
                              family = stats::gaussian,
                              model_args = list()) {
  final_formula <- check_formula(final_formula, final_data, add_main = TRUE)

  full_model <- create_model(
    final_formula,
    final_data,
    model_type,
    family,
    model_args
  )

  model_overview <- stats::coef(summary(full_model)) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "predictor")
  c(categorical_vars,
    term_factors,
    numeric_var) %<-% get_term_factors(full_model,
                                       model_type,
                                       model_overview)
  categorical_vars <- categorical_vars |>
    dplyr::rename(trait_term = "predictor") |>
    dplyr::mutate(predictor = .data$trait_term,
                  main_effect = .data$term)
  term_factors <- term_factors |>
    dplyr::rename(trait_term = "interaction") |>
    dplyr::bind_rows(categorical_vars[, c("term", "trait_term",
                                          "predictor", "main_effect")])

  data <- stats::model.frame(full_model)
  m_matrix <- stats::model.matrix(full_model)
  m_matrix <- m_matrix[, attr(m_matrix, "assign") != 0L, drop = FALSE] |>
    as.data.frame()
  formula <- stats::formula(full_model)

  cat_main_added <- c()
  for (term in unique(term_factors$term)) {
    d <- paste(". ~ . -", term)
    formula <- stats::update(stats::as.formula(formula), d)
    if (term %in% colnames(data)) data <- dplyr::select(data,
                                                        -dplyr::all_of(c(term)))

    cat_main_effects <- term_factors$main_effect[
      (term_factors$term == term) &
        (term_factors$predictor != term_factors$main_effect)
    ]
    for (cat_main in cat_main_effects) {
      if (!cat_main %in% cat_main_added) {
        factor_formula <- stats::as.formula(paste0("~ 0 + ", cat_main))
        factor_matrix <- stats::model.matrix(factor_formula, data = final_data)
        data <- data |>
          dplyr::bind_cols(factor_matrix)
        cat_main_added <- append(cat_main_added, cat_main)
      }
    }

    for (trait in unique(term_factors$trait_term[term_factors$term == term])) {
      d <- paste(". ~ . +", trait)
      formula <- stats::update(stats::as.formula(formula), d)
    }
  }

  num_model <- create_model(
    formula,
    data,
    model_type,
    family,
    model_args
  )

  num_model
}
