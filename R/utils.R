#' Pivots correlations to long format
#'
#' @param correlations
#'  The result of [corrplot::cor.mtest()] as a dataframe
#' @param values_to
#'  Column name for values
#' @returns
#'  A pivoted dataframe with correlations
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

#' Sort and filter correlations
#'
#' @param correlations_w_p
#'  Dataframe with information on correlations between variables with p-values
#' @param threshold
#'  The threshold at which to consider two variables autocorrelated
#' @returns
#'  A sorted dataframe of autocorrelations with a correlation coefficient
#'  equal to or larger than threshold
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

#' Format correlation information for removal of variables
#'
#' @param autocorrelations
#'  Dataframe with autocorrelated variables with p-values
#' @param coefficients
#'  List of variables to consider for removal
#' @returns
#'  A dataframe with indexed variables and a dataframe with pairs of
#'    autocorrelated variables
cor_prep_autocor <- function(autocorrelations, coefficients) {
  coefficients_df <- data.frame(
    idx = seq_along(coefficients),
    coefficient = coefficients
  )

  autocorrelations <- autocorrelations |>
    dplyr::left_join(coefficients_df, by = c("coefficientA" = "coefficient")) |>
    dplyr::left_join(coefficients_df, by = c("coefficientB" = "coefficient"),
                     suffix = c("1", "2")) |>
    dplyr::mutate(
      idx_smaller = pmin(.data$idx1, .data$idx2),
      idx_bigger = pmax(.data$idx1, .data$idx2)
    ) |>
    dplyr::arrange(dplyr::desc(.data$idx_bigger)) |>
    dplyr::select(-c("idx1", "idx2"))

  list(coefficients_df, autocorrelations)
}

#' Determine appropriate model family
#'
#' @param data
#'  Data to consider when matching model family
#' @param lhs
#'  Left-hand-side of model formula (i.e., response term)
#' @returns
#'  A character vector with valid families given the data and
#'    a notice explaining the choice
determine_model_family <- function(data, lhs) {
  response_col <- paste(lhs, collapse = "")

  if (grepl("cbind", response_col)) {
    return(check_cbind_model_family(data, paste(lhs)))
  }

  if (!response_col %in% colnames(data)) {
    message("Transformed response for model family check.")
    data[response_col] <- with(data, eval(parse(text = lhs)))
  }

  response_data <- data[!is.na(data[[response_col]]), response_col]

  check_response_data_format(response_col, response_data)

  is_num <- is.numeric(response_data)

  if (is.logical(response_data) ||
        all(response_data %in% c(0, 1))) {
    valid_families <- c("binomial")
    notice <- paste("For your data, binomial appears to be an apprioriate",
                    "distribution. If observations represent grouped binomial",
                    "outcomes, provide trials/weights or use",
                    "cbind(successes, failures). Please also check for",
                    "overdispersion; if present, quasibinomial may be more",
                    "appropriate",
                    collape = " ")
  } else if (is_num && (min(response_data) >= 0 &&
                          max(response_data) <= 1)) {
    valid_families <- c("gaussian", "quasibinomial")
    notice <- paste("Your values appear to be proportions.",
                    "quasibinomial may be appropriate",
                    "if values represent proportions from",
                    "binomial trials and trial sizes are supplied as",
                    "weights. Otherwise gaussian may be more",
                    "appropriate.", collapse = " ")
  } else if (min(response_data) >= 0 && all(response_data %% 1 == 0)) {
    valid_families <- c("poisson")
    notice <- paste("Poisson assumes non-negative count data",
                    "and mean-variance equality. Check for",
                    "overdispersion; if present, quasipoisson or",
                    "negative binomial may be more appropriate.",
                    "Note, however, that we currently only",
                    "cover poisson distributions.", collapse = " ")
  } else {
    valid_families <- c("gaussian")
    notice <- paste("Gaussian is suggested for continuous outcomes.",
                    "Consider transformations or other families if",
                    "residual diagnostics indicate non-normality,",
                    "heteroscedasticity, or bounded support.",
                    collapse = " ")
  }

  list(valid_families = valid_families,
       notice = notice)
}

#' Family to character
#'
#' @param family
#'  The family to turn into a string
get_family_character <- function(family) {
  if (inherits(family, "family")) {
    family <- family$family
  } else if (typeof(family) == "closure") {
    family <- family()$family
  }

  family
}

#' Extract interactions
#'
#' @param predictors
#'  A list of predictors extracted from a formula
#' @returns
#'  A long dataframe with interactions and their separated effects
extract_interactions <- function(predictors) {
  predictor_info <- data.frame(list(predictor = predictors)) |>
    dplyr::mutate(
      is_interaction = dplyr::case_when(
        stringr::str_detect(.data$predictor, "(^[^(].+:)|(:.+[^)]$)") ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::filter(.data$is_interaction) |>
    dplyr::select(c("predictor")) |>
    dplyr::mutate(main_effect = stringr::str_split(.data$predictor, ":")) |>
    tidyr::unnest_longer("main_effect") |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

  predictor_info
}

#' Compare models
#'
#' @param evaluation_methods
#'  Character vector with methods to use for model evaluation.
#'  Allowed evaluation methods: "aic", "aicc", "bic", or "anova".
#'  Default: c("anova")
#' @param models
#'  List with models to compare
#' @param direction
#'  Direction of model selection used for anova interpretation
#' @param old_model_assessments
#'  The evaluation results of the previous model
#' @param delta
#'  Minimal value by which the two models must differ for the newer model to
#'    be rejected
#' @returns
#'  Information on whether or not the model has improved and the results
#'    of the model's evaluation
compare_models <- function(evaluation_methods,
                           models,
                           direction,
                           old_model_assessments,
                           delta = 2) {
  if (direction == "backward") models_idx <- list("small" = 1, "large" = 2)
  else models_idx <- list("small" = 2, "large" = 1)

  evaluators <- list(
    "aic" = stats::AIC,
    "aicc" = MuMIn::AICc,
    "bic" = stats::BIC
  )

  evaluation_results <- list()
  for (eval_m in evaluation_methods[evaluation_methods != "anova"]) {
    res <- evaluators[[eval_m]](models[[1]])
    evaluation_results[eval_m] <- res
    diff <- res - old_model_assessments[[eval_m]]
    if (diff > delta) return(
      list("has_improved" = FALSE,
           "assessments" = evaluation_results)
    )
  }

  if ("anova" %in% evaluation_methods) {
    anova_res <- stats::anova(models[[models_idx$small]],
                              models[[models_idx$large]])
    p_col <- colnames(anova_res)[[grep("Pr\\(", colnames(anova_res))]]
    if ((direction == "backward") &&
          (anova_res[2, p_col] < .05)) {
      return(list("has_improved" = FALSE,
                  "assessments" = evaluation_results))
    } else if ((direction == "forward") &&
                 (anova_res[2, p_col] >= .05)) {
      return(list("has_improved" = FALSE,
                  "assessments" = evaluation_results))
    }
  }

  list("has_improved" = TRUE,
       "assessments" = evaluation_results)
}

#' Get p-values of model coefficients
#'
#' @param model
#'  A model to get p-values for
#' @param model_type
#'  The type of model to get p-values for. Can be either "glm", "lm", "glmer",
#'    "lmer", "gam"
#' @param type
#'  Type of Anova
#' @returns
#'  Dataframe with p-values of predictors
get_model_p_values <- function(model,
                               model_type,
                               type = 3) {
  if (model_type %in% c("glm", "lm")) {
    p_values <- stats::anova(model) |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")
  } else if (model_type %in% c("gam")) {
    p_values <- stats::anova(model)$p.table |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")
  } else if (model_type %in% c("lmer", "glmer")) {
    p_values <- car::Anova(model,
                           type = type,
                           test.statistic = "Chisq") |>
      tibble::rownames_to_column(var = "predictor")
  } else { # nls
    p_values <- stats::coef(summary(model))
  }

  p_values
}

#' Get removable terms
#'
#' @param model
#'  A model to determine removable terms for
#' @param model_type
#'  The type of model to get removable terms for. Can be either "glm", "lm",
#'    "glmer", "lmer", "gam"
#' @returns
#'  A dataframe with removable terms and the name of the p-value column
get_removable_terms <- function(model,
                                model_type) {
  if (model_type == "gam") {
    pterms <- summary(model)$pTerms.table |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")
    adjustable_terms <- summary(model)$s.table |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor") |>
      dplyr::bind_rows(pterms)

    p_col <- "p-value"
  } else if (model_type == "nls") {
    adjustable_terms <- stats::coef(summary(model)) |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")

    p_col <- colnames(adjustable_terms)[[grep("Pr\\(",
                                              colnames(adjustable_terms))]]
  } else {
    family <- stats::family(model)$family

    stats_test <- get_stats_test(model_type, family)
    adjustable_terms <- stats::drop1(model, test = stats_test) |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")

    p_col <- colnames(adjustable_terms)[[grep("Pr\\(",
                                              colnames(adjustable_terms))]]
  }
  adjustable_terms <- adjustable_terms |>
    dplyr::filter((.data$predictor != "<none>") &
                    (.data[[p_col]] >= 0.05)) |>
    dplyr::arrange(dplyr::desc(.data[[p_col]]))

  list(adjustable_terms, p_col)
}

#' Get addable terms
#'
#' @param formula
#'  The full formula to get terms to add from
#' @param model
#'  A model to determine removable terms for
#' @param model_type
#'  The type of model to get removable terms for. Can be either "glm", "lm",
#'    "glmer", "lmer", "gam"
#' @returns
#'  A dataframe with terms to add and the name of the p-value column
get_addable_terms <- function(formula,
                              model,
                              model_type) {
  if (model_type == "gam") {
    current_terms <- stats::terms(stats::formula(model),
                                  specials = c("s", "te", "ti", "t2"))
    upper_terms <- stats::terms(stats::as.formula(formula),
                                specials = c("s", "te", "ti", "t2"))
    adjustable_terms <- stats::add.scope(current_terms, upper_terms)
    adjustable_terms <- data.frame(list(predictor = adjustable_terms))

    p_col <- "predictor"
  } else {
    family <- stats::family(model)$family
    stats_test <- get_stats_test(model_type, family)

    adjustable_terms <- tryCatch(
      stats::add1(model,
                  reformulas::nobars(formula),
                  test = stats_test) |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "predictor"),
      error = function(e) {
        if (grepl("no terms in scope", e$message)) {
          tibble::tibble("predictor" = character(),
                         "Pr(" = numeric())
        } else {
          stop(e$message)
        }
      }
    )

    p_col <- colnames(adjustable_terms)[[grep("Pr\\(",
                                              colnames(adjustable_terms))]]
  }

  adjustable_terms <- adjustable_terms |>
    dplyr::filter(.data$predictor != "<none>") |>
    dplyr::arrange(.data[[p_col]])

  list(adjustable_terms, p_col)
}

#' Get appropriate statistical test
#'
#' @param model_type
#'  The type of model to get removable terms for. Can be either "glm", "lm",
#'    "glmer", "lmer", "gam"
#' @param family
#'  The model family. Default: gaussian
#' @returns
#'  The appropriate statistical test to be used with add1/drop1
get_stats_test <- function(model_type,
                           family = stats::gaussian) {
  if ((model_type == "glm") &&
        (family %in% c("binomial", "poisson"))) {
    stats_test = "LRT"
  } else if (model_type %in% c("lmer", "glmer")) {
    stats_test = "Chisq"
  } else {
    stats_test = "F"
  }

  stats_test
}

#' New model is better than old
#'
#' @param assess1
#'  The evaluation results of model 1
#' @param assess2
#'  The evaluation results of model 2
#' @returns
#'  Boolean that describes whether model 1 is better than model 2
new_model_is_better <- function(assess1, assess2) {
  evaluation_results <- c()
  for (eval_m in names(assess1)) {
    stat1 <- assess1[[eval_m]]
    stat2 <- assess2[[eval_m]]
    diff <- stat1 - stat2
    evaluation_results <- append(evaluation_results, diff < 0)
  }

  new_is_better <- length(evaluation_results[evaluation_results])
  total <- length(evaluation_results)
  new_is_better >= (total / 2)
}

#' Remove autocorrelated predictors from formula
#'
#' @param formula
#'  Formula to remove predictors from
#' @param predictors
#'  Predictors to remove
#' @returns
#'  Updated formula without autocorrelated predictors
remove_autocor_predictors <- function(formula,
                                      predictors) {
  terms <- attr(stats::terms.formula(formula), "term.labels")
  interactions <- extract_interactions(terms)

  for (pred in predictors) {
    d <- paste(". ~ . -", pred)
    formula <- stats::update(formula, d)
    terms <- attr(stats::terms.formula(formula), "term.labels")
    interactions <- extract_interactions(terms)
    interacts <- interactions$predictor[interactions$main_effect == pred] |>
      unique()
    for (interact in interacts) {
      d <- paste(". ~ . -", interact)
      formula <- stats::update(formula, d)
    }
  }

  formula
}

#' Extract columns mentioned in formula
#'
#' @param formula
#'  A model formula
#' @param data
#'  Data with columns mentioned in formula
#' @returns
#'  Returns all data columns mentioned in formula
formula_related_cols <- function(formula, data) {
  terms <- attr(stats::terms.formula(formula), "term.labels")
  interactions <- extract_interactions(terms)
  main_effects <- interactions$main_effect
  all_main <- unique(append(main_effects,
                            terms[!(terms %in% interactions$predictor)]))
  numeric_columns <- data |>
    dplyr::select_if(function(col) !(is.character(col) || is.factor(col)))
  all_main[all_main %in% colnames(numeric_columns)]
}

#' Extract categorical variables
#'
#' @param model_frame
#'  The [model.frame()] of a model
#' @returns
#'  A dataframe of categorical variables with their levels
extract_categorical_vars <- function(model_frame) {
  character_columns <- model_frame |>
    dplyr::select_if(function(col) is.character(col) || is.factor(col))

  character_columns <- character_columns |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ paste0(dplyr::cur_column(), as.character(.x))
      )
    )

  lapply(character_columns, function(col) {
    as.data.frame(list(trait = unique(col)))
  }) |>
    dplyr::bind_rows(.id = "var") |>
    rbind(tibble::tibble(trait = character(),
                         var = character()))
}

#' Make level combinations of factors within interaction
#'
#' @param x
#'  A dataframe group
#' @param y
#'  The grouping variable(s)
#' @returns
#'  All combinations of factor levels belonging to an interaction
make_interaction_combinations <- function(x, y) {
  main_effects <- unique(x$main_effect)
  trait_values <- lapply(
    main_effects,
    function(effect) {
      x$predictor[x$main_effect == effect & !is.na(x$predictor)]
    }
  )
  names(trait_values) <- main_effects

  combinations <- expand.grid(
    trait_values,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  combinations$interaction <- do.call(paste, c(combinations, sep = ":"))
  combinations$num_coeffs <- ncol(combinations) - 1

  combinations |>
    tidyr::pivot_longer(-c("interaction", "num_coeffs"),
                        names_to = "term",
                        values_to = "predictor") |>
    dplyr::select(-c("term"))
}

#' Add reference factor levels
#'
#' @param model_overview
#'  Overview of model coefficients with estimates and p-values
#' @param categorical_vars
#'  Dataframe with categorical variables and factor levels
#' @param interactions
#'  Dataframe with interactions and underlying main effects
#' @returns
#'  Updated model overview with information on interactions and categorical
#'    coefficients.
add_reference_factors <- function(model_overview,
                                  categorical_vars,
                                  interactions) {
  ref_interactions <- interactions |>
    dplyr::group_by(.data$term) |>
    dplyr::filter((!all(.data$is_numeric)) &
                    (!.data$interaction %in% model_overview$predictor)) |>
    dplyr::ungroup() |>
    dplyr::select(c("interaction", "term")) |>
    dplyr::rename(predictor = "interaction") |>
    dplyr::mutate(is_ref = TRUE) |>
    dplyr::distinct(.data$predictor, .keep_all = TRUE)

  interactions_sub <- interactions |>
    dplyr::select(c("interaction", "term")) |>
    dplyr::rename(predictor = "interaction") |>
    dplyr::filter(!.data$predictor %in% ref_interactions$predictor) |>
    dplyr::distinct(.data$predictor, .keep_all = TRUE)

  ref_categorical_vars <- categorical_vars |>
    dplyr::select(c("predictor", "term")) |>
    dplyr::mutate(is_ref = TRUE) |>
    dplyr::filter((!.data$predictor %in% model_overview$predictor))

  categorical_vars_sub <- categorical_vars |>
    dplyr::select(c("predictor", "term")) |>
    dplyr::filter(!.data$predictor %in% ref_categorical_vars$predictor)


  if (length(grep("Pr\\(", colnames(model_overview))) == 0) {
    model_overview$p_value <- 1
    p_col <- "p_value"
  } else {
    p_col <- colnames(model_overview)[[
      grep("Pr\\(", colnames(model_overview))
    ]]
  }

  model_overview <- model_overview |>
    dplyr::left_join(interactions_sub, by = "predictor") |>
    dplyr::left_join(categorical_vars_sub, by = "predictor") |>
    dplyr::mutate(term = ifelse(is.na(.data$term.x),
                                ifelse(is.na(.data$term.y),
                                       .data$predictor,
                                       .data$term.y),
                                .data$term.x)) |>
    dplyr::select(-c("term.x", "term.y")) |>
    dplyr::mutate(is_ref = FALSE) |>
    dplyr::bind_rows(ref_interactions) |>
    dplyr::bind_rows(ref_categorical_vars) |>
    dplyr::mutate(is_cat = ifelse(
      (.data$predictor %in% interactions$interaction[interactions$is_numeric]) |
        (.data$predictor %in% categorical_vars$predictor), TRUE, FALSE
    ),
    Estimate = ifelse(is.na(.data$Estimate), 0, .data$Estimate),
    term = ifelse(is.na(.data$term), .data$predictor, .data$term),
    is_interaction = dplyr::case_when(
      .data$predictor %in% interactions$interaction ~ TRUE,
      TRUE ~ FALSE
    )) |>
    dplyr::arrange(dplyr::desc(.data$term), .data$is_ref) |>
    dplyr::mutate(
      formatted_pred = dplyr::case_when(
        !.data$is_ref & .data$is_cat ~ paste0("italic(", .data$predictor, ")"),
        .data$is_ref ~ paste0("bold(", .data$predictor, ")"),
        TRUE ~ .data$predictor
      ),
      formatted_pred = factor(
        .data$formatted_pred,
        levels = .data$formatted_pred
      ),
      significance = dplyr::case_when(
        .data[[p_col]] < 0.001 ~ "***",
        .data[[p_col]] < 0.01 ~ "**",
        .data[[p_col]] < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  formatted_labels <- lapply(model_overview$formatted_pred, function(x) {
    if (grepl("^italic", x) || grepl("^bold", x)) {
      parse(text = as.character(x))[[1]]
    } else {
      as.character(x)
    }
  })

  list(model_overview, formatted_labels)
}

#' Get all factor levels as terms
#'
#' Lists all factor level terms, including in interactions
#' @param model
#'  A model to extract terms from
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm" or "glm"
#' @param model_overview
#'  Overview of model coefficients with estimates and p-values
#' @returns
#'  A dataframe with interactions covering at least one categorical variable
#'    with their factor levels and a dataframe with non-interaction
#'    categorical variables with their factor levels
get_term_factors <- function(model,
                             model_type,
                             model_overview) {
  categorical_vars <- extract_categorical_vars(stats::model.frame(model))
  categorical_vars <- categorical_vars |>
    dplyr::rename(predictor = "trait", term = "var") |>
    dplyr::filter(!.data$term %in% model_overview$predictor) |>
    dplyr::left_join(model_overview, by = "predictor")

  if (model_type %in% c("lmer", "glmer")) {
    formula <- reformulas::nobars(
      stats::formula(model)
    )
  } else {
    formula <- stats::formula(model)
  }
  terms <- attr(stats::terms.formula(formula), "term.labels")
  interactions <- extract_interactions(terms)

  interactions <- interactions |>
    dplyr::rename(term = "predictor") |>
    dplyr::left_join(categorical_vars[, c("term", "predictor")],
                     by = c("main_effect" = "term"),
                     relationship = "many-to-many") |>
    dplyr::mutate(is_numeric = is.na(.data$predictor),
                  predictor = ifelse(is.na(.data$predictor),
                                     .data$main_effect,
                                     .data$predictor))

  if (nrow(interactions) > 0) {
    interactions_w_traits <- interactions |>
      dplyr::group_by(.data$term) |>
      dplyr::group_modify(make_interaction_combinations) |>
      dplyr::ungroup()
    interactions <- interactions_w_traits |>
      dplyr::left_join(interactions, by = c("term", "predictor"))
  } else {
    interactions <- tibble::tibble("term" = character(),
                                   "main_effect" = character(),
                                   "predictor" = character(),
                                   "interaction" = character(),
                                   "is_numeric" = logical())
  }

  numeric_vars <- setdiff(setdiff(terms,
                                  categorical_vars$term),
                          interactions$term)

  list(categorical_vars, interactions, numeric_vars)
}
