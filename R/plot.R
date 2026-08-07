#' Plot post model selection inference result
#'
#' @param psi_info
#'  Dataframe with raw and corrected p-values of coefficients
#' @param p_threshold
#'  P-value threshold for significance evaluation
#' @param label_size
#'  Size of labels within post-selection inference plot.
#'  Default: 2.5
#' @returns
#'  Plot with raw and corrected p-values of coefficients
plot_psi <- function(psi_info,
                     p_threshold,
                     label_size = 2.5) {
  psi_info <- psi_info |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.numeric),
        ~ dplyr::case_when(
          is.na(.x) ~ 1,
          .x == 0 ~ 1e-12,
          TRUE ~ .x
        )
      ),
      raw_sig = ifelse(.data$raw_p_value <= p_threshold, "Yes", "No"),
      corr_sig = ifelse(.data$corrected_p_value <= p_threshold, "Yes", "No")
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = psi_info,
      ggplot2::aes(
        x = 1,
        xend = 2,
        y = log10(.data$raw_p_value),
        yend = log10(.data$corrected_p_value),
        color = abs(log10(.data$corrected_p_value) - log10(.data$raw_p_value))
      ),
      linewidth = 0.8, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = psi_info,
      ggplot2::aes(x = 1,
                   y = log10(.data$raw_p_value),
                   fill = .data$raw_sig),
      shape = 21, color = "black", size = 3
    ) +
    ggplot2::geom_point(
      data = psi_info,
      ggplot2::aes(x = 2,
                   y = log10(.data$corrected_p_value),
                   fill = .data$corr_sig),
      shape = 21, color = "black", size = 3
    ) +
    ggrepel::geom_text_repel(
      data = psi_info,
      ggplot2::aes(x = 1,
                   y = log10(.data$raw_p_value),
                   label = .data$predictor),
      size = label_size, hjust = 0, nudge_x = 0.1, segment.size = 0.25,
      box.padding = 0.3, direction = "y", max.overlaps = Inf
    ) +
    ggplot2::geom_hline(
      yintercept = log10(p_threshold), linetype = "dashed", color = "gray40"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2), labels = c("Raw", "Corrected")
    ) +
    ggplot2::scale_color_gradient(
      low = "gray70", high = "purple", name = "|\u0394 log10(p)|"
    ) +
    ggplot2::scale_fill_manual(
      values = c("Yes" = "#1b9e77", "No" = "#d95f02"), name = "Significant"
    ) +
    ggplot2::labs(
      x = "p-value", y = expression(log[10](italic(p))),
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

#' Plot model
#'
#' Plots model features such as estimates, effect sizes, and plots
#'  visualizing the relationship between the response variable and predictors
#' @param model
#'  A model for which to plot features
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "gam"
#' @param quality_assessment
#'  The mode of model quality assessment. Either "baseR" or "performance".
#'    Default: "baseR"
#' @param plot_relationships
#'  Whether to plot regression, effect size, and estimates.
#'  Default: TRUE
#' @param test
#'  Either "t.test" or "wilcox".
#'  Used to calculate statistics for regression plots of categorical variables.
#'  Default: "wilcox"
#' @param plot_type
#'  Either "boxplot" or "violin".
#'  Used to plot regression plots for categorical variables.
#'  Default: "boxplot"
#' @param plot_curve
#'  Whether to plot [ggplot2::geom_smooth()] in regression plots. Default: TRUE
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param point_position
#'  Position adjustment for model feature plots. See the position paramter of
#'    [ggplot2::geom_point()] for more information
#' @return
#'  A list of plots (estimate, regression, and effect size)
#'  and statistics for categorical variables.
#' @examples
#'  # setup
#'  data("plants")
#'  # generate an example glm model for plotting
#'  final_model <- create_model(
#'   sexual_seed_prop ~ altitude +
#'     solar_radiation +
#'     annual_mean_temperature +
#'     isothermality +
#'     habitat +
#'     ploidy +
#'     solar_radiation:isothermality,
#'   data = plants,
#'   model_type = "glm",
#'   family = "quasibinomial"
#'  )
#'  # plot the model's features including plotting of estimates,
#'  # numeric and categorical variables using boxplots where possible
#'  p <- plot_model(
#'   final_model,
#'   model_type = "glm",
#'   plot_type = "boxplot",
#'   point_position = "jitter",
#'   test = "wilcox",
#'   round_p = 3
#'  )
#' @export
plot_model <- function(model,
                       model_type,
                       quality_assessment = "baseR",
                       plot_relationships = TRUE,
                       test = "wilcox",
                       plot_type = "boxplot",
                       plot_curve = TRUE,
                       round_p = 5,
                       point_position = "jitter") {
  model_plots <- list()

  model_plots$quality_check <- assess_basic_model_quality(
    model,
    quality_assessment,
    model_type
  )

  if (model_type == "gam") {
    p_table <- summary(model)$p.table
    model_overview <- p_table |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")
  } else {
    model_overview <- stats::coef(summary(model)) |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "predictor")
  }

  m_matrix <- stats::model.matrix(model)
  m_frame <- stats::model.frame(model)
  if (model_type %in% c("lmer", "glmer")) {
    formula <- reformulas::nobars(
      stats::formula(model)
    )
  } else {
    formula <- model$formula
  }
  response_str <- deparse1(formula.tools::lhs(formula))
  c(categorical_vars,
    interactions,
    numeric_vars) %<-% prepare_plot_data(model,
                                         model_type,
                                         model_overview)

  if (model_type == "gam") {
    s_table <- summary(model)$s.table
    numeric_vars <- numeric_vars[!numeric_vars %in% rownames(s_table)]
  }

  c(model_overview,
    formatted_labels) %<-% add_reference_factors(model_overview,
                                                 categorical_vars,
                                                 interactions)

  if (plot_relationships) {
    model_plots$estimates <- plot_estimates(model_overview,
                                            formatted_labels)
    model_plots$effect_sizes <- plot_effect_sizes(model_overview,
                                                  formatted_labels)
    categorical_plots <- plot_categorical_vars(categorical_vars,
                                               response_str,
                                               m_frame,
                                               model_overview,
                                               test,
                                               plot_type)
    family <- stats::family(model)$family
    numeric_plots <- plot_numeric_vars(numeric_vars,
                                       response_str,
                                       m_frame,
                                       m_matrix,
                                       model_overview,
                                       model_type,
                                       family,
                                       plot_curve,
                                       round_p,
                                       point_position)
    interaction_plots <- plot_interactions(interactions,
                                           response_str,
                                           m_frame,
                                           m_matrix,
                                           model_overview,
                                           model_type,
                                           family,
                                           plot_curve,
                                           round_p,
                                           point_position)
    model_plots$categorical_variables <- categorical_plots
    model_plots$relationships <- c(numeric_plots,
                                   interaction_plots)
  }

  model_plots
}

#' Prepare data for plotting
#'
#' @param model
#'  A model for which to plot features
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "gam"
#' @param model_overview
#'  An overview of model coefficients with estimates and p-values
#' @returns
#'  Categorical variables, interactions, and numeric variables
prepare_plot_data <- function(model,
                              model_type,
                              model_overview) {
  model_overview <- model_overview |>
    dplyr::mutate(
      effect_direction = dplyr::case_when(
        .data$Estimate < 0 ~ "negative",
        TRUE ~ "positive"
      )
    )

  get_term_factors(model, model_type, model_overview)
}

#' Assess basic model quality
#'
#' @param model
#'  A model for which to plot features
#' @param quality_assessment
#'  The mode of model quality assessment. Either "baseR" or "performance"
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "gam"
#' @returns
#'  Either a grid with plots or NA if quality_assessment differed from "baseR"
#'    and "performance"
assess_basic_model_quality <- function(
  model,
  quality_assessment,
  model_type
) {
  if (quality_assessment == "performance") {
    if (model_type == "gam") {
      checked_model <- withr::with_pdf(NULL, {
        plot(performance::check_model(
          model,
          type = "normal"
        ))
      })
    } else if (stats::family(model)$family == "binomial") {
      withCallingHandlers(
        checked_model <- withr::with_pdf(NULL, {
          plot(performance::check_model(
            model,
            type = "discrete_both"
          ))
        }),
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
        checked_model <- withr::with_pdf(NULL, {
          plot(performance::check_model(model))
        }),
        message = function(w) {
          if (grepl("unknown labels", w$message)) {
            tryInvokeRestart("muffleMessage")
          }
        }
      )
      checked_model[[3]] <- checked_model[[3]] +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }
  } else if (quality_assessment == "baseR") {
    checked_model <- withr::with_pdf(NULL, {
      graphics::par(mfrow = c(2, 2), mar = rep(2, 4))
      graphics::plot(model)
      grDevices::recordPlot()
    })
  } else {
    checked_model <- NA
  }

  checked_model
}

#' Plot estimates
#'
#' @param model_overview
#'  Overview of model coefficients with estimates
#' @param formatted_labels
#'  A list with formatted labels (italic and bold for factor levels)
#' @returns
#'  An estimate plot
plot_estimates <- function(model_overview,
                           formatted_labels) {
  p <- ggplot2::ggplot(
    model_overview,
    ggplot2::aes(x = .data$Estimate, y = .data$formatted_pred)
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

#' Plot effect sizes
#'
#' @param model_overview
#'  Overview of model coefficients with estimates
#' @param formatted_labels
#'  A list with formatted labels (italic and bold for factor levels)
#' @returns
#'  An effect size plot
plot_effect_sizes <- function(model_overview,
                              formatted_labels) {
  model_overview <- model_overview |>
    dplyr::mutate(
      Estimate_abs = abs(.data$Estimate),
      Est_sum = sum(.data$Estimate_abs),
      `Effect size` = dplyr::case_when(
        .data$Estimate == 0 ~ 0,
        .data$Estimate > 0 ~ (.data$Estimate_abs / .data$Est_sum) * 100,
        .data$Estimate < 0 ~ (.data$Estimate_abs / .data$Est_sum) * -100
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$`Effect size`))

  p <- ggplot2::ggplot(
    model_overview,
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
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_x_discrete(labels = formatted_labels)
  p
}

#' Plot categorical variables
#'
#' @param categorical_vars
#'  A dataframe with categorical variables and corresponding factor levels
#' @param response
#'  A string representing the response variable
#' @param m_frame
#'  A model frame. See [stats::model.frame()]
#' @param model_overview
#'  An overview of model coefficients withe estimates and p-values
#' @param test
#'  Either "t.test" or "wilcox".
#'  Used to calculate statistics for regression plots of categorical variables.
#'  Default: "wilcox"
#' @param plot_type
#'  Either "boxplot" or "violin".
#'  Used to plot regression plots for categorical variables.
#'  Default: "boxplot"
#' @returns
#'  Plots of categorical variables alongside results of statistical tests
plot_categorical_vars <- function(categorical_vars,
                                  response,
                                  m_frame,
                                  model_overview,
                                  test = "wilcox",
                                  plot_type = "boxplot") {
  stat_results <- list()
  plots <- list()
  for (cat_var in unique(categorical_vars$term)) {
    c(m_frame_w_letters,
      stat_result) %<-% run_stats(m_frame,
                                  response,
                                  cat_var,
                                  test = test)
    stat_results[[cat_var]] <- stat_result

    model_overview <- model_overview |>
      dplyr::filter((!.data$is_interaction) & (.data$is_cat)) |>
      dplyr::mutate(
        trait = stringr::str_remove(.data$predictor, paste0("^", .data$term))
      )

    m_frame_w_letters <- m_frame_w_letters |>
      dplyr::left_join(model_overview[, c("trait", "is_ref", "significance")],
                       by = stats::setNames(
                         object = "trait",
                         nm = cat_var
                       )) |>
      dplyr::mutate(significance = dplyr::case_when(
        (.data$significance == "ns") & (.data$is_ref) ~ "ref",
        TRUE ~ .data$significance
      ))

    y_lim_min <- min(m_frame_w_letters[!is.na(m_frame_w_letters[[response]]),
                                       response])
    y_lim_max <- max(m_frame_w_letters[!is.na(m_frame_w_letters[[response]]),
                                       response])

    p <- ggplot2::ggplot(
      data = m_frame_w_letters,
      ggplot2::aes(x = .data[[cat_var]], y = .data[[response]])
    ) +
      ggplot2::scale_color_viridis_d(option = "G", end = 0.9) +
      ggplot2::scale_fill_viridis_d(option = "G", end = 0.9) +
      ggplot2::scale_y_continuous(
        limits = c(y_lim_min, y_lim_max * 1.15),
        expand = ggplot2::expansion(mult = c(0, .05))
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data[[cat_var]], label = .data$letter),
        y = y_lim_max * 1.15,
        check_overlap = TRUE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data[[cat_var]], label = .data$significance),
        y = y_lim_max * 1.05,
        check_overlap = TRUE
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0)
      ) +
      ggplot2::guides(fill = "none", color = "none") +
      ggplot2::labs(
        title = stringr::str_interp(
          "Regression of coefficient: ${cat_var}"
        )
      )

    if (plot_type == "boxplot") {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(color = .data[[cat_var]], fill = .data[[cat_var]])
      )
    } else {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(color = .data[[cat_var]], fill = .data[[cat_var]])
      )
    }

    plots[[cat_var]] <- p
  }

  list(
    plots = plots,
    stat_results = stat_results
  )
}

#' Plot numeric variables
#'
#' @param numeric_vars
#'  A vector with names of numeric variables
#' @param response
#'  A string representing the response variable
#' @param m_frame
#'  A model frame. See [stats::model.frame()]
#' @param m_matrix
#'  A model frame. See [stats::model.matrix()]
#' @param model_overview
#'  An overview of model coefficients withe estimates and p-values
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "gam"
#' @param family
#'  The model family
#' @param plot_curve
#'  Whether to plot [ggplot2::geom_smooth()] in regression plots. Default: TRUE
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param point_position
#'  Position adjustment for model feature plots. See the position paramter of
#'    [ggplot2::geom_point()] for more information. Default: "jitter"
#' @param group_var
#'  A grouping variable for coloring
#' @returns
#'  Plots of numeric variables
plot_numeric_vars <- function(numeric_vars,
                              response,
                              m_frame,
                              m_matrix,
                              model_overview,
                              model_type,
                              family,
                              plot_curve = TRUE,
                              round_p = 5,
                              point_position = "jitter",
                              group_var = NA) {
  plots <- list()
  is_interaction <- typeof(group_var) == "logical"
  if (is_interaction) {
    group_var <- "3644fc2ee6434fc78b66cfedea549a36"
    m_frame[group_var] <- "trait"
  }

  for (numeric_var in numeric_vars) {
    m_frame[numeric_var] <- m_matrix[, numeric_var]

    significance <- model_overview[
      model_overview["predictor"] == numeric_var,
      "significance"
    ]
    if ("Estimate" %in% colnames(model_overview)) {
      estimate <- round(
        model_overview[
          model_overview["predictor"] == numeric_var,
          "Estimate"
        ],
        digits = round_p
      )
    } else {
      estimate = "NA"
    }
    est_len <- length(as.character(estimate))
    x_max <- max(m_frame[[numeric_var]])
    x_min <- min(m_frame[[numeric_var]])
    if (estimate > 0 || is.na(estimate)) {
      box_pos <- x_min + (abs((abs(x_max) - abs(x_min))) * (est_len * .2))
    } else {
      box_pos <- x_max - (abs((abs(x_max) - abs(x_min))) * (est_len * .2))
    }
    plot_label <- stringr::str_interp("estimate = ${estimate}\n${significance}")

    if (is_interaction) {
      plot_title <- stringr::str_interp(
        "Regression of coefficient: ${numeric_var}"
      )
    } else {
      plot_title <- stringr::str_interp(
        "Regression of coefficient: ${numeric_var}:${group_var}"
      )
    }

    p <- ggplot2::ggplot(data = m_frame,
                         ggplot2::aes(x = .data[[numeric_var]],
                                      y = .data[[response]])) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[group_var]],
                                       fill = .data[[group_var]]),
                          pch = 21, position = point_position) +
      ggplot2::scale_color_viridis_d(option = "G", end = 0.9) +
      ggplot2::scale_fill_viridis_d(option = "G", end = 0.9) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = .5)
      ) +
      ggplot2::annotate("label",
                        x = box_pos, y = Inf,
                        label = plot_label,
                        vjust = 1, hjust = .5) +
      ggplot2::labs(title = plot_title)

    if (is_interaction) {
      p <- p +
        ggplot2::guides(fill = "none", color = "none")
    }

    if (plot_curve && significance != "ns") {
      p <- p +
        ggplot2::geom_smooth(
          method = model_type,
          method.args = list(family = family),
          formula = y ~ x,
          se = FALSE
        )
    }

    plots[[numeric_var]] <- p
  }

  plots
}

#' Plot interactions
#'
#' @param interactions
#'  Dataframe with interactions
#' @param response
#'  A string representing the response variable
#' @param m_frame
#'  A model frame. See [stats::model.frame()]
#' @param m_matrix
#'  A model frame. See [stats::model.matrix()]
#' @param model_overview
#'  An overview of model coefficients withe estimates and p-values
#' @param model_type
#'  Model type to be used as character string.
#'  Options: "lm", "glm", "gam"
#' @param family
#'  The model family
#' @param plot_curve
#'  Whether to plot [ggplot2::geom_smooth()] in regression plots. Default: TRUE
#' @param round_p
#'  Convenience parameter for automatic rounding of p-values. Default: 5
#' @param point_position
#'  Position adjustment for model feature plots. See the position paramter of
#'    [ggplot2::geom_point()] for more information. Default: "jitter"
#' @returns
#'  Plots of interactions between numeric variables or between a numeric and a
#'    categorical variable
plot_interactions <- function(interactions,
                              response,
                              m_frame,
                              m_matrix,
                              model_overview,
                              model_type,
                              family,
                              plot_curve = TRUE,
                              round_p = 5,
                              point_position = "jitter") {
  interactions <- check_plot_interactions(interactions) |>
    dplyr::group_by(.data$term) |>
    dplyr::mutate(only_numeric = all(.data$is_numeric))

  numeric_vars <- unique(interactions$term[interactions$only_numeric])
  numeric_plots <- plot_numeric_vars(numeric_vars,
                                     response,
                                     m_frame,
                                     m_matrix,
                                     model_overview,
                                     model_type,
                                     family,
                                     plot_curve,
                                     round_p,
                                     point_position)

  mixed_vars <- interactions[!interactions$only_numeric, ]
  mixed_plots <- list()
  for (mixed_numeric_term in unique(mixed_vars$term)) {
    numeric_var <- unique(mixed_vars$main_effect[
      (mixed_vars$term == mixed_numeric_term) &
        (mixed_vars$is_numeric)
    ])
    cat_var <- unique(mixed_vars$main_effect[
      (mixed_vars$term == mixed_numeric_term) &
        (!mixed_vars$is_numeric)
    ])

    p <- plot_numeric_vars(c(numeric_var),
                           response,
                           m_frame,
                           m_matrix,
                           model_overview,
                           model_type,
                           family,
                           FALSE,
                           round_p,
                           point_position,
                           group_var = cat_var)
    mixed_plots[[mixed_numeric_term]] <- p[[numeric_var]]
  }

  c(numeric_plots, mixed_plots)
}
