library(tidyr)
library(dplyr)
library(tibble)
library(corrplot)
library(stringr)


#' This is data to be included in my package
#'
#' @name dataset_info
#' @docType data
#' @author Lara Kösters \email{lkoesters@@bgc-jena.mpg.de}
#' @keywords data
NULL

#' Autocorrelation removal
#'
#' Remove autocorrelations
#' @param df A (wide) dataframe with predictors and response as columns
#' @param coefficients A list of predictors sorted by their relevance (most to least)
#' @param automatic_removal A boolean to decide whether to automatically remove autocorrelations or leave it to the user to do so manually
#' @param autocorrelation_threshold The threshold deciding (together with the p-values) whether two variables are autocorrelated
#' @param correlation_method The method used for correlation calculation
#' @return An edited list of coefficients where autocorrelated variables are removed (according to sorted list of coefficients)
#' @examples 
#' data("dataset_info");
#' coefficients <- c('val_min_gene_length', 'val_max_gene_length',
#'                   'train_min_gene_length', 'train_max_gene_length',
#'                   'train_num_samples', 'val_num_samples',
#'                   'mean_gene_length', 'val_mean_gene_length',
#'                   'train_mean_gene_length');
#' remove_autocorrelations(df, coefficients, automatic_removal = TRUE);
#' @export
remove_autocorrelations <- function(df, coefficients, automatic_removal=TRUE,
                                    autocorrelation_threshold=0.8,
                                    correlation_method='pearson') {
  correlations <- as.data.frame(cor(df[, coefficients], method=correlation_method))
  # compute p-values
  correlation_p_values <- as.data.frame(cor.mtest(correlations)$p)
  
  correlations_l <- correlations %>%
    rownames_to_column(var='coefficientA') %>%
    pivot_longer(!coefficientA, names_to = "coefficientB", values_to = "correlation")
  
  correlations_p_values_l <- correlation_p_values %>%
    rownames_to_column(var='coefficientA') %>%
    pivot_longer(!coefficientA, names_to = "coefficientB", values_to = "p_value")
  
  correlations_complete <- merge(correlations_l, correlations_p_values_l, by=c('coefficientA', 'coefficientB'))
  
  correlations_complete <- correlations_complete %>%
    filter(coefficientA != coefficientB) %>%
    mutate(col1 = pmin(coefficientA, coefficientB), col2 = pmax(coefficientA, coefficientB), 
           comparison = paste(col1, col2)) %>%
    distinct(comparison, .keep_all = TRUE) %>%
    select(!c(col1, col2, comparison))
  
  autocorrelations <- correlations_complete[(correlations_complete$cor >= autocorrelation_threshold) & (correlations_complete$p_value < 0.05),]
  
  if (nrow(autocorrelations) > 0) {
    if (automatic_removal) {
      coefficients_df <- data.frame(idx = 1:length(coefficients),
                                    coefficient = coefficients)
      
      autocorrelations <- merge(autocorrelations, coefficients_df[c('coefficient', 'idx')], by.x='coefficientA', by.y='coefficient')
      autocorrelations <- merge(autocorrelations, coefficients_df[c('coefficient', 'idx')], by.x='coefficientB', by.y='coefficient', suffixes = c('1', '2'))
      autocorrelations <- autocorrelations %>%
        rowwise() %>%
        mutate(idx_smaller = sort(c(idx1, idx2))[[1]],
               idx_bigger = sort(c(idx1, idx2))[[2]],
               comb = paste(idx_smaller, idx_bigger, sep='|'))
      
      for (i in 1:nrow(autocorrelations)) {
        # check if we have the following constellation: A=B, B=C, A!=C
        autocor_row <- autocorrelations[order(autocorrelations$idx_bigger, decreasing = TRUE),][i,]
        c <- autocor_row$idx_bigger
        b <- autocor_row$idx_smaller
        
        b_to_a <- autocorrelations[autocorrelations$idx_bigger == b,]
        if (nrow(b_to_a) > 0) {
          for (a in b_to_a$idx_smaller) {
            if(nrow(autocorrelatins[(autocorrelations$idx_bigger == autocor_row$idx_bigger) & (autocorrelations$idx_smaller == a),]) == 0) {
              # A != C but A==B and B==C: remove B
              coefficients <- coefficients[-b]
            } else {
              coefficients <- coefficients[-c]
            }
          }
        } else {
          coefficients <- coefficients[-c]
        }
      }
    } else {
      f <- str_interp("autocorrelations_${format(Sys.time(), '%m%d%Y_%H%M')}.tsv")
      write.table(autocorrelations, f, 
                  sep = '\t', quote = FALSE, row.names = FALSE)
      stop(str_interp("Some of your variables are autocorrelated. Check ${f} for more info"))
    }
  }
  
  return(coefficients)
}

#' Model simplification
#'
#' Simplify models
#' @param glm_model Starting glm model that will be simplified
#' @param aic_bic Whether to apply AIC and BIC difference in addition to ANOVA to determine if model is better/worse than ancestor
#' @return A log with models + statistics as well as final model
#' @examples
#' data("dataset_info");
#' coefficients <- c('val_min_gene_length', 'val_max_gene_length',
#'                   'train_min_gene_length', 'train_max_gene_length',
#'                   'train_num_samples', 'val_num_samples',
#'                   'mean_gene_length', 'val_mean_gene_length',
#'                   'train_mean_gene_length');
#' frm <- as.formula(paste("species_confusion", "~", paste(coefficients, collapse='+')))
#' glm_model <- eval(bquote(glm(as.formula(frm), data=.(dataset_info), family = .(glm_family))))
#' simplify_model(glm_model, aic_bic=TRUE);
#' @export
simplify_model <- function(glm_model, aic_bic=FALSE, model_method='binomial') {
  # setup
  simplification_history <- list(list("summary" = summary(glm_model),
                                      "model" = glm_model,
                                      "final" = FALSE))
  final_model <- NA
  
  p_values_model <- data.frame(p_values = coef(summary(glm_model))[,4]) %>%
    rownames_to_column(var='coefficients') %>%
    arrange(desc(p_values)) %>%
    filter(coefficients != '(Intercept)')

  if (p_values_model[1, 1] < 0.1) {
    # All p-values are below p<0.1! Stop looping and accept the last model
    marginal_sign_vars <- p_values_model[(p_values_model$p_values < 0.1) & (p_values_model$p_values >= 0.05),]
    sign_vars <- p_values_model[p_values_model$p_values < 0.05,]
    simplification_history[[1]]$anovas <- res_anova
    simplification_history[[1]]$significant_variables <- sign_vars$coefficients
    simplification_history[[1]]$marginally_significant_variables <- marginal_sign_vars$coefficients
    final_model <- simplification_history[[1]]
  } else {
    
    if (!grepl('quasi', model_method) & aic_bic) {
      aic_res <- AIC(glm_model)
      bic_res <- BIC(glm_model)
      simplification_history[[1]]$aic_res <- aic_res
      simplification_history[[1]]$bic_res <- bic_res
    }
    
    for (i in 2:length(p_values_model$coefficients)) {
      last_model <- simplification_history[[i-1]]$model
      glm_model <- update(last_model, paste("~ . -", p_values_model[1, 1]))
      p_values_model <- data.frame(p_values = coef(summary(glm_model))[,4]) %>%
        rownames_to_column(var='coefficients') %>%
        arrange(desc(p_values)) %>%
        filter(coefficients != '(Intercept)')
      model_summary <- summary(glm_model)
      simplification_history[[i]] <- list("summary" = model_summary,
                                        "model" = glm_model)

      anova_res <- anova(glm_model, last_model)
      simplification_history[[i]]$anova_res <- anova_res
      no_more_simplification <- anova_res$Deviance[2] > 0
      
      if (!grepl('quasi', model_method) & aic_bic) {
        aic_res <- AIC(glm_model)
        bic_res <- BIC(glm_model)
        simplification_history[[i]]$aic_res <- aic_res
        simplification_history[[i]]$bic_res <- bic_res
        
        aic_no_more_simplification <- (simplification_history[[i-1]]$aic_res - aic_res) < 0
        bic_no_more_simplification <- (simplification_history[[i-1]]$bic_res - bic_res) < 0
        no_more_simplification <- all(no_more_simplification, aic_no_more_simplification, bic_no_more_simplification)
      }

      if (no_more_simplification > 0 & p_values_model[1, 2] < 0.1) {
        marginal_sign_vars <- p_values_model[(p_values_model$p_values < 0.1) & (p_values_model$p_values >= 0.05),]
        sign_vars <- p_values_model[p_values_model$p_values < 0.05,]
        
        simplification_history[[i]]$significant_variables <- sign_vars$coefficients
        simplification_history[[i]]$marginally_significant_variables <- marginal_sign_vars$coefficients
        final_model <- simplification_history[[i]]
        break
      }
    }
  }

  # TODO: potential memory problem because data gets stored in glm_model object
  return(list(history = simplification_history, final_model = final_model))
}


#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables that do not signicantly predict response variable
#' @param df A (wide) dataframe with predictors and response as columns
#' @param response_cols A list of response columns sorted by their relevance (most to least)
#' @param predictor_cols A list of predictors sorted by their relevance (most to least)
#' @param glm_family The family used for glm calculation
#' @param autocorrelation_threshold The threshold deciding (together with the p-values) whether two variables are autocorrelated
#' @param automatic_removal Whether to automatically remove autocorrelations
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param plot Whether to plot correlations for predictors with significant impact
#' @param aic_bic Whether to apply AIC and BIC difference in addition to ANOVA to determine if model is better/worse than ancestor
#' @param autocorrelation_method The method used for correlation calculation
#' @return An edited list of coefficients where autocorrelated variables are removed (according to sorted list of coefficients)
#' @examples 
#' data("dataset_info");
#' coefficients <- c('val_min_gene_length', 'val_max_gene_length',
#'                   'train_min_gene_length', 'train_max_gene_length',
#'                   'train_num_samples', 'val_num_samples',
#'                   'mean_gene_length', 'val_mean_gene_length',
#'                   'train_mean_gene_length');
#' optimize_model(dataset_info, "species_confusion", coefficients, automatic_removal = TRUE);
#' @export
optimize_model <- function(df, response_cols, predictor_cols, glm_family='binomial', 
                           autocorrelation_threshold=0.8, automatic_removal=FALSE, 
                           round_p=5, plot=FALSE, aic_bic=FALSE, autocorrelation_method='pearson') {
  model_histories <- list()
  final_model_overviews <- list()
  
  for (i1 in 1:length(response_cols)) {
    response_col <- response_cols[[i1]]
      
    frm <- as.formula(paste(response_col, "~", paste(predictor_cols, collapse='+')))
    glm_model <- eval(bquote(glm(as.formula(frm), data=.(df), family = .(glm_family))))
    predictor_cols <- remove_autocorrelations(df, predictor_cols,
                                         automatic_removal = automatic_removal,
                                         autocorrelation_threshold = autocorrelation_threshold,
                                         correlation_method = correlation_method)
    simplified_model_info <- simplify_model(glm_model, aic_bic=aic_bic)
    model_histories <- simplified_model_info$history
    final_model <- simplified_model_info$final_model
    
    all_vars <- c(final_model$significant_variables, final_model$marginally_significant_variables)
    
    if (length(all_vars) != 0) {
      estimates <- as.data.frame(coef(summary(final_model$model)))
      estimates$predictor <- rownames(estimates)
      estimates$response <- response_col
      estimates <- estimates[estimates$predictor != '(Intercept)',]
      
      for (i in 1:length(all_vars)) {
        sign_var <- all_vars[[i]]
        if (plot) {
          plot_vars(df, sign_var, response_col, final_model$model, all_vars)
        }
        final_model_overview <- data.frame('response' = c(response_col),
                                           'predictor' = c(sign_var))
        
        final_model_overview <- merge(final_model_overview, estimates, 
                                      by=c('response', 'predictor'), 
                                      all.x = TRUE)
        final_model_overviews[[i + ((i1-1)*length(response_cols))]] <- final_model_overview
      }
      # im final model overview ist kein R Wert... wollen wir den?
    }
  }
  
  models_overview <- do.call(rbind, final_model_overviews)
  
  p_value_col <- tail(colnames(models_overview)[-1], n=1)
  
  models_overview <- models_overview %>%
    mutate_if(is.numeric, round, digits=round_p) %>%
    mutate(significance = case_when(
      !!ensym(p_value_col) < 0.0001 ~ '****',
      !!ensym(p_value_col) < 0.001 ~ '***',
      !!ensym(p_value_col) < 0.01 ~ '**',
      !!ensym(p_value_col) < 0.05 ~ '*',
      !!ensym(p_value_col) < 0.1 ~ '',
      TRUE ~ ''
    ))
  
  return(list(
    optimization_history = model_histories,
    summary = models_overview
  ))
}
