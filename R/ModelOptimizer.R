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
                                    correlation_method='spearman') {
  
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
  
  print(correlations_complete)
  
  autocorrelations <- correlations_complete[((correlations_complete$cor >= autocorrelation_threshold) | 
                                               (correlations_complete$cor <= -autocorrelation_threshold)) & 
                                              (correlations_complete$p_value < 0.05),]
  autocorrelations <- autocorrelations %>%
    add_column(note = NA)
  
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
            if(nrow(autocorrelations[(autocorrelations$idx_bigger == autocor_row$idx_bigger) & (autocorrelations$idx_smaller == a),]) == 0) {
              # A != C but A==B and B==C: remove B
              coefficient_b <- coefficients_df[b, 'coefficient']
              coefficients <- coefficients[!coefficients == coefficient_b]
              autocorrelations[i, 'note'] <- str_interp('removed ${coefficient_b}')
            } else {
              coefficient_c <- coefficients_df[c, 'coefficient']
              coefficients <- coefficients[!coefficients == coefficient_c]
              autocorrelations[i, 'note'] <- str_interp('removed ${coefficient_c}')
            }
          }
        } else {
          coefficient_c <- coefficients_df[c, 'coefficient']
          coefficients <- coefficients[!coefficients == coefficient_c]
          autocorrelations[i, 'note'] <- str_interp('removed ${coefficient_c}')
        }
      }
    } else {
      f <- str_interp("autocorrelations_${format(Sys.time(), '%m%d%Y_%H%M')}.tsv")
      write.table(autocorrelations, f, 
                  sep = '\t', quote = FALSE, row.names = FALSE)
      stop(str_interp("Some of your variables are autocorrelated. Check ${f} for more info"))
    }
  }
  
  return(list("predictor_cols" = coefficients, "autocorrelations" = autocorrelations[,c('coefficientA', 'coefficientB', 'correlation', 'p_value', 'note')]))
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
simplify_model <- function(glm_model, aic_bic=FALSE, model_method='binomial',
                           model_summaries=FALSE) {
  # setup
  simplification_history <- list(list("model" = glm_model,
                                      "final" = FALSE))
  
  if (model_summaries) { simplification_history[[1]]$summary = summary(glm_model) }
  
  final_model <- NA
  
  p_values_model <- data.frame(p_values = coef(summary(glm_model))[,4]) %>%
    rownames_to_column(var='coefficients') %>%
    arrange(desc(p_values)) %>%
    filter(coefficients != '(Intercept)')
  
  if (p_values_model[1, 1] < 0.1) {
    # All p-values are below p<0.1! Stop looping and accept the last model
    marginal_sign_vars <- p_values_model[(p_values_model$p_values < 0.1) & (p_values_model$p_values >= 0.05),]
    sign_vars <- p_values_model[p_values_model$p_values < 0.05,]
    
    simplification_history[[i]]$anovas <- res_anova
    simplification_history[[i]]$significant_variables <- sign_vars$coefficients
    simplification_history[[i]]$marginally_significant_variables <- marginal_sign_vars$coefficients
    
    final_model <- simplification_history[[1]]
  } else {
    
    if (!grepl('quasi', model_method) & aic_bic) {
      aic_res <- AIC(glm_model)
      bic_res <- BIC(glm_model)
      simplification_history[[1]][c('aic_res', 'bic_res')] <- c(aic_res, bic_res)
    }
    
    for (i in 2:length(p_values_model$coefficients)) {
      last_model <- simplification_history[[i-1]]$model
      
      glm_model <- update(last_model, paste("~ . -", p_values_model[1, 1]))
      p_values_model <- data.frame(p_values = coef(summary(glm_model))[,4]) %>%
        rownames_to_column(var='coefficients') %>%
        arrange(desc(p_values)) %>%
        filter(coefficients != '(Intercept)')
      model_summary <- summary(glm_model)
      simplification_history[[i]] <- list("model" = glm_model)
      if (model_summaries) { simplification_history[[i]]$summary = summary(glm_model) }
      
      anova_res <- anova(glm_model, last_model)
      simplification_history[[i]]$anova_res <- anova_res
      no_more_simplification <- anova_res$Deviance[2] > 0
      
      if (!grepl('quasi', model_method) & aic_bic) {
        aic_res <- AIC(glm_model)
        bic_res <- BIC(glm_model)
        simplification_history[[i]][c('aic_res', 'bic_res')] <- c(aic_res, bic_res)
        
        aic_no_more_simplification <- (simplification_history[[i-1]]$aic_res - aic_res) < 0
        bic_no_more_simplification <- (simplification_history[[i-1]]$bic_res - bic_res) < 0
        no_more_simplification <- all(no_more_simplification, aic_no_more_simplification, bic_no_more_simplification)
      }
      
      if (no_more_simplification > 0 & p_values_model[1, 2] < 0.1) {
        print(coef(summary(glm_model)))
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
#' @param resp_preds A dictionary of response columns and respective predictors sorted by their relevance (most to least)
#' @param glm_family The family used for glm calculation
#' @param autocorrelation_threshold The threshold deciding (together with the p-values) whether two variables are autocorrelated
#' @param automatic_removal Whether to automatically remove autocorrelations
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param plot Whether to plot correlations for predictors with significant impact
#' @param aic_bic Whether to apply AIC and BIC difference in addition to ANOVA to determine if model is better/worse than ancestor
#' @param correlation_method The method used for correlation calculation
#' @param glm_iter The number of max iterations used for the glm calculation
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

# resp_preds <- c('response' = c('predictor_col1', 'predictor_col2'))
# 
# resp_preds <- c('response' = c(c('predictor_col1' = c('^2', 'log')), 'predictor_col2'))
# 
# transforms = c('log_10', 'log_e', 'arcsin', '')
# 
# transform_response <- 'nutzer definierte Liste mit transforms für die response var'
# 
# (x1 + x2 + x3)^2 <- 'dann rechnet glm das automatisch'
# quadratic <- 'Liste '


optimize_model <- function(df, resp_preds, glm_family=c('binomial'), 
                           autocorrelation_threshold=0.8, automatic_removal=FALSE, 
                           round_p=5, plot=FALSE, aic_bic=FALSE, correlation_method='spearman',
                           model_summaries=FALSE, glm_iter=1000, transform_response=FALSE,
                           interactions=FALSE, quadratic) {
  model_histories <- list()
  final_model_overviews <- list()
  autocorrelations_overview <- list()
  
  for (i1 in 1:length(resp_preds)) {
    response_col <- names(resp_preds)[[i1]]
    predictor_cols <- resp_preds[[response_col]]
    df_sub <- df[!is.na(df[[response_col]]),]
    print(response_col)
    process_code <- str_interp('${response_col}_${i1}')
    
    autocorrelation_res <- remove_autocorrelations(df_sub, predictor_cols,
                                                   automatic_removal = automatic_removal,
                                                   autocorrelation_threshold = autocorrelation_threshold,
                                                   correlation_method = correlation_method)
    
    predictor_cols <- autocorrelation_res$predictor_cols
    autocorrelations <- autocorrelation_res$autocorrelations
    if (nrow(autocorrelations) > 0) {
      autocorrelations$response_col <- response_col
    }
    autocorrelations_overview[[process_code]] <- autocorrelations
    
    frm <- as.formula(paste(response_col, "~", paste(predictor_cols, collapse='+')))
    dat_cocc <- df_sub[, c('intergeneric_confusion_species', 'train_num_samples_ds_species',
                       'train_median_species_ds_gene_length', 'val_median_species_ds_gene_length',
                       'train_min_species_ds_gene_length', 'train_max_species_ds_gene_length')]
    glm_model <- eval(bquote(glm(as.formula(frm), data=.(df_sub), family = .(glm_family), 
                                 control = glm.control(maxit = .(glm_iter)))))
    
    simplified_model_info <- simplify_model(glm_model, aic_bic=aic_bic, 
                                            model_summaries=model_summaries)
    model_histories[[process_code]] <- simplified_model_info$history
    final_model <- simplified_model_info$final_model
    
    if (all(is.na(final_model))) {
      return(list(
        optimization_history = model_histories,
        summary = NA,
        autocorrelations = autocorrelations
      ))
    }
    
    all_sign_vars <- c(final_model$significant_variables, final_model$marginally_significant_variables)
    all_sign_vars <- all_sign_vars[!is.na(all_sign_vars)]
    
    if (length(all_sign_vars) != 0) {
      models_overview <- as.data.frame(coef(summary(final_model$model))) %>%
        rownames_to_column(var='predictor') %>%
        mutate(response = response_col, .before = predictor) %>%
        filter(predictor != '(Intercept)')
      
      if (plot) { plot_glm_vars(df_sub, response_col, all_sign_vars, glm_family=glm_family,
                glm_iter=glm_iter) }
      
      final_model_overviews[[process_code]] <- models_overview
    }
  }
  
  models_overview <- do.call(rbind, final_model_overviews)
  autocorrelations <- do.call(rbind, autocorrelations_overview)
  
  p_value_col <- tail(colnames(models_overview)[-1], n=1)
  
  models_overview <- models_overview %>%
    mutate_if(is.numeric, round, digits=round_p) %>%
    mutate(significance = case_when(
      !!ensym(p_value_col) < 0.0001 ~ '****',
      !!ensym(p_value_col) < 0.001 ~ '***',
      !!ensym(p_value_col) < 0.01 ~ '**',
      !!ensym(p_value_col) < 0.05 ~ '*',
      TRUE ~ 'ns'
    ))
  
  return(list(
    optimization_history = model_histories,
    summary = models_overview,
    autocorrelations = autocorrelations
  ))
}


#' GLM response~predictor plotting
#'
#' Plots response~predictor relationship
#' @param df A (wide) dataframe with predictors and response as columns
#' @param response The response columns within the df
#' @param predictors A list of predictors (i.e. column names)
#' @param glm_family The family used for glm calculation
#' @param glm_iter The number of max iterations used for the glm calculation
#' @return Nothing. Plots relationships in right bottom panel
#' @examples
#' data("dataset_info");
#' coefficients <- c('val_min_gene_length', 'val_max_gene_length',
#'                   'train_min_gene_length', 'train_max_gene_length',
#'                   'train_num_samples', 'val_num_samples',
#'                   'mean_gene_length', 'val_mean_gene_length',
#'                   'train_mean_gene_length');
#' plot_glm_vars(dataset_info, "species_confusion", coefficients);
#' @export
plot_glm_vars <- function(df, response, predictors, glm_family='binomial', 
                      glm_iter=1000) {
  for (predictor in predictors) {
    frm <- as.formula(paste(response, "~", predictor))
    glm_model <- eval(bquote(glm(.(frm), data=.(df), 
                                 family = .(glm_family), 
                                 control = glm.control(maxit = .(glm_iter)))))
    
    estimates <- as.data.frame(coef(summary(glm_model)))
    
    min_val <- min(df[[predictor]])
    max_val <- max(df[[predictor]])
    
    intercept <- glm_model$coefficients[['(Intercept)']]
    estimate <- round(estimates[rownames(estimates) == predictor, 'Estimate'], 2)
    p_value <- round(estimates[rownames(estimates) == predictor, ncol(estimates)], 5)
    sign <- case_when(
      p_value < 0.0001 ~ '****',
      p_value < 0.001 ~ '***',
      p_value < 0.01 ~ '**',
      p_value < 0.05 ~ '*',
      TRUE ~ 'ns'
    )
    
    pred1 <- data.frame(predictor = seq(from = min_val,
                                        to = max_val, by = 2))
    colnames(pred1) <- c(predictor)
    
    pred <- predict(glm_model, newdata = pred1, type = "response")
    
    range_x <- max_val - min_val
    y_max <- max(df[[response]])
    range_y <- y_max - min(df[[response]])
    
    text_x = min_val + (range_x * .01)
    text_y = max(df[[response]]) - (range_y * .05)
    
    # TODO: how to plot as example
    plot(x = df[[predictor]], y = df[[response]], 
         xlab = predictor, ylab=str_interp("${response} values"), 
         main=str_interp('${predictor}_${response}'), pch=16, lwd=2, lty=2)
    lines(pred1[[predictor]], pred, lwd=3)
    rect(
      xleft = text_x - range_x * 0.02,
      ybottom = text_y - range_y * 0.08,
      xright = text_x + range_x * 0.32,
      ytop = text_y + range_y * 0.08,
      col = adjustcolor("#91BAB6", alpha.f = 0.3),
      border = NA
    )
    text(
      x = text_x,
      y = text_y,
      labels = str_interp('estimate=${estimate}\np=${p_value}, ${sign}'),
      adj=0
    )
  }
}