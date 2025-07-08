#' This is data to be included in my package
#'
#' @name dataset_info
#' @docType data
#' @author Lara Kösters \email{lkoesters@@bgc-jena.mpg.de}
#' @author Kevin Karbstein \email{kkarb@@bgc-jena.mpg.de}
#' @keywords data
NULL

#' Autocorrelation identification and removal
#'
#' Identify autocorrelated variables. Automatically remove autocorrelated variables if specified
#' @param df Dataframe with response and predictors as columns
#' @param coefficients List of predictors sorted by their relevance (most to least)
#' @param ... Arguments given directly to 'cor' method
#' @param automatic_removal Determines whether to automatically remove autocorrelations
#' @param autocorrelation_threshold Threshold at which two variables are considered autocorrelated
#' @param correlation_method Method used for correlation calculation
#' @return Named list with a) a vector containing all removed predictors (empty if none were removed), and b) a dataframe containing autocorrelations and information on removed variables
remove_autocorrelations <- function(df, coefficients, ..., automatic_removal=TRUE,
                                    autocorrelation_threshold=0.7,
                                    correlation_method='pearson') {
  removed_coefficients = vector()
  cor_args = list(x=df[,coefficients], method=correlation_method)
  cor_args = c(cor_args, list(...))
  correlations <- as.data.frame(do.call(cor, cor_args))
  # compute p-values
  correlation_p_values <- as.data.frame(cor.mtest(correlations)$p)
  
  correlations_l = correlations %>%
    rownames_to_column(var='coefficientA') %>%
    pivot_longer(!coefficientA, names_to = "coefficientB", values_to = "correlation")
  
  correlations_p_values_l = correlation_p_values %>%
    rownames_to_column(var='coefficientA') %>%
    pivot_longer(!coefficientA, names_to = "coefficientB", values_to = "p_value")
  
  correlations_complete = merge(correlations_l, correlations_p_values_l, by=c('coefficientA', 'coefficientB'))
  
  correlations_complete = correlations_complete %>%
    filter(coefficientA != coefficientB) %>%
    mutate(col1 = pmin(coefficientA, coefficientB), col2 = pmax(coefficientA, coefficientB), 
           comparison = paste(col1, col2)) %>%
    distinct(comparison, .keep_all = TRUE) %>%
    select(!c(col1, col2, comparison))
  
  autocorrelations = correlations_complete[((correlations_complete$cor >= autocorrelation_threshold) | 
                                               (correlations_complete$cor <= -autocorrelation_threshold)) & 
                                              (correlations_complete$p_value < 0.05),]
  autocorrelations = autocorrelations %>%
    add_column(note = NA)
  
  if (nrow(autocorrelations) > 0) {
    if (automatic_removal) {
      coefficients_df = data.frame(idx = 1:length(coefficients),
                                    coefficient = coefficients)
      
      autocorrelations = merge(autocorrelations, coefficients_df[c('coefficient', 'idx')], by.x='coefficientA', by.y='coefficient')
      autocorrelations = merge(autocorrelations, coefficients_df[c('coefficient', 'idx')], by.x='coefficientB', by.y='coefficient', suffixes = c('1', '2'))
      autocorrelations = autocorrelations %>%
        rowwise() %>%
        mutate(idx_smaller = sort(c(idx1, idx2))[[1]],
               idx_bigger = sort(c(idx1, idx2))[[2]],
               comb = paste(idx_smaller, idx_bigger, sep='|')) %>%
        arrange(desc(idx_bigger))
      
      for (i in 1:nrow(autocorrelations)) {
        # C is the least important and part of comparison
        # B is the more important of comparison
        # A is more important than B
        
        # check for following constellation: A==B, B==C, A!=C
        autocor_row = autocorrelations[i,]
        c = autocor_row$idx_bigger
        b = autocor_row$idx_smaller
        
        b_to_a <- autocorrelations[autocorrelations$idx_bigger == b,]
        if (nrow(b_to_a) > 0) {  # check if there's at least 1 A
          for (a in b_to_a$idx_smaller) {  # for each variable that is correlated to and more important than B
            coefficient_a <- coefficients_df[a, 'coefficient']
            coefficient_b <- coefficients_df[b, 'coefficient']
            coefficient_c <- coefficients_df[c, 'coefficient']
            
            if(nrow(autocorrelations[(autocorrelations$idx_bigger == autocor_row$idx_bigger) & (autocorrelations$idx_smaller == a),]) == 0) {
              # A!=C but A==B and B==C: remove B
              coefficients = coefficients[!coefficients == coefficient_b]
              removed_coefficients = append(removed_coefficients, coefficient_b)
              autocorrelations[i, 'note'] = str_interp('${coefficient_a}!=${coefficient_c}, but ${coefficient_a}==${coefficient_b} and ${coefficient_b}==${coefficient_c}; removed ${coefficient_b}')
            } else {
              coefficients = coefficients[!coefficients == coefficient_c]
              removed_coefficients = append(removed_coefficients, coefficient_c)
              autocorrelations[i, 'note'] = str_interp('removed ${coefficient_c}')
            }
          }
        } else {
          coefficient_c = coefficients_df[c, 'coefficient']
          coefficients = coefficients[!coefficients == coefficient_c]
          removed_coefficients = append(removed_coefficients, coefficient_c)
          autocorrelations[i, 'note'] = str_interp('removed ${coefficient_c}')
        }
      }
    } else {
      f = str_interp("autocorrelations_${format(Sys.time(), '%m%d%Y_%H%M')}.tsv")
      write.table(autocorrelations, f, 
                  sep = '\t', quote = FALSE, row.names = FALSE)
      stop(str_interp("Some of your variables are autocorrelated. Check ${f} for more info"))
    }
  }
  
  res = autocorrelations[,c('coefficientA', 'coefficientB', 'correlation', 'p_value', 'note')]
  return(list("removed_predictor_cols" = removed_coefficients, "autocorrelations" = autocorrelations[,c('coefficientA', 'coefficientB', 'correlation', 'p_value', 'note')]))
}

#' Model variables categorical check
#'
#' Helper method used for expanding the model summary by adding information on whether or not model variables are categorical
#' @param df Dataframe with response and predictors as columns
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @param col Column to be used for check. Default is 'coefficients'
#' @return Language object with if/else statements to be applied to model summary with information on whether model variables are categorical
#' @examples
#' data("plants");
#' 
#' setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor))
#' @export
setup_categorical_check <- function(df, categorical_vars, col=quote(coefficients)) {
  if (length(categorical_vars) > 0) {
    categorical_check = c()
    for (cat_var in categorical_vars) {
      cat_vals = as.character(unique(df[[cat_var]]))
      categorical_check = append(categorical_check, purrr::map(cat_vals,
                                                                ~quo(str_detect(!!col, regex(paste0(!!cat_var, !!.x, "(?![a-z0-9_\\-#])"))) ~ !!cat_var)))
    }
  }
  else categorical_check = list(rlang::expr(TRUE ~ NA), rlang::expr(TRUE ~ NA))
  return(categorical_check)
}

#' Add model assessments
#'
#' Add evaluations of model given model and evaluation methods (anova/aic/aicc/bic)
#' @param regression_model Regression model that is assessed
#' @param evaluation_methods List of methods to use for model assessment. Options are: anova, aic, aicc, and bic
#' @param is_lmer Boolean stating whether model is of type lmer
#' @return Regression model with added assessments based on assessment methods provided
add_assessments <- function(regression_model, evaluation_methods, is_lmer=FALSE) {
  if ('aic' %in% evaluation_methods) {
    if (is_lmer) attr(regression_model, 'aic') = AIC(regression_model)
    else regression_model$aic = AIC(regression_model)
  }
  if ('aicc' %in% evaluation_methods) {
    if (is_lmer) attr(regression_model, 'aicc') = AICc(regression_model)
    else regression_model$aicc = AICc(regression_model)
  }
  if ('bic' %in% evaluation_methods) {
    if (is_lmer) attr(regression_model, 'bic') = BIC(regression_model)
    else regression_model$bic = BIC(regression_model)
  }
  
  return(regression_model)
}

#' Regression model was improved
#'
#' Determine if model 1 better fits the data than model 2 given two models and evaluation methods (anova/aic/aicc/bic)
#' @param regression_model Regression model 1
#' @param old_regression_model Regression model 2
#' @param evaluation_methods List of methods to use for model assessment. Options are: anova, aic, aicc, and bic
#' @param direction Mode of stepwise model improvement. Either 'forward' (i.e., forward selection), or 'backward' (i.e., backward simplification)
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @return Boolean; TRUE if model 1 better fits the data than model 2
model_improved <- function(regression_model, old_regression_model, evaluation_methods, direction, model_type) {
  if ('anova' %in% evaluation_methods) {
    if (model_type == 'gam') anova_res = anova.gam(regression_model, old_regression_model)
    else anova_res = anova(regression_model, old_regression_model)
    
    p_col = colnames(anova_res)[[grep('Pr\\(', colnames(anova_res))]]
    if ((direction == 'backward') && (anova_res[2, p_col] < .05)) return(FALSE)
    else if ((direction == 'forward') && (anova_res[2, p_col] > .05)) return(FALSE)
  }
  
  is_lmer = model_type == 'lmer' || model_type == 'glmer'
  
  if ('aic' %in% evaluation_methods) {
    if (is_lmer) {
      aic_new = attr(regression_model, 'aic')
      aic_old = attr(old_regression_model, 'aic')
    }
    else {
      aic_new = regression_model$aic
      aic_old = old_regression_model$aic
    }
    
    if (aic_new >= aic_old) return(FALSE)
  }
  if ('aicc' %in% evaluation_methods) {
    if (is_lmer) {
      aicc_new = attr(regression_model, 'aicc')
      aicc_old = attr(old_regression_model, 'aicc')
    }
    else {
      aicc_new = regression_model$aicc
      aicc_old = old_regression_model$aicc
    }
    
    if (aicc_new >= aicc_old) return(FALSE)
  }
  if ('bic' %in% evaluation_methods) {
    if (is_lmer) {
      bic_new = attr(regression_model, 'bic')
      bic_old = attr(old_regression_model, 'bic')
    }
    else {
      bic_new = regression_model$bic
      bic_old = old_regression_model$bic
    }
    
    if (bic_new >= bic_old) return(FALSE)
  }
  
  return(TRUE)
}

#' Expansion of model summary
#'
#' Expand model summary by adding column names and information on interaction groups
#' @param df Dataframe with response and predictors as columns
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param regression_model Regression model
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @return Expanded summary of a regression model with added information on base variables and variable type
expanded_model_summary <- function(df, model_type, regression_model, categorical_vars) {
  if (model_type == 'gam') model_sum = summary(regression_model)$p.table
  else model_sum = coef(summary(regression_model))
  
  if ('Pr(>|z|)' %in% colnames(model_sum)) model_sum = model_sum[,'Pr(>|z|)']
  else if ('Pr(>|t|)' %in% colnames(model_sum)) model_sum = model_sum[,'Pr(>|t|)']
  else stop(str_interp("Not able to calculate p-values."))
  
  p_values_model = data.frame(p_values = model_sum) %>%
    rownames_to_column(var='coefficients') %>%
    filter(coefficients != '(Intercept)') %>% # rownames(.)
    mutate(is_interaction = case_when(str_detect(coefficients, ':') ~ TRUE,
                                      TRUE ~ FALSE)) %>%
    mutate(inter_var1 = case_when(is_interaction ~ str_extract(coefficients, "^[^:]+")),
           inter_var2 = case_when(is_interaction ~ str_extract(coefficients, "[^:]+$"))) %>%
    mutate(inter_var1 = case_when(!!!setup_categorical_check(df, categorical_vars, quote(inter_var1))),
           inter_var2 = case_when(!!!setup_categorical_check(df, categorical_vars, quote(inter_var2)))) %>%
    mutate(inter_var1 = case_when(is_interaction & is.na(inter_var1) ~ str_extract(coefficients, "^[^:]+"),
                                  TRUE ~ inter_var1),
           inter_var2 = case_when(is_interaction & is.na(inter_var2) ~ str_extract(coefficients, "[^:]+$"),
                                  TRUE ~ inter_var2)) %>%
    mutate(cat_var = case_when(!!!setup_categorical_check(df, categorical_vars))) %>%
    mutate(var_type = ifelse(is.na(cat_var), 'cont', 'cat')) %>%
    mutate(cat_var = ifelse(is.na(cat_var), coefficients, cat_var)) %>%
    mutate(cat_var = case_when(is_interaction ~ paste0(inter_var1, ':', inter_var2),
                               TRUE ~ cat_var),
           inter_var1 = ifelse(is.na(inter_var1), cat_var, inter_var1),
           inter_var2 = ifelse(is.na(inter_var2), cat_var, inter_var2)) %>%
    mutate(cat_var = ifelse(grepl(')[0-9]+$', coefficients), gsub('[0-9]+$', '', coefficients), cat_var)) # TEST
  
  p_values_model = p_values_model %>%
    rowwise() %>%
    mutate(smallest_p_value = min(p_values_model[(p_values_model['inter_var1'] == cat_var) | (p_values_model['inter_var2'] == cat_var) | (p_values_model['cat_var'] == cat_var), 'p_values']))
  
  non_sign_interactions = p_values_model[(p_values_model$p_values > 0.1) & p_values_model$is_interaction, 'coefficients']

  if (nrow(non_sign_interactions) > 0) {
    p_values_model = p_values_model %>%
      arrange(desc(is_interaction), desc(smallest_p_value))
  } else {
    p_values_model = p_values_model %>%
      arrange(desc(smallest_p_value))
  }
  
  return(p_values_model)
}

#' Generate regression model
#'
#' Helper method to generate a regression model based on a given model type, family and formula
#' @param df Dataframe with response and predictors as columns
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param model_family The family used for glm calculation (default: gaussian)
#' @param ... Arguments given to model call
#' @return List with regression model and used model arguments
generate_regression_model <- function(df, model_type, term, model_family, ...) {
  errors=c()
  
  if (model_type == 'lmer') model_args = list(formula=term, data=df, REML=FALSE)
  else if ((model_type == 'lm') || (model_type == 'lmer') || (model_type == 'nls')) model_args = list(formula=term, data=df)
  else if (model_type == 'nlmer') model_args = list(data=df)
  else model_args = list(formula=term, data=df, family=model_family)
  model_args = c(model_args, list(...))
  if (model_type == 'glm') regression_model = withCallingHandlers(do.call(glm, model_args),
                                                                  warning = function(w){
                                                                    if (grepl("non-integer #successes in a binomial glm", w$message)) {
                                                                      errors = append(errors, w$message)
                                                                      tryInvokeRestart("muffleWarning")
                                                                    }
                                                                    else warning(w$message)
                                                                  })
  else if (model_type == 'lm') regression_model = do.call(lm, model_args)
  else if (model_type == 'glmer') regression_model = try(do.call(glmer, model_args), silent=TRUE)
  else if (model_type == 'lmer') regression_model = do.call(lmer, model_args)
  else if (model_type == 'gam') regression_model = do.call(gam, model_args)
  else if (model_type == 'nlmer') {
    parameter_names = names(model_args)
    parameter_names[parameter_names == 'non_linear'] = 'model'
    names(model_args) = parameter_names
    regression_model = do.call(nlme, model_args)
  } else if (model_type == 'nls') regression_model = do.call(nls, model_args)
  else stop(str_interp("Unknown model type ${model_type}. Please choose either (g/n)(l/a)m(er)."))

  if ('try-error' %in% class(regression_model)) {
    e = attr(regression_model, 'condition')
    if (grepl("converge", e$message)) {
      warning(str_interp("${e$message}. Please change the random effect or x-variable set."))
      errors = append(errors, str_interp("${e$message}. Please change the random effect or x-variable set."))
      return(list(regression_model=regression_model, model_args=model_args, errors=errors))
    }
  }
  
  return(list(regression_model=regression_model, model_args=model_args, errors=errors))
}

#' Forward model selection
#'
#' Forward model selection
#' @param df Dataframe with response and predictors as columns
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param evaluation_methods Methods to be used for model evaluation (options: anova, aic, aicc, bic, default: anova)
#' @param ... Parameters to be directly used with model call
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @param model_family The family used for glm calculation (default: gaussian)
#' @param trace Store and return model selection history (default: FALSE)
#' @param omit.na Either 'overall' or 'stepwise'. If 'overall', NAs are removed before modeling. If 'stepwise', NAs are removed per step based on the variables in the current formula
#' @return List containing the final regression model, the significant and marginally significant model variables and the selection history if trace is TRUE
forward_selection <- function(df, model_type, term, evaluation_methods, categorical_vars, ...,
                                   model_family='binomial', trace=FALSE, omit.na='overall') {
  errors = c()
  single_vars = c()
  interactions = c()
  
  split_frm = find_call(term, return='fixed')

  for (expr in split_frm) {
    if (str_detect(deparse(expr), ':')) interactions = append(interactions, deparse(expr))
    else single_vars = append(single_vars, deparse(expr))
  }
  
  if (omit.na == 'overall') df = remove_nas(df, term)
  
  frm = update(as.formula(term), '. ~ 1')
  regression_model = generate_regression_model(df, model_type, term, model_family, ...)
  model_args = regression_model$model_args
  regression_model = regression_model$regression_model
  errors = append(errors, regression_model$errors)
  
  if (trace) history = list()
  
  blank_start = TRUE
  for (vars in list(single_vars, interactions)) {
    simplify = TRUE
    while (simplify) {
      res = mo_step(df, regression_model, vars, 'forward', evaluation_methods, categorical_vars, model_type, blank_start=blank_start, stepwise_omit=(omit.na == 'stepwise'), model_args=model_args)
      errors = append(errors, res$errores)
      if (blank_start) blank_start = FALSE
      simplify = res$simplify
      regression_model = res$regression_model
      vars = res$vars
      if (trace && simplify) history[[str_interp("${length(history)+1}_${res$var}")]] = regression_model
    }
  }
  
  model_summary = expanded_model_summary(df, model_type, regression_model, categorical_vars)
  main_effects = model_summary[model_summary$is_interaction, c('inter_var1', 'inter_var2')]
  main_effects = c(main_effects$inter_var1, main_effects$inter_var2)
  
  for (main_effect in main_effects) {
    if (!(main_effect %in% model_summary$coefficients)) {
      if (model_type == 'lmer') {
        # update formula
        model_args$formula = update(formula(regression_model), paste('. ~ . +', main_effect))
        
        # update data
        if (omit.na == 'stepwise') model_args$data = remove_nas(model_args$data, model_args$formula)
        
        regression_model = do.call(lmer, model_args)
      }
      else regression_model = withCallingHandlers(update(regression_model, paste('~ . +', main_effect)),
                                                  warning = function(w){
                                                    if (grepl("non-integer #successes in a binomial glm", w$message)) {
                                                      errors = append(errors, w$message)
                                                      tryInvokeRestart("muffleWarning")
                                                    }
                                                    else warning(w$message)
                                                  })
      if (trace) history[[str_interp("${length(history)+1}_${main_effect}")]] = regression_model
    }
  }
  
  out = list(final_model = regression_model,
             significant_variables = model_summary[model_summary$p_values < 0.05, 'coefficients'],
             marginally_significant_variables = model_summary[(model_summary$p_values < 0.1) & (model_summary$p_values >= 0.05), 'coefficients'])
  
  if (trace) out$history = history
  out$errors = errors
  
  return(out)
}

#' Take step in model simplification/selection
#'
#' Used to either remove or add a variable to an already existing regression model
#' @param df Dataframe with response and predictors as columns
#' @param regression_model The current regression model
#' @param vars List of model variables that are yet to be added/removed
#' @param direction Mode of stepwise model improvement. Either 'forward' (i.e., forward selection), or 'backward' (i.e., backward simplification)
#' @param evaluation_methods Methods to be used for model evaluation (options: anova, aic, aicc, bic, default: anova)
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param blank_start Used for first run in forward model selection. Creates empty model
#' @param stepwise_omit If TRUE, NAs are removed per step based on the variables in the current formula
#' @param model_args Arguments given directly to model call
#' @return List containing a) information on whether model was simplified/expanded, b) the (new) current regression model, c) updated list of variables to be added/removed, d) added/removed variable, and e) summary of the model
mo_step <- function(df, regression_model, vars, direction, evaluation_methods, categorical_vars, model_type,
                    blank_start=FALSE, stepwise_omit=FALSE, model_args=NA) {
  errors = c()
  is_lmer = model_type == 'lmer' || model_type == 'glmer'
  final_model = NA
  adjusted_var = NA
  model_summary = NA
  min_p_value = Inf
  reason = NA
  simplified = FALSE
  
  interactions = unique(vars[str_detect(vars, ':')])
  for (var in vars) {  # loop important for forward simplification
    # add or remove var from model
    if (direction == 'forward') {
      if (blank_start) {
        random_effects = lapply(findbars(formula(regression_model)), function(y) paste0('(', deparse(y), ')'))
        if (length(random_effects) > 0) addon = paste(paste(random_effects, sep=' + '), ' + ')
        else addon = ''
        d = paste('. ~ ', addon)
      }
      else d = '. ~ . + '
    } else d = '. ~ . -'
    
    # check for interactions when removing single vars and remove them as well
    var_m = var
    if (direction == 'backward' && !str_detect(var, ':')) {
      for (interaction in interactions) {
        if (str_detect(interaction, fixed(str_interp('${var}:'))) || str_detect(interaction, fixed(str_interp(':${var}')))) {
          var_m = paste(var_m, '-', interaction)
        }
      }
    }
    
    # update model
    if (direction == 'forward') {
      # stepwise: more restrictive model is regression_model_updated
      if (is_lmer) {
        # update formula
        model_args$formula = update(formula(regression_model), paste(d, var_m))
        
        # update data
        if (stepwise_omit) model_args$data = remove_nas(model_args$data, model_args$formula)
        
        regression_model_updated = do.call(lmer, model_args)
        
        if (stepwise_omit) {
          model_args$formula = formula(regression_model)
          regression_model = do.call(lmer, model_args)
        }
      } else {
        regression_model_updated = withCallingHandlers(update(regression_model, paste(d, var_m), data=df),
                                                       warning = function(w){
                                                         if (grepl("non-integer #successes in a binomial glm", w$message)) {
                                                           errors = append(errors, w$message)
                                                           tryInvokeRestart("muffleWarning")
                                                         }
                                                         else warning(w$message)
                                                       })
        
        # if interaction, make sure order is same as before
        if (str_detect(var, ':')) {
          if (!str_detect(Reduce(paste, deparse(regression_model_updated$formula)), var)) {
            vars = vars[vars != var]
            var_spl = as.vector(unlist(str_split(var, ':')))
            var = paste0(var_spl[[2]], ':', var_spl[[1]])
            vars = c(vars, var)
          }
        }
        
        if (stepwise_omit) {
          tryCatch(
            expr = {regression_model = update(regression_model, '. ~ .', data=regression_model_updated$model)},
            error = function(e) {regression_model = update(regression_model, '. ~ .', data=regression_model_updated$data)}
          )
        }
      }
    } else {
      # stepwise: more restrictive model is regression_model
      if (is_lmer) {
        # update formula
        model_args$formula = update(formula(regression_model), paste(d, var_m))
        
        # update data
        if (stepwise_omit) model_args$data = remove_nas(model_args$data, formula(regression_model))
        
        if (model_type == 'glmer') regression_model_updated = try(do.call(glmer, model_args), silent=TRUE)
        else regression_model_updated = do.call(lmer, model_args)
        
        if ('try-error' %in% class(regression_model_updated)) {
          e = attr(regression_model_updated, 'condition')
          if (grepl("converge", e$message)) {
            warning(str_interp("${e$message}. Please change the random effect or x-variable set."))
            errors = append(errors, str_interp("${e$message}. Please change the random effect or x-variable set."))
            return(list(simplify=simplified, regression_model=regression_model, vars=vars, model_summary=model_summary, var=adjusted_var, reason=reason, errors=errors))
          }
        }
      }
      else {
        regression_model_updated = withCallingHandlers(update(regression_model, paste(d, var_m), data=regression_model$model),
                                                          warning = function(w){
                                                            if (grepl("non-integer #successes in a binomial glm", w$message)) {
                                                              errors = append(errors, w$message)
                                                              tryInvokeRestart("muffleWarning")
                                                            }
                                                            else warning(w$message)
                                                          })
      }
    }
    model_summary = expanded_model_summary(df, model_type, regression_model_updated, categorical_vars)
    regression_model_updated = add_assessments(regression_model_updated, evaluation_methods, model_type == 'lmer')
    
    # first check if variables are significant
    if (((direction == 'forward') && (model_summary[1, 'smallest_p_value'] > 0.1)) || 
        !blank_start && !model_improved(regression_model_updated, regression_model, evaluation_methods, direction, model_type)) {
      reason = 'the p-value was not significant'
      next
    }
    this_p_value = model_summary[model_summary$cat_var == var,][['smallest_p_value']]
    if (length(this_p_value) > 1) this_p_value = this_p_value[1]
    
    if ((direction == 'backward') || (this_p_value < min_p_value)) {
      final_model = regression_model_updated
      adjusted_var = var
      min_p_value = this_p_value
      simplified = TRUE
    } else {
      reason = 'the p-value was not significant'
    }
    
    if (direction == 'backward') break
  }
  
  if (simplified) return(list(simplify=simplified, regression_model=final_model, vars=vars[vars != adjusted_var], var=adjusted_var, model_summary=model_summary, errors=errors))
  else return(list(simplify=simplified, regression_model=regression_model, vars=vars, model_summary=model_summary, var=adjusted_var, reason=reason, errors=errors))
}


#' Model simplification
#'
#' Simplify or expand a regression model given the data and model parameters
#' @param df Dataframe with response and predictors as columns
#' @param model_type The model to be used (options: (g/n)(l/a)m(er))
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param evaluation_methods Methods to be used for model evaluation (options: anova, aic, aicc, bic)
#' @param ... Parameters to be directly used with model call
#' @param model_family The family used for glm calculation (default: gaussian)
#' @param direction Mode of stepwise model improvement. Either 'forward' (i.e., forward selection), or 'backward' (i.e., backward simplification), or 'both' (default: both)
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @param backward_simplify_model If FALSE, the model and information on significant and marginally significant variables are returned without any improvements
#' @param trace Store and return model selection history (default: FALSE)
#' @param omit.na Either 'overall' or 'stepwise'. If 'overall', NAs are removed before modeling. If 'stepwise', NAs are removed per step based on the variables in the current formula
#' @return List of results from forward and/or backward model selection. Backward/forward result will have the following structure: List containing the final regression model, the significant and marginally significant model variables and the selection history if trace is TRUE
#' @examples
#' data("plants");
#' 
#' simplify_model(plants, "glm", quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + 
#' isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + 
#' solar_radiation:isothermality + annual_mean_temperature:isothermality),
#' c("anova"), model_family='quasibinomial', direction='backward',
#' categorical_vars=c('habitat', 'ploidy'), backward_simplify_model=TRUE, trace=TRUE, omit.na='overall')
#' @export
simplify_model <- function(df, model_type, term, evaluation_methods, ..., model_family='binomial', 
                           direction='both', categorical_vars=NA, backward_simplify_model=TRUE, trace=FALSE, omit.na='overall') {
  models = list()
  errors = c()
  
  if (model_type == 'nlmer') {
    generated_model = generate_regression_model(df, model_type, term, model_family, ...)
    errors = append(errors, generated_model$errors)
    model_args = generated_model$model_args
    regression_model = generated_model$regression_model
    regression_model = add_assessments(regression_model, evaluation_methods, FALSE)
    anova_res = anova.lme(regression_model)
    regression_model$anova = anova_res
    regression_model$aic = AIC(regression_model)
    regression_model$bic = BIC(regression_model)
    regression_model$aicc = AICc(regression_model)
    
    models$nlmer = list(final_model = regression_model)
    if (trace) models$nlmer$history = NA
    return(models)
  }
  
  if (direction == 'both' | direction == 'forward') {
    models$forward = forward_selection(df, model_type, term, evaluation_methods, categorical_vars, ..., model_family=model_family, trace=trace, omit.na=omit.na)
    errors = append(errors, models$forward$erorrs)
  }
  if (direction == 'both' | direction == 'backward') {
    models$backward = backward_simplification(df, model_type, term, evaluation_methods, categorical_vars, ..., model_family=model_family, simplify_model=backward_simplify_model, trace=trace, omit.na=omit.na)
    errors = append(errors, models$backward$erorrs)
  }
  
  return(models)
}

#' Remove NAs
#'
#' Small helper method to remove NAs from dataframe given model formula
#' @param df Dataframe with response and predictors as columns
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @return List containing the final regression model, the significant and marginally significant model variables and the selection history if trace is TRUE
remove_nas <- function(df, term) {
  columns_in_formula = as.character(find_call(term, return='atomic', df_cols=colnames(df)))
  response = terms(as.formula(term))[[2]]
  cols_w_y = append(as.character(response), columns_in_formula)
  df = df[,cols_w_y]
  df = na.omit(df)
}

#' Backward regression model simplification
#'
#' Simplify the regression model by eliminating model variables one by one
#' @param df Dataframe with response and predictors as columns
#' @param model_type The model to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param evaluation_methods Methods to be used for model evaluation (options: anova, aic, aicc, bic, default: anova)
#' @param categorical_vars List of categorical variables within the model (base names, i.e., names of columns)
#' @param ... Parameters to be directly used with model call
#' @param model_family The family used for glm calculation (default: gaussian)
#' @param simplify_model If FALSE, the model and information on significant and marginally significant variables are returned without any improvements
#' @param trace Store and return model selection history (default: FALSE)
#' @param omit.na Either 'overall' or 'stepwise'. If 'overall', NAs are removed before modeling. If 'stepwise', NAs are removed per step based on the variables in the current formula
#' @return List containing the final regression model, the significant and marginally significant model variables and the selection history if trace is TRUE
backward_simplification <- function(df, model_type, term, evaluation_methods, categorical_vars, ..., model_family='binomial',
                          simplify_model=TRUE, trace=FALSE, omit.na='overall') {
  if (omit.na == 'overall') df = remove_nas(df, term)
  
  errors = c()
  generated_model = generate_regression_model(df, model_type, term, model_family, ...)
  model_args = generated_model$model_args
  regression_model = generated_model$regression_model
  
  if (trace) {
    history = list()
    history[[str_interp("${length(history)+1}_starting_model")]] = regression_model
  }
  out = list('final_model' = NA)
  model_summary = expanded_model_summary(df, model_type, regression_model, categorical_vars)
  
  if ((model_summary[1, 'smallest_p_value'] < 0.1) || (!simplify_model)) {
    # All p-values are below p<0.1 --> accept the original model
    out$significant_variables = model_summary[model_summary$p_values < 0.05, 'coefficients']
    out$marginally_significant_variables = model_summary[(model_summary$p_values < 0.1) & (model_summary$p_values >= 0.05), 'coefficients']
    out$final_model = regression_model
  } else {
    regression_model = add_assessments(regression_model, evaluation_methods, model_type == 'lmer')
    simplify = TRUE
    while (simplify) {
      res = mo_step(df, regression_model, unique(model_summary[model_summary['smallest_p_value'] > 0.1,][['cat_var']]), 'backward', evaluation_methods, categorical_vars, model_type, model_args = model_args)
      simplify = res$simplify
      errors = append(errors, res$errors)
      if (simplify) {
        regression_model = res$regression_model
        if (omit.na == 'stepwise') {
          current_term = formula(regression_model)
          regression_model = generate_regression_model(df, model_type, current_term, model_family, ...)
          regression_model = add_assessments(regression_model, evaluation_methods, model_type == 'lmer')
          model_summary = expanded_model_summary(df, model_type, regression_model, categorical_vars)
        } else {
          model_summary = res$model_summary
        }
        if (trace) history[[str_interp("${length(history)+1}_${res$var}")]] = regression_model
      }
    }
    var_cnt = length(find_call(term, return='fixed'))
    if ((var_cnt > 2) && (formula(regression_model) == term)) {
      warning(str_interp("We initialised the model successfully, but weren't able to simplify it, because ${res$reason}. Please check overfitting and multicollinearity among predictors to exclude potential modeling issues"))
    }
    
    out$final_model = regression_model
    out$significant_variables = model_summary[model_summary$p_values < 0.05, 'coefficients']
    out$marginally_significant_variables = model_summary[(model_summary$p_values < 0.1) & (model_summary$p_values >= 0.05), 'coefficients']
    if (trace) out$history = history
  }
  
  out$errors = errors
  return(out)
}


#' Split and search model formula
#'
#' Extract part of model formula
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param return Can be 'all', 'atomic', 'interactions', 'main_effects', 'negative', or 'match'. Default: 'all'
#' @param pred_full Only needed for 'match' and 'negative'. The term to search for within the formula
#' @param df_cols Only needed if return == 'atomic'
#' @return The extracted part(s) of the model formula matching a given pattern. Can be either a call or a vector
find_call <- function(term, pred_full=NA, return='all', df_cols=NA) {
  expr_terms = terms(as.formula(term))
  term_labels = colnames(attr(expr_terms, 'factors'))
  
  if (return == 'all') {  # return as is
    out = term_labels
  } else if (return == 'fixed') {
    randoms = findbars(as.formula(term))
    vars = term_labels
    response = lhs(as.formula(term))
    if (length(response) > 1) {
      response = deparse(lhs(response))
      random_idx = sapply(vars, function(y) grepl(y, randoms, fixed=TRUE))
      out = vars[(vars != response) & !random_idx]
    } else {
      response = deparse(response)
      out = vars[vars != response]
    }
  } else if (return == 'random') {
    return(as.character(findbars(as.formula(term))))
  } else if (return == 'atomic') {  # return column names of fixed effects, i.e., remove transforms
    randoms = findbars(as.formula(term))
    vars = all.vars(as.formula(term))
    response = lhs(as.formula(term))
    if (length(response) > 1) {
      response = deparse(lhs(response))
    } else {
      response = deparse(response)
    }
    
    out = vars[vars != response]
    # out = c()
    # for (term_label in term_labels) {
    #   out = append(out, str_extract(term_label, regex('[^()\\^\\|\\s]+(?=\\)|$|\\^|,)')))
    # }
  } else if (return == 'interactions') {  # subset for interactions
    out = term_labels[str_detect(term_labels, ':')]
  } else if (return == 'main_effects') {  # subset for main effects
    out = term_labels[!str_detect(term_labels, ':')]
  } else if (return == 'negative') {  # get term labels that don't match search term
    out = term_labels[!str_detect(term_labels, fixed(pred_full))]
  } else if (return == 'match') {  # get term labels matching search term
    out = term_labels[str_detect(term_labels, fixed(pred_full))]
  }
  
  out_v = list()
  for (x in out) out_v = append(out_v, parse(text=x)[[1]])
  return(unique(out_v))
}

#' Determine regression model family
#'
#' Determine model family automatically based on response column in dataframe
#' @param df Dataframe with response and predictors as columns
#' @param response_frm Response formula
#' @return Model family to be used with regression model
determine_model_family = function(df, response_frm) {
  response_col = paste(response_frm)
  df[response_col] = with(df, eval(parse(text=response_frm)))
  
  response = df[!is.na(df[[response_col]]),][[response_col]]
  is_num = is.numeric(response)
  is_int = is.integer(response)
  
  if (is_int || is_num) {
    if ((min(response) < 0) || (max(response) > 100)) {
      model_family = c('gaussian')
    } else if (is_int || (is_num && ((min(response) < 0) || (max(response) > 1)))) {
      model_family = c('gaussian', 'poisson')
    } else if (is_num && (min(response) >= 0) && (max(response) <= 1)) {
      model_family = c('gaussian', 'poisson', 'quasibinomial')
    } else {
      model_family = c()
    }
  } else if (is.logical(response)) {
    model_family = c('binomial')
  } else {
    model_family = c()
  }
  
  return(model_family)
}

#' A minor version of expanded model summary
#'
#' Expand model summary by adding column names and information on interaction groups
#' @param model_summary Initial model summary
#' @param response_frm Response formula
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param categorical_check Output of [setup_categorical_check()]
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param df Dataframe with response and predictors as columns
#' @return Expanded summary of a regression model with added information on base variables, variable types and significance levels
#' @examples
#' data("plants")
#' 
#' model_family = 'quasibinomial'
#' categorical_check = setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor))
#' final_model = glm(sexual_seed_prop ~ altitude + solar_radiation + 
#' annual_mean_temperature + isothermality + habitat + ploidy + 
#' solar_radiation:isothermality, family=model_family, data=plants)
#' models_overview = coef(summary(final_model))
#' slightly_expanded_model_summary(models_overview, quote(sexual_seed_prop), quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + solar_radiation:isothermality + annual_mean_temperature:isothermality), categorical_check, 3, plants)
#' @export
slightly_expanded_model_summary <- function(model_summary, response_frm, term, categorical_check, round_p, df) {
  models_overview = as.data.frame(model_summary) %>%
    rownames_to_column(var='predictor') %>%
    filter(predictor != '(Intercept)') %>% # rownames(.)
    mutate(response = paste(response_frm), .before = "predictor") %>%
    mutate(predictor_col = case_when(!!!categorical_check), .after = "predictor") %>%
    mutate(var_type = ifelse(is.na(predictor_col), 'numeric', 'categorical'), .after = "predictor_col") %>%
    mutate(predictor_col = ifelse(is.na(predictor_col), predictor, predictor_col), .after = "predictor")
  p_value_col = tail(colnames(models_overview)[-1], n=1)
  
  split_frm = find_call(term, return='atomic', df_cols=colnames(df))
  
  for (el in split_frm) {
    ptrn = str_interp('(?<!:)${deparse(el)}(?![a-z0-9_\\-#:])')
    for (pred in models_overview[,'predictor']) {
      if (str_detect(pred, regex(ptrn, ignore_case = TRUE))) {
        models_overview[models_overview['predictor'] == pred, 'predictor_col'] = deparse(el)
      }
    }
  }
  
  models_overview = models_overview %>%
    mutate_if(is.numeric, round, digits=round_p) %>%
    mutate(significance = case_when(
      !!ensym(p_value_col) < 0.001 ~ '***',
      !!ensym(p_value_col) < 0.01 ~ '**',
      !!ensym(p_value_col) < 0.05 ~ '*',
      TRUE ~ 'ns'
    ))
    
  return(models_overview)
}


#' Parent function for model optimization
#'
#' Optimize model by removing autocorrelations and variables that do not significantly predict response variable
#' @param df Dataframe with response and predictors as columns
#' @param term The formula to be used with the model. Can be either quote() or formula()
#' @param ... Arguments given directly to model call
#' @param autocorrelation_cols Sorted list of columns in dataframe to be considered when removing autocorrelations. Should be sorted by priority. Last element gets eliminated first
#' @param automatic_removal Whether to automatically remove autocorrelations
#' @param autocorrelation_threshold Threshold at which two variables are considered autocorrelated
#' @param correlation_method The method used for correlation calculation
#' @param cor_use Parameter 'use' for 'cor' method. Describes handling of missing values
#' @param model_type Model type to be used (options: (g/n)(l/a)m(er), default: glm)
#' @param model_family Model family used for glm calculation (default: gaussian)
#' @param evaluation_methods Methods to be used for model evaluation (options: anova, aic, aicc, bic, default: anova)
#' @param simplification_direction Mode of stepwise model improvement. Either 'forward' (i.e., forward selection), or 'backward' (i.e., backward simplification), or 'both' (default: both)
#' @param omit.na Either 'overall' or 'stepwise'. If 'overall', NAs are removed before modeling. If 'stepwise', NAs are removed per step based on the variables in the current formula
#' @param scale_predictor Whether to apply scaling to predictor variables
#' @param plot_quality_assessment Module to use for plots for model quality assessment (options: performance, baseR)
#' @param plot_relationships Whether to plot a) regression, b) effect size, and c) estimates
#' @param jitter_plots Whether geom_point plots should use to jitter
#' @param plot_type Either 'boxplot' or 'violin'. Used to plot regression plots for categorical variables (default: 'boxplot')
#' @param stat_test Either 't.test' or 'wilcox'. Used to calculate statistics for regression plots of categorical variables (default: 'wilcox')
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param backward_simplify_model If FALSE, the model and information on significant and marginally significant variables are returned without backward simplification (default: TRUE)
#' @param trace Store and return model selection history (default: FALSE)
#' @return List with a) information on autocorrelated variables and b) final simplified/expanded models with further information (see function [simplify_model()] for further details)
#' @examples
#' data("plants");
#' optimize_model(plants, quote(sexual_seed_prop ~ altitude + latitude_gps_n + longitude_gps_e + (solar_radiation + annual_mean_temperature + isothermality)^2 + I(isothermality^2) + habitat + ploidy), 
#' autocorrelation_cols = c("solar_radiation", "annual_mean_temperature", "isothermality", "altitude", "latitude_gps_n", "longitude_gps_e"), 
#' automatic_removal=TRUE, autocorrelation_threshold = 0.8, correlation_method="spearman", model_type = "glm", model_family = "quasibinomial", 
#' evaluation_methods=c("anova"), simplification_direction="backward", omit.na="overall", scale_predictor=TRUE, plot_quality_assessment="performance", 
#' round_p=3, cor_use="complete.obs", plot_relationships=TRUE, jitter_plots=TRUE, plot_type="boxplot", stat_test="wilcox", backward_simplify_model=TRUE, trace=TRUE)
#' @export
optimize_model <- function(df, term, ..., autocorrelation_cols=NA,
                           automatic_removal=FALSE, autocorrelation_threshold=0.7, correlation_method='pearson', cor_use='complete.obs',
                           model_type='glm', model_family='gaussian', evaluation_methods=c('anova'), simplification_direction='both', omit.na='overall',
                           scale_predictor=FALSE, plot_quality_assessment='baseR', plot_relationships=FALSE, 
                           jitter_plots=FALSE, plot_type='boxplot', stat_test='wilcox',
                           round_p=5, backward_simplify_model=TRUE, trace=FALSE) {
  
  errors = c()
  # check for autocorrelations
  if (!all(is.na(autocorrelation_cols))) {
    are_cols_numeric = unlist(lapply(autocorrelation_cols, function(pr) is.numeric(df[,pr])))
    if (!all(are_cols_numeric)) {
      warning("There were non-numeric columns given as columns to check for autocorrelations. We're gonna remove them for you, but please don't do it again. ;)")
      autocorrelation_cols = autocorrelation_cols[are_cols_numeric]
    }
    
    autocorrelations_and_preds = remove_autocorrelations(df, autocorrelation_cols, use=cor_use,
                                                          automatic_removal = automatic_removal,
                                                          autocorrelation_threshold = autocorrelation_threshold,
                                                          correlation_method = correlation_method)
    autocorrelations = autocorrelations_and_preds$autocorrelations
  } else {
    autocorrelations = NA
  }
  
  # some basic checks
  if (((model_type == 'lm') || (model_type == 'lmer')) && (model_family != 'gaussian')) {
    warning("When opting for lm and lmer models, you should always specify gaussian as the model family. But don't ya worry, I'll change the family for you. ;)")
    model_family = 'gaussian'
  }
  if (!('anova' %in% evaluation_methods) && (grepl('quasi', model_family))) {
    warning("Anova is the only method (that we support) that works with quasibinomial distributions. Imma change it for ya real quick, but please remember this for your future endeavors. ;)")
    evaluation_methods = c('anova')
  }
  
  # check that interactions + main effects are in formula (i.e., x:y + x + y)
  main_effects = as.character(find_call(term, return='main_effects'))
  interactions = find_call(term, return='interactions')
  need_to_add = FALSE
  
  extract_main_effects = function(interaction) {
    interaction_w_tr = interaction[[1]] != ':'
    
    if (!interaction_w_tr) {
      return(str_split(deparse(interaction), ':')[[1]])
    } else {
      return(extract_main_effects(interaction[[2]]))
    }
  }
  
  term = as.formula(term)
  for (interaction in interactions) {
    effects = extract_main_effects(interaction)
    
    for (effect in effects) {
      new_term = paste('~ . +', effect)
      if (!(effect %in% main_effects)) {
        term = update(term, new_term)
        need_to_add = TRUE
      }
    }
  }
  
  if (need_to_add) warning("Please remember to add all main effects as separate items in your formula when including interactions. We've added all main effects for you this time. Might not next time, though, so better include 'em yourself. ;)")
  
  response_frm = term[[2]]
  predictor_frm = term[[3]]
  
  # remove autocorrelations from formula and gather categorical variables
  split_frm = find_call(term, return = 'all')
  
  categorical_vars = vector()
  numerical_vars = vector()
  remaining_pred_call = vector()
  for (frm_part in split_frm) {
    frm_str = deparse(frm_part)

    if ((frm_str %in% colnames(df)) && (is.factor(df[,frm_str]))) {
      categorical_vars = append(categorical_vars, frm_str)
    } else if ((frm_str %in% colnames(df)) && (is.numeric(df[,frm_str]))) {
      numerical_vars = append(numerical_vars, frm_str)
    }
    
    if (!all(is.na(autocorrelation_cols)) && (nrow(autocorrelations) > 0)) {
      ptrn = str_interp('(${paste(autocorrelations_and_preds$removed_predictor_cols, collapse="|")})(?![a-z0-9_\\-#])')
      if (!str_detect(frm_str, regex(ptrn, ignore_case = TRUE))) {
        remaining_pred_call = append(remaining_pred_call, frm_part)
      }
    }
  }
  
  if (!all(is.na(autocorrelation_cols)) && (nrow(autocorrelations) > 0)) {
    predictor_frm = remaining_pred_call[[1]]
    for (i in 2:length(remaining_pred_call)) { predictor_frm = call('+', predictor_frm, remaining_pred_call[[i]]) }
  }
  
  if (scale_predictor) {
    for (numerical_var in numerical_vars) {
      df[str_interp("${numerical_var}_aB3cD5eF6G")] = df[,numerical_var]
      df[numerical_var] = as.vector(scale(df[,numerical_var]))
    }
  }
  
  var_cnt = length(find_call(term, return = 'atomic'))
  if (nrow(df) / var_cnt <= 4) {
    warning(str_interp("There are probably too many variables (${var_cnt}) in comparison to datapoints (${nrow(df)}). Please add data, or remove x variables to reduce overfitting, multicollinearity, convergence issues, statistical significance and interpretability."))
  }
  
  term = substitute(response ~ predictor, list(response = response_frm, predictor = predictor_frm))
  
  if (model_type != 'nlmer') {
    possible_model_family = determine_model_family(df, response_frm)
    if (model_family == 'automatic') {
      if (length(possible_model_family) > 1) {
        print(str_interp("Your response variable allows for ${paste(possible_model_family, collapse=' and ')} distribution. Which one would you prefer?"))
        model_family = readline()
      } else if (length(possible_model_family) == 0) stop("No appropriate family found. Please check data types of columns.")
      else {
        model_family = possible_model_family[1]
        print(str_interp("Continuing with ${model_family} distribution."))
      }
    } else if (!(model_family %in% possible_model_family)) {
      warning(str_interp("Chosen distribution '${model_family}' does not match response values. Would recommend ${paste(possible_model_family, collapse=', or ')}. Please check."))
    }
  }
  
  simplified_models = simplify_model(df, model_type, term, evaluation_methods, ..., 
                                     model_family=model_family,
                                     categorical_vars=categorical_vars,
                                     direction=simplification_direction,
                                     backward_simplify_model=backward_simplify_model, 
                                     trace=trace, omit.na=omit.na)
  
  errors = simplified_models$errors
  for (error in unique(errors)) warning(error)
  models_with_info = list()
  
  for (direction in names(simplified_models)) {
    simplified_model_info = simplified_models[[direction]]
    model_histories = simplified_model_info$history
    final_model = simplified_model_info$final_model
    
    if (scale_predictor) {
      for (numerical_var in numerical_vars) {
        df[numerical_var] = df[,str_interp("${numerical_var}_aB3cD5eF6G")] 
        df[numerical_var] = df[numerical_var] %>%
          select(!c(!!numerical_var))
      }
    }
    
    if ((!is.list(final_model)) && (!typeof(final_model) == 'S4')) {
      simplified_models[[direction]] = simplified_model_info
      next
    }
    
    all_sign_vars = c(simplified_model_info$significant_variables, simplified_model_info$marginally_significant_variables)
    all_sign_vars = all_sign_vars[!is.na(all_sign_vars)]
    
    if (length(all_sign_vars) != 0) {
      # CHECK MODEL
      model_plots = list()
      stat_results = list()
      if (plot_quality_assessment == 'performance') {
        if (model_type == 'gam') plot(check_model(final_model, residual_type = 'normal'))
        else if (model_family == 'binomial') withCallingHandlers(plot(check_model(final_model, type = 'discrete_both')),
                                                                 warning = function(w) {
                                                                   if (grepl("stat_density", w$message)) {
                                                                     warning(str_interp("${w$message}. Please choose another model family."))
                                                                     tryInvokeRestart("muffleWarning")
                                                                   }
                                                                   else {
                                                                     warning(w$message)
                                                                   }
                                                                 })
        else plot(check_model(final_model))
        checked_model = recordPlot()
        dev.off()
      } else if (plot_quality_assessment == 'baseR') {
        par(mfrow = c(2, 2), mar = rep(2, 4))
        plot(final_model)
        checked_model = recordPlot()
        dev.off()
      } else {
        checked_model = NA
      }
      model_plots$model_check = checked_model
      categorical_check = setup_categorical_check(df, categorical_vars, quote(predictor))
      
      if (model_type == 'gam') {
        models_overview = summary(final_model)$p.table
        models_overview_s = summary(final_model)$s.table
        models_overview_s = slightly_expanded_model_summary(models_overview_s, response_frm, term, categorical_check, round_p, df)
      }
      else models_overview = coef(summary(final_model))
      
      models_overview = slightly_expanded_model_summary(models_overview, response_frm, term, categorical_check, round_p, df)
      models_overview = models_overview %>%
        mutate(effect_direction = case_when(
          Estimate < 0 ~ 'negative',
          TRUE ~ 'positive'
        ))
      
      if (plot_relationships && (nrow(models_overview) > 0)) {
        plot_out = plot_model_features(final_model, models_overview, model_type, model_family, plot_type, jitter_plots, stat_test, round_p, !backward_simplify_model)
        model_plots = append(model_plots, plot_out)
      }
      if (plot_relationships && (model_type == 'gam')) {
        plot_out = plot_model_features(final_model, models_overview_s, model_type, model_family, plot_type, jitter_plots, stat_test, round_p, !backward_simplify_model)
        model_plots = append(model_plots, plot_out)
      }
    } else {
      simplified_models[[direction]] = simplified_model_info
      next
    }
    
    models_overview = models_overview %>%
      rename(Response = response,
             Predictor = predictor,
             `Variable Type` = var_type,
             Significance = significance,
             `Effect Direction` = effect_direction) %>%
      select(-predictor_col)
    
    if (model_type == 'gam') models_overview = list(p.table = models_overview, s.table = models_overview_s)
    
    simplified_models[[direction]] = list(
      overview = models_overview,
      final_model = simplified_model_info$final_model,
      plots = model_plots
    )
    
    if (trace) simplified_models[[direction]]$history = simplified_model_info$history
  }
  
  return(list(
    autocorrelations = autocorrelations,
    models_with_info = simplified_models
  ))
}

#' Plot model features
#'
#' Plot estimate, regression and effect size
#' @param regression_model The (simplified) regression model
#' @param models_overview Custom model summary (see [slightly_expanded_model_summary()])
#' @param model_type Type of regression model (e.g., 'glm')
#' @param model_family Model family used for glm calculation (default: gaussian)
#' @param plot_type Whether to plot as boxplot or violin plot for categorical variables (default: boxplot)
#' @param jitter_plots Whether to use jitter when generating geom_point plots (default: FALSE)
#' @param test Test to be used for pairwise significance testing for categorical variables. Either wilcox or t.test (default: wilcox)
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param remove_insignificant Used to exclude insignificant relationships from both estimate and effect size plots and omit curves from regression plots (default: FALSE)
#' @return A list of plots (estimate, regression, and effect size) and statistics for categorical variables
#' @examples
#' data("plants");
#' 
#' model_type = 'glm'
#' model_family = 'quasibinomial'
#' categorical_check = setup_categorical_check(plants, c('habitat', 'ploidy'), col=quote(predictor))
#' final_model = glm(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + isothermality + habitat + ploidy + solar_radiation:isothermality, family=model_family, data=plants)
#' 
#' models_overview = coef(summary(final_model))
#' models_overview = slightly_expanded_model_summary(models_overview, quote(sexual_seed_prop), quote(sexual_seed_prop ~ altitude + solar_radiation + annual_mean_temperature + isothermality + I(isothermality^2) + habitat + ploidy + solar_radiation:annual_mean_temperature + solar_radiation:isothermality + annual_mean_temperature:isothermality), categorical_check, 3, plants)
#' models_overview = models_overview |>
#' dplyr::mutate(effect_direction = dplyr::case_when(
#' Estimate < 0 ~ 'negative',
#' TRUE ~ 'positive'
#' ))
#' 
#' plot_model_features(final_model, models_overview, model_type, model_family, plot_type='boxplot', jitter_plots=TRUE, test='wilcox', round_p=3, remove_insignificant=FALSE)
#' @export
plot_model_features <- function(regression_model, models_overview, model_type, model_family, plot_type='boxplot', jitter_plots=FALSE, test='wilcox', round_p=5, remove_insignificant=FALSE) {
  if (remove_insignificant) models_overview_ns = models_overview[models_overview$significance != 'ns',]
  else models_overview_ns = models_overview
  
  model_plots = list(
    regression_plots = plot_regression(regression_model, models_overview, model_type, model_family, plot_type, jitter_plots, test, round_p, !remove_insignificant)
  )
  
  if ('Estimate' %in% colnames(models_overview)) {
    model_plots[['estimate_plot']] = plot_estimate(models_overview_ns, regression_model, model_type)
    model_plots[['effect_size_plot']] = plot_effect_size(models_overview_ns)
  }
  
  return(model_plots)
}

#' Plot regression model estimate
#'
#' Plot estimate, regression and effect size
#' @param models_overview Custom model summary (see [slightly_expanded_model_summary()])
#' @param regression_model The (simplified) regression model
#' @param model_type Type of regression model (e.g., 'glm')
#' @return A plot of the estimates listed within the model summary
plot_estimate <- function(models_overview, regression_model, model_type) {
  if ((model_type == 'lmer') || (model_type == 'glmer')) {
    df = model.frame(regression_model)
  } else {
    df = regression_model$model
  }
  
  new_rows = lapply(unique(models_overview[models_overview['var_type'] == 'categorical', 'predictor_col']), function(cat_var) {
    cat_var_states = unique(df[[cat_var]])
    reference = setdiff(cat_var_states, sub(cat_var, '', models_overview$predictor))
    data.frame(
      predictor = paste0(cat_var, reference),
      predictor_col = cat_var,
      Estimate = 0,
      `Std. Error` = 0,
      significance = 'ref',
      stringsAsFactors = FALSE
    )
  })
  models_overview = bind_rows(models_overview, do.call(rbind, new_rows))
  models_overview = models_overview %>%
    arrange(desc(predictor_col), desc(predictor)) %>%
    mutate(formatted_predictor = case_when((significance != 'ref') & (var_type == 'categorical') ~ paste0('italic(', predictor, ')'),
                                           TRUE ~ predictor),
           formatted_predictor = factor(formatted_predictor, levels=formatted_predictor))
 
  formatted_labels = sapply(models_overview$formatted_predictor, function(x) {
    if (grepl("^italic", x)) {
      parse(text = as.character(x))
    } else {
      as.character(x)
    }
  })
  
  p = ggplot(models_overview, aes(x=Estimate, y=formatted_predictor)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(hjust = 0)) +
    labs(y='Predictor') +
    scale_y_discrete(labels = formatted_labels) +
    geom_vline(xintercept=0, colour="black", linewidth=.5) +
    geom_point(aes(color=Estimate)) +
    geom_errorbar(aes(xmin=Estimate-`Std. Error`, xmax=Estimate+`Std. Error`, color=Estimate), width=.0) +
    scale_color_continuous_diverging(rev=TRUE)
  return(p)
}

#' Plot model regression (categorical)
#'
#' Plot categorical regression of model variables
#' @param cat_var Categorical variable to plot
#' @param df Dataframe object extracted from model
#' @param models_overview Custom model summary (see [slightly_expanded_model_summary()])
#' @param test Test to be used for pairwise significance testing for categorical variables. Either wilcox or t.test
#' @param response_str Response as string
#' @param response_col Response column name
#' @param plot_type Whether to plot as boxplot or violin plot for categorical variables
#' @param regression_plots A list storing regression plots
#' @param stat_results A list storing results of statistics
#' @return Regression plots of variables listed within the model summary + statistics
plot_regression_categorical <- function(cat_var, df, models_overview, test, response_str, response_col, plot_type, regression_plots, stat_results) {
  models_overview_group = models_overview[models_overview$predictor_col == cat_var,] %>%
    mutate(factor_group = str_remove(predictor, predictor_col))
  
  new_rows = lapply(unique(models_overview_group[models_overview_group['var_type'] == 'categorical', 'predictor_col']), function(cat_var) {
    cat_var_states = unique(df[[cat_var]])
    reference = setdiff(cat_var_states, sub(cat_var, '', models_overview_group$predictor))
    data.frame(
      response = models_overview_group[1, 'response'],
      predictor = paste0(cat_var, reference),
      predictor_col = cat_var,
      Estimate = 0,
      significance = 'ref',
      stringsAsFactors = FALSE
    )
  })
  models_overview_group = bind_rows(models_overview_group, do.call(rbind, new_rows))
  
  df = merge(df, models_overview_group, by.x=cat_var, by.y='factor_group', all.x=TRUE, all.y=FALSE)
  
  cat_var = parse(text=cat_var)[[1]]
  response_sym = parse(text=response_str)[[1]]
  
  if (nrow(models_overview_group) > 1) {  # multiple vars
    stat_out = run_stats(df, response_col, cat_var, test=test)
    df = stat_out$df
    stat_results[[cat_var]] = stat_out$stat_results
    
    df = df %>% mutate(significance = ifelse(is.na(significance), 'ref', significance))
    p = ggplot(data=df, aes(x=!!cat_var, y=!!response_sym))
    if (plot_type == 'boxplot') p = p + geom_boxplot(aes(color=!!cat_var, fill=!!cat_var))
    else p = p + geom_violin(aes(color=!!cat_var, fill=!!cat_var))
    
    p = p +
      scale_color_viridis_d(option = "G") +
      scale_fill_viridis_d(option = "G") +
      geom_text(aes(x=!!cat_var, label=letter), y=max(df[!is.na(df[response_str]),response_str])*1.1, check_overlap = TRUE) +
      geom_text(aes(x=!!cat_var, label=significance), y=max(df[!is.na(df[response_str]),response_str]*1.05), check_overlap = TRUE) +
      ylim(min(df[!is.na(df[response_str]),response_str]), max(df[!is.na(df[response_str]),response_str])*1.1)
  } else {  # single vars
    p = ggplot(data=df, aes(x=!!cat_var, y=!!response_sym))
    if (plot_type == 'boxplot') p = p + geom_boxplot(color='#21918c')
    else p = p + geom_violin(color='#21918c')
  }
  
  p = p + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust=1))
  regression_plots[[cat_var]] = p
  
  return(list(regression_plots=regression_plots, stat_results=stat_results))
}

#' Plot model regression (numeric & mixed)
#'
#' Plot numeric/mixed regression of model variables
#' @param predictor_full Numeric/mixed variable to plot
#' @param regression_model The (simplified) regression model
#' @param df Dataframe object extracted from model
#' @param term The formula to be used with the model
#' @param models_overview Custom model summary (see [slightly_expanded_model_summary()])
#' @param response_str Response as string
#' @param regression_plots A list storing regression plots
#' @param jitter Whether to use jitter when generating dotplots
#' @param plot_curve Whether to plot regression plot curve using geom_smooth
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param model_family Model family used for glm calculation
#' @param model_type Type of regression model (e.g., 'glm')
#' @return Regression plots of variables listed within the model summary
plot_regression_numeric <- function(predictor_full, regression_model, df, term, models_overview, response_str, regression_plots, jitter, plot_curve, round_p, model_family, model_type) {
  if ('Estimate' %in% colnames(models_overview)) estimate = round(models_overview[models_overview['predictor'] == predictor_full, 'Estimate'], digits = round_p)
  else estimate = 'NA'
  sign = models_overview[models_overview['predictor'] == predictor_full, 'significance']
  
  # get predictor columns (i.e. either one or two depending on whether there's an interaction)
  predictor_cols = vector()
  for (predictor_col in colnames(df)) {
    if ('Estimate' %in% colnames(models_overview)) ptrn = str_interp('\\Q${predictor_col}\\E(?=$|\\:)')  # TODO/test pattern may need more work
    else ptrn = str_interp('\\(?\\Q${predictor_col}\\E\\)?(?=$|\\:)')  # TODO: check if this can be used in all cases
    
    if (str_detect(predictor_full, regex(ptrn))) {
      predictor_cols = append(predictor_cols, predictor_col)
    }
  }
  
  # check if columns are numeric
  cols_are_numeric = sapply(predictor_cols, function(pr) is.numeric(df[,pr]))
  numeric_col_exists = any(cols_are_numeric)
  
  # transform df variables according to formula
  if ((length(predictor_cols) > 2) && (!all(cols_are_numeric))) {
    warning(str_interp("We don't allow for plotting of >2 interactions for categorical variables. Skipping plotting for ${predictor_full}."))
    return()
  } else if (length(predictor_cols) == 2) {  # interaction
    # === split formula call to get interaction
    complete_call = find_call(term, predictor_full, return='match')[[1]]
    interaction_transformed = complete_call[[1]] != ':'
    
    # === transform individual variables
    preds_w_tr <- vector()
    for (i in 2:3) {  # 2 and 3 are the positions in the call :)
      pred_w_tr = complete_call[[i]]
      df[deparse(pred_w_tr)] = with(df, eval(pred_w_tr))
      preds_w_tr[[length(preds_w_tr)+1]] = deparse(pred_w_tr)
    }
    
    # === compute interaction if continuous
    if (all(cols_are_numeric)) {  # only continuous variables
      col_name = if (interaction_transformed) {str_interp('${preds_w_tr[[1]]}_${preds_w_tr[[2]]}')} else {predictor_full}
      df[col_name] = df[preds_w_tr[[1]]] * df[preds_w_tr[[2]]]
    }
    
    # === check if transform needs to be applied for interaction
    if (all(cols_are_numeric) && interaction_transformed) {  # transform around interaction
      interaction_call = deparse(complete_call[[2]])
      tr_inter_str = gsub(interaction_call, col_name, deparse(complete_call, width.cutoff=500))
      df[predictor_full] = with(df, eval(parse(text=tr_inter_str)))
    } else if (length(complete_call) != 3) {  # 3 == no transform around interaction, different value needs to be checked
      warning(str_interp("Internal error. Check interaction call. No plot generated for ${predictor_full}"))
    }
  } 
  else if (str_detect(predictor_full, regex('s\\(.+\\)'))) {  # no interaction + s() transform
    df[str_interp("${predictor_cols[1]}_smooth")] = predict(regression_model, type = "terms")[, predictor_full]
    predictor_full = str_interp("${predictor_cols[1]}_smooth")
  }
  
  if ((length(predictor_cols) == 1) | (all(cols_are_numeric))) {  # either one variable or only continuous variables
    if (length(predictor_cols) == 1) x = predictor_full
    else x = col_name
    
    x = parse(text=x)[[1]]
    response_str = parse(text=response_str)[[1]]
    
    if (str_detect(deparse(x), regex('I\\(.+\\^2\\)'))) {
      x_base = str_extract(deparse(x), regex('(?<=I\\().+(?=\\^2\\))'))
      new_x = str_interp("${x_base}_squared_${format(Sys.time(), format='%Y%m%d%H%M')}")
      df = df %>%
        rename(!!new_x := deparse(x))
      x = parse(text = new_x)[[1]]
      df[[new_x]] = as.numeric(df[[new_x]])
    } else if (str_detect(deparse(x), ':')) {
      new_x = gsub(':', '_', deparse(x))
      df = df %>%
        rename(!!new_x := deparse(x))
      x = parse(text = new_x)[[1]]
      df[[new_x]] = as.numeric(df[[new_x]])
    }
    
    p = ggplot(data=df, aes(x=!!x, y=!!response_str)) +
      geom_point(color='white', fill='#21918c', pch=21, position=jitter)
    
    if (plot_curve && sign != 'ns') p = p + geom_smooth(method = model_type, method.args = list(family=model_family), formula = y ~ x, se = FALSE, color='#21918c', fill='#21918c') # formula = y ~ x
  } else if (numeric_col_exists) {  # >1 variable with one continuous
    x = preds_w_tr[sapply(preds_w_tr, function(pr) is.numeric(df[,pr]))]
    sign_vars_sorted = c(preds_w_tr[x], preds_w_tr[!x])
    
    x = parse(text=sign_vars_sorted[1])[[1]]
    color = parse(text=sign_vars_sorted[2])[[1]]
    
    p = ggplot(data=df, aes(x=!!x, y=!!response_str)) +
      geom_point(aes(color=!!color), position=!!jitter) +
      scale_color_viridis_d(option = "G") +
      scale_fill_viridis_d(option = "G")
    
    if (plot_curve && sign != 'ns') p = p + geom_smooth(method = model_type, formula = y ~ x, method.args = list(family=model_family), 
                                                se = TRUE, aes(color=!!color, fill=!!color))
  }
  
  x_max = max(df[[deparse(x)]])
  x_min = min(df[[deparse(x)]])
  est_len = length(as.character(estimate))
  
  if (estimate > 0 | is.na(estimate)) {
    x = x_min + ((x_max - x_min) * (.1 + est_len * .01))
  } else {
    x = x_max - ((x_max - x_min) * (.1 + est_len * .01))
  }
  
  p = p +
    theme_minimal() +
    theme(legend.position="bottom") +
    guides(fill="none") +
    annotate("label", x=x, y = Inf, label = str_interp("estimate=${estimate}\n${sign}"), vjust=1, hjust=.5)
  regression_plots[[predictor_full]] = p
  
  return(regression_plots)
}

#' Plot model regression
#'
#' Plot regression of model variables
#' @param regression_model The (simplified) regression model
#' @param models_overview Custom model summary (see [expanded_model_summary()])
#' @param model_type Type of regression model (e.g., 'glm')
#' @param model_family Model family used for glm calculation (default: gaussian)
#' @param plot_type Whether to plot as boxplot or violin plot for categorical variables (default: boxplot)
#' @param jitter_plots Whether to use jitter when generating dotplots (default: FALSE)
#' @param test Test to be used for pairwise significance testing for categorical variables. Either wilcox or t.test (default: wilcox)
#' @param round_p Convenience parameter for automatic rounding of p-values
#' @param plot_curve Whether to plot regression plot curve using geom_smooth
#' @return Regression plots of variables listed within the model summary
plot_regression <- function(regression_model, models_overview, model_type, model_family, plot_type='boxplot', jitter_plots=FALSE, test='wilcox', round_p=5, plot_curve=FALSE) {
  if (jitter_plots) { jitter = 'jitter' } else { jitter = 'identity' }
  if ((model_type == 'lmer') || (model_type == 'glmer')) {
    df = model.frame(regression_model)
  } else {
    df = regression_model$model # TEST with lm
  }
  
  term = formula(regression_model)
  response_col = parse(text=all.vars(term)[1])[[1]]
  # all_predictor_cols <- all.vars(term)[c(-1)]
  predictors_frm = term[[3]]
  
  regression_plots = list()
  stat_results = list()
  
  # transform response
  # TEST: probably not needed when using regression_model$model/model.frame
  # TEST: won't work if we have 3-way formula... possibility of 3-way formula?
  response_frm = term[[2]]
  response_str = paste(response_frm)
  df[response_str] = with(df, eval(parse(text=response_frm)))
  
  # plot categorical vars (if no interaction with numeric)
  for (cat_var in unique(models_overview[models_overview$var_type == 'categorical', 'predictor_col'])) {
    res = plot_regression_categorical(cat_var, df, models_overview, test, response_str, response_col, plot_type, regression_plots, stat_results)
    regression_plots = res$regression_plots
    stat_results = res$stat_results
  }
  
  # for continuous or mixed variables
  for (predictor_full in models_overview[models_overview$var_type == 'numeric', 'predictor']) {
    regression_plots = plot_regression_numeric(predictor_full, regression_model, df, term, models_overview, response_str, regression_plots, jitter, plot_curve, round_p, model_family, model_type)
  }
  
  return(list(regression_plots=regression_plots, stat_results=stat_results))
}

#' Plot regression model effect sizes
#'
#' Plot effect sizes of model variables
#' @param models_overview Custom model summary (see function [expanded_model_summary()])
#' @return A plot of the effect sizes of the variables listed within the model summary
plot_effect_size = function(models_overview) {
  models_overview = models_overview %>%
    mutate(Estimate_abs = abs(Estimate),
           Est_sum = sum(Estimate_abs),
           `Effect size` = case_when(
             Estimate > 0 ~ (Estimate_abs / Est_sum) * 100,
             Estimate < 0 ~ (Estimate_abs / Est_sum) * -100
           )) %>%
    arrange(desc(`Effect size`))
  
  p = ggplot(models_overview, aes(y=`Effect size`, x=predictor, fill=`Effect size`)) +
    geom_bar(stat = 'identity') +
    scale_fill_continuous_diverging(rev=TRUE) +
    theme_minimal() +
    ylab('Effect size [%]') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  return(p)
}

#' Function for calculating statistics for given categorical variable
#'
#' Calculate pairwise statistics for given categorical variable and add letters indicating significant differences
#' @param df Dataframe with response and predictors as columns
#' @param response Character string representing the response/y variable (column name)
#' @param predictor Character string representing the predictor/x variable (column name)
#' @param test Test to be used for pairwise significance testing for categorical variables. Either wilcox or t.test.
#' @return Output of the statistical test and the updated dataframe with letters indicating significant differences between traits of a categorical variable
run_stats <- function(df, response, predictor, test='wilcox') {
  add_letters <- function(df, stat_results, test) {
    df_sorted = df %>%
      group_by(!!predictor) %>%
      mutate(resp_mean = mean(!!response)) %>%
      ungroup() %>%
      distinct(!!predictor, .keep_all = TRUE) %>%
      arrange(desc(!!response)) %>%
      mutate(idx = row_number())
  
    if (any(stat_results$p_value < 0.05)) {
      stat_results = merge(stat_results, df_sorted[c(deparse(predictor), 'idx')], by.x='var1', by.y=deparse(predictor))
      stat_results = merge(stat_results, df_sorted[c(deparse(predictor), 'idx')], by.x='var2', by.y=deparse(predictor), 
                            suffixes = c('1', '2'))
      # for accessibility
      stat_results = stat_results %>%
        rowwise() %>%
        mutate(comb = paste(sort(c(idx1, idx2)), collapse='|'))
      
      letter = 'A'
      for (i in 1:nrow(df_sorted)) {
        if (i == 1) { # highest accuracy
          df_sorted[1, 'letter'] = letter
        } else if (i != nrow(df_sorted)) { # everything between highest and lowest accuracy
          # capital letter in comment is the model we're assigning the letter to at the moment ;)
          if (stat_results[stat_results$comb == paste(i-1, i, sep='|'),]$p_value < 0.05) {
            # i and i-1 are significantly different (a-B lettering)
            new_letter = LETTERS[match(letter, LETTERS)+1]
            df_sorted[i, 'letter'] = new_letter
            letter <- new_letter
          } else if (stat_results[stat_results$comb == paste(i, i+1, sep='|'),]$p_value < 0.05) {
            # a-A-b
            df_sorted[i, 'letter'] = letter
          } else if (stat_results[stat_results$comb == paste(i-1, i+1, sep='|'),]$p_value < 0.05) {
            # a-AB-b
            new_letter <- LETTERS[match(letter, LETTERS)+1]
            df_sorted[i, 'letter'] = paste(letter, new_letter, sep='')
            letter <- new_letter
          } else {
            # a-A-a
            df_sorted[i, 'letter'] = letter
          }
        } else { # lowest accuracy
          if (stat_results[stat_results$comb == paste(i-1, i, sep='|'),]$p_value < 0.05) {
            # i and i-1 are significantly different (a-B lettering)
            df_sorted[i, 'letter'] = LETTERS[match(letter, LETTERS)+1]
          } else {
            # a-A-end
            df_sorted[i, 'letter'] = letter
          }
        }
      }
      # add letters to results dataframe
      df = merge(df, df_sorted[c(deparse(predictor), 'letter')], by=deparse(predictor))
    } else {
      df['letter'] = 'A'
    }
    return(df)
  }
  
  if (test == 'wilcox') {
    stat_result = withCallingHandlers(as.data.frame(pairwise.wilcox.test(df[[response]], df[[predictor]])$p.value),
                        warning = function(w) {
                          if (grepl("cannot compute exact p-value with ties", w$message)) {
                            tryInvokeRestart("muffleWarning")
                            as.data.frame(pairwise.wilcox.test(df[[response]], df[[predictor]], exact=FALSE)$p.value)
                          }
                          else warning(w$message)
                        })
    stat_results = stat_result %>%
      rownames_to_column(var='var1') %>%
      pivot_longer(cols = -c('var1'),
                   names_to = 'var2',
                   values_to = 'p_value') %>%
      filter((var1 != var2) & (!is.na(p_value)))
  } else if (test == 't.test') {
    stat_results = as.data.frame(pairwise.t.test(df[[response]], df[[predictor]])$p.value) %>%
      rownames_to_column(var='var1') %>%
      pivot_longer(cols = -c('var1'),
                   names_to = 'var2',
                   values_to = 'p_value') %>%
      filter((var1 != var2) & (!is.na(p_value)))
  } else {
    warning("We only allow wilcox and t.test for statistical testing.")
  }
  
  df = add_letters(df[,c(deparse(response), deparse(predictor), 'significance')], stat_results, test)
  return(list(df=df, stat_results=stat_results))
}