test_that("dataset info model optimization without errors", {
  data("dataset_info")
  coefficients <- c('val_min_gene_length', 'val_max_gene_length',
                    'train_min_gene_length', 'train_max_gene_length',
                    'train_num_samples', 'val_num_samples',
                    'mean_gene_length', 'val_mean_gene_length',
                    'train_mean_gene_length')
  resp_preds <- list('genus_confusion'=coefficients, 
                     'species_confusion'=coefficients[!coefficients == 'val_num_samples'])
  expect_no_error(optimize_model(dataset_info, resp_preds,
                                    glm_family='quasibinomial',
                                    automatic_removal=TRUE))
  
  expect_no_error(optimize_model(dataset_info, resp_preds,
                                 glm_family='quasibinomial',
                                 automatic_removal=TRUE, plot=TRUE))
})