simplify_model <- function(
    formula,
    family,
    data,
    model_type,
    model_args,
    direction) {
  
  switch(
    direction,
    "backward" = simplify_backward(formula,
                                   family,
                                   data,
                                   model_type,
                                   model_args
                                   ),
    "forward" = simplify_forward(),
    "both" = {
      list(simplify_backward(formula,
                             family,
                             data,
                             model_type,
                             model_args),
           simplify_forward())
    }
  )
}

simplify_forward <- function() {
  TRUE
}

simplify_backward <- function(
    formula,
    family,
    data,
    model_type,
    model_args) {
  
  # CREATE MODEL
  regression_model <- create_model(formula,
                                   family,
                                   data,
                                   model_type,
                                   model_args)
}