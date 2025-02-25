% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LazyModeler.R
\name{not_so_expanded_model_summary}
\alias{not_so_expanded_model_summary}
\title{A minor version of expanded model summary}
\usage{
not_so_expanded_model_summary(
  model_summary,
  response_frm,
  term,
  categorical_check,
  round_p,
  df
)
}
\arguments{
\item{model_summary}{Initial model summary}

\item{response_frm}{Response formula}

\item{term}{The formula to be used with the model. Can be either quote() or formula()}

\item{categorical_check}{Output of setup_categorical_check()}

\item{round_p}{Convenience parameter for automatic rounding of p-values}

\item{df}{Dataframe with response and predictors as columns}
}
\value{
Expanded summary of a regression model with added information on base variables, variable types and significance levels
}
\description{
Expand model summary by adding column names and information on interaction groups
}
