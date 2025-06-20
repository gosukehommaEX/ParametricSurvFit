% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractCovarianceMatrix.R
\name{ExtractCovarianceMatrix}
\alias{ExtractCovarianceMatrix}
\title{Extract Covariance Matrix from Single Parametric Survival Model}
\usage{
ExtractCovarianceMatrix(
  dataset,
  distribution,
  parameter_name = "Overall Survival",
  stratify_reference_level = NULL,
  include_shape = FALSE,
  format_output = TRUE,
  create_kable = FALSE
)
}
\arguments{
\item{dataset}{A data frame created by DataParametricSurv() function containing
survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT}

\item{distribution}{Character string specifying which distribution to fit.
Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"}

\item{parameter_name}{Character string for the parameter name in output table.
Default is "Overall Survival"}

\item{stratify_reference_level}{Character string specifying which level of STRATIFY
should be used as reference (coded as 1). If NULL (default), uses the first level}

\item{include_shape}{Logical indicating whether to include shape parameters
in the covariance matrix. Only applicable for distributions with shape parameters
(gamma, gengamma, gompertz). Default is FALSE}

\item{format_output}{Logical indicating whether to return formatted output.
If TRUE (default), returns formatted data frame. If FALSE, returns raw matrix}

\item{create_kable}{Logical indicating whether to create formatted kable output.
If TRUE, returns kableExtra formatted table with merged cells. If FALSE (default), returns data frame}
}
\value{
A data frame or matrix containing covariance matrix information
}
\description{
This function extracts the variance-covariance matrix from a single fitted
parametric survival model with explicit control over shape parameter inclusion.
}
\details{
The function handles different distribution types:
\itemize{
  \item \strong{No Shape Parameters}: exp, weibull, lnorm, llogis
  \item \strong{Has Shape Parameters}: gamma (shape), gengamma (Q), gompertz (shape)
}

For distributions with shape parameters, use \code{include_shape} to control
whether the shape parameter is included in the covariance matrix calculation.

When \code{create_kable = TRUE}, the output includes merged cells for
Distribution and Parameter columns, similar to FitSurvMods formatting.
}
\examples{
\dontrun{
# Exponential distribution (no shape parameter)
cov_exp <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "exp"
)

# Weibull distribution with formatted kable output
cov_weibull_kable <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "weibull",
  create_kable = TRUE
)

# Generalized Gamma without shape in covariance matrix
cov_gengamma_no_shape <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "gengamma",
  include_shape = FALSE
)

# Generalized Gamma with shape in covariance matrix
cov_gengamma_with_shape <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "gengamma",
  include_shape = TRUE
)

# Custom stratification reference level with kable
cov_custom_ref <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "weibull",
  stratify_reference_level = "Male",  # For SEX stratification
  create_kable = TRUE
)
}
}
