#' Create Dataset for Parametric Survival Analysis
#'
#' This function reads ADaM datasets and creates a combined dataset suitable for
#' parametric survival analysis. It allows flexible specification of analysis
#' population, survival parameter, and stratification factors.
#'
#' @param adsl_path Character string specifying the file path to the ADSL dataset
#' @param adtte_path Character string specifying the file path to the ADTTE dataset
#' @param population Character string specifying the analysis population flag.
#'   Must be one of "SAFFL", "ITTFL", or "RANDFL". Default is "ITTFL"
#' @param variable Character string specifying the parameter code from PARAMCD
#'   column in ADTTE dataset. Default is "OS" (Overall Survival)
#' @param stratify_by Character string specifying the column name for stratification
#'   factor from ADSL dataset, or NULL for no stratification. Default is "SEX"
#'
#' @return A data frame containing the following columns:
#'   \item{SUBJID}{Subject identifier}
#'   \item{ARM}{Treatment arm}
#'   \item{STRATIFY}{Stratification factor (renamed from specified column), or "Overall" if no stratification}
#'   \item{SURVTIME}{Survival time}
#'   \item{CNSR}{Censoring indicator (1=censored, 0=event)}
#'   \item{EVENT}{Event indicator (1=event, 0=censored)}
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Reads SAS datasets using haven::read_sas()
#'   \item Filters data based on specified population flag and parameter code
#'   \item Creates an EVENT variable as 1 - CNSR for survival analysis
#'   \item Converts ARM to factor variable
#'   \item If stratify_by is specified, converts stratification factor to factor variable
#'   \item If stratify_by is NULL, creates a single "Overall" group for all subjects
#' }
#'
#' @importFrom haven read_sas
#' @importFrom dplyr select left_join filter mutate rename
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Create dataset for overall survival analysis with ITT population
#' surv_data <- DataParametricSurv(
#'   adsl_path = "path/to/adsl.sas7bdat",
#'   adtte_path = "path/to/adtte.sas7bdat",
#'   population = "ITTFL",
#'   variable = "OS",
#'   stratify_by = "SEX"
#' )
#'
#' # Example of using the stratification variable name in FitSurvMods
#' surv_data_region <- DataParametricSurv(
#'   adsl_path = "path/to/adsl.sas7bdat",
#'   adtte_path = "path/to/adtte.sas7bdat",
#'   population = "SAFFL",
#'   variable = "PFS",
#'   stratify_by = "REGION"
#' )
#'
#' # Use the original variable name in the output table
#' results <- FitSurvMods(surv_data_region, stratify_name = "REGION")
#' }
DataParametricSurv <- function(adsl_path,
                               adtte_path,
                               population = "ITTFL",
                               variable = "OS",
                               stratify_by = "SEX") {

  # Validate population parameter
  valid_populations <- c("SAFFL", "ITTFL", "RANDFL")
  if (!population %in% valid_populations) {
    stop("Population must be one of: ", paste(valid_populations, collapse = ", "))
  }

  # Read ADaM datasets
  adam_adsl <- haven::read_sas(adsl_path)
  s_adam_adtte <- haven::read_sas(adtte_path)

  # Check if stratify_by column exists in ADSL (only if not NULL)
  if (!is.null(stratify_by) && !stratify_by %in% names(adam_adsl)) {
    stop("Stratification variable '", stratify_by, "' not found in ADSL dataset")
  }

  # Check if population flag exists in ADTTE
  if (!population %in% names(s_adam_adtte)) {
    stop("Population flag '", population, "' not found in ADTTE dataset")
  }

  # Check if variable exists in PARAMCD
  if (!variable %in% s_adam_adtte$PARAMCD) {
    stop("Variable '", variable, "' not found in PARAMCD column of ADTTE dataset")
  }

  # Define variables for NSE to avoid R CMD check notes
  SUBJID <- ARM <- PARAMCD <- CNSR <- EVENT <- STRATIFY <- NULL
  SURVTIME <- NULL

  # Select columns based on whether stratification is used
  if (is.null(stratify_by)) {
    # No stratification - select only SUBJID and ARM
    adsl_selected <- adam_adsl %>%
      dplyr::select(SUBJID, ARM)
  } else {
    # With stratification - select SUBJID, ARM, and stratification variable
    adsl_selected <- adam_adsl %>%
      dplyr::select(SUBJID, ARM, !!rlang::sym(stratify_by))
  }

  # Create dataset including necessary information
  dataset <- adsl_selected %>%
    dplyr::left_join(
      s_adam_adtte %>%
        dplyr::select(SUBJID, PARAMCD, SURVTIME, CNSR, !!rlang::sym(population)),
      by = 'SUBJID'
    ) %>%
    dplyr::filter(
      PARAMCD == variable &
        !!rlang::sym(population) == 'Y'
    ) %>%
    dplyr::mutate(
      ARM = as.factor(ARM),
      EVENT = 1 - CNSR
    )

  # Handle stratification variable
  if (is.null(stratify_by)) {
    # Create overall group when no stratification
    dataset <- dataset %>%
      dplyr::mutate(STRATIFY = factor("Overall")) %>%
      dplyr::select(SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT)
  } else {
    # Use specified stratification variable
    dataset <- dataset %>%
      dplyr::mutate(!!rlang::sym(stratify_by) := as.factor(!!rlang::sym(stratify_by))) %>%
      dplyr::rename(STRATIFY = !!rlang::sym(stratify_by)) %>%
      dplyr::select(SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT)
  }

  return(dataset)
}
