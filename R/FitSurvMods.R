#' Fit Multiple Parametric Survival Models and Create Formatted Output
#'
#' This function fits multiple parametric survival distributions to survival data
#' and creates a formatted table output using kableExtra package. It performs
#' stratified analysis by treatment arm and specified stratification factors.
#'
#' @param dataset A data frame created by DataParametricSurv() function containing
#'   survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT
#' @param distributions Character vector specifying which distributions to fit.
#'   Default includes all available distributions: "exp", "weibull", "lnorm",
#'   "llogis", "gompertz", "gengamma", "gamma"
#' @param table_caption Character string for table caption. Default is
#'   "Parametric Survival Model Results"
#' @param format_output Logical indicating whether to return formatted kable output.
#'   If TRUE (default), returns formatted table. If FALSE, returns raw data frame
#' @param stratify_name Character string specifying the original stratification
#'   variable name for column header. If NULL (default), uses "STRATIFY"
#'
#' @return If format_output = TRUE: A formatted kable object ready for display.
#'   If format_output = FALSE: A data frame with the following columns:
#' \describe{
#'   \item{Distribution}{Distribution name}
#'   \item{ARM}{Treatment arm}
#'   \item{STRATIFY}{Stratification factor level}
#'   \item{Parameter}{Parameter name}
#'   \item{Coefficient Estimate}{Parameter estimate}
#'   \item{SE}{Standard error}
#'   \item{95% Confidence Limits}{Formatted confidence interval}
#' }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Fits specified parametric distributions using flexsurv::flexsurvreg()
#'   \item Extracts parameters using ExtractParams() function
#'   \item Combines results across all distributions and strata
#'   \item Creates formatted output table using kableExtra
#' }
#'
#' The formatted table includes:
#' \itemize{
#'   \item Grouped headers by distribution
#'   \item Alternating row colors for readability
#'   \item Centered column alignment
#'   \item Hover effects for interactive viewing
#' }
#'
#' @importFrom dplyr group_by reframe mutate select arrange bind_rows mutate_if
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv
#' @importFrom kableExtra kable kable_styling collapse_rows column_spec
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Create survival dataset first
#' surv_data <- DataParametricSurv(
#'   adslph_path = "path/to/adslph.sas7bdat",
#'   adtte_path = "path/to/adtte.sas7bdat",
#'   population = "ITTFL",
#'   variable = "OS",
#'   stratify_by = "SEX"
#' )
#'
#' # Fit all available distributions with formatted output
#' results_table <- FitSurvMods(surv_data)
#' print(results_table)
#'
#' # Fit specific distributions only
#' selected_results <- FitSurvMods(
#'   dataset = surv_data,
#'   distributions = c("exp", "weibull", "lnorm"),
#'   format_output = TRUE,
#'   table_caption = "Selected Parametric Models for Overall Survival"
#' )
#' print(selected_results)
#'
#' # Fit models with custom stratification column name
#' results_custom <- FitSurvMods(
#'   dataset = surv_data,
#'   distributions = c("exp", "weibull", "gamma"),
#'   table_caption = "Parametric Models by SEX",
#'   stratify_name = "SEX"
#' )
#'
#' # Fit models for PFS data
#' pfs_data <- DataParametricSurv(
#'   adslph_path = "path/to/adslph.sas7bdat",
#'   adtte_path = "path/to/adtte.sas7bdat",
#'   population = "ITTFL",
#'   variable = "PFS",
#'   stratify_by = "REGION"
#' )
#'
#' pfs_results <- FitSurvMods(
#'   dataset = pfs_data,
#'   distributions = c("exp", "weibull", "gamma"),
#'   table_caption = "Parametric Models for Progression-Free Survival"
#' )
#' }
FitSurvMods <- function(dataset,
                        distributions = c('exp', 'weibull', 'lnorm', 'llogis', 'gompertz', 'gengamma', 'gamma'),
                        format_output = TRUE,
                        table_caption = "Parametric Survival Model Results",
                        stratify_name = NULL) {

  # Validate input dataset
  required_cols <- c("SUBJID", "ARM", "STRATIFY", "SURVTIME", "CNSR", "EVENT")
  missing_cols <- setdiff(required_cols, names(dataset))
  if (length(missing_cols) > 0) {
    stop("Dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate distributions
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  invalid_distributions <- setdiff(distributions, valid_distributions)
  if (length(invalid_distributions) > 0) {
    stop("Invalid distributions specified: ", paste(invalid_distributions, collapse = ", "),
         "\nValid options are: ", paste(valid_distributions, collapse = ", "))
  }

  # Define variables for NSE to avoid R CMD check notes
  ARM <- STRATIFY <- fit_model <- param_data <- Parameters <- NULL
  Estimates <- SEs <- LCLs <- UCLs <- Distribution <- Parameter <- NULL

  # Initialize empty list to store results
  all_results <- list()

  # Loop through each distribution
  for (dist in distributions) {

    # Fit model for each group
    result_dist <- dataset %>%
      dplyr::group_by(ARM, STRATIFY) %>%
      dplyr::reframe(
        # Model information
        Distribution = dist,
        # Fit model
        fit_model = list(
          flexsurv::flexsurvreg(survival::Surv(SURVTIME, EVENT) ~ 1, dist = dist)
        ),
        # Extract parameters using the helper function
        param_data = list(ExtractParams(fit_model[[1]], dist))
      ) %>%
      # Unnest parameter data
      dplyr::mutate(
        Parameters = purrr::map(param_data, ~ .x[['Parameters']]),
        Estimates  = purrr::map(param_data, ~ .x[['Estimates']]),
        SEs  = purrr::map(param_data, ~ .x[['SEs']]),
        LCLs = purrr::map(param_data, ~ .x[['LCLs']]),
        UCLs = purrr::map(param_data, ~ .x[['UCLs']])
      ) %>%
      dplyr::select(-fit_model, -param_data) %>%
      # Convert to long format
      tidyr::unnest(c(Parameters, Estimates, SEs, LCLs, UCLs)) %>%
      # Create confidence interval text
      dplyr::mutate(
        `95% Confidence Limits` = paste0(round(LCLs, 2), ', ', round(UCLs, 2))
      ) %>%
      # Select and rename final columns
      dplyr::select(
        Distribution, ARM, STRATIFY,
        Parameter = Parameters,
        `Coefficient Estimate` = Estimates,
        SE = SEs,
        `95% Confidence Limits`
      )

    # Add to results list
    all_results[[dist]] <- result_dist
  }

  # Combine all results
  final_result <- dplyr::bind_rows(all_results) %>%
    dplyr::arrange(Distribution, ARM, STRATIFY, Parameter) %>%
    dplyr::mutate_if(is.numeric, round, 3)

  # Set column name for stratification variable
  stratify_col_name <- if (is.null(stratify_name)) "STRATIFY" else stratify_name
  colnames(final_result)[colnames(final_result) == "STRATIFY"] <- stratify_col_name

  # Return raw data frame if format_output is FALSE
  if (!format_output) {
    return(final_result)
  }

  # Create formatted kable output
  formatted_table <- final_result %>%
    kableExtra::kable(
      caption = table_caption,
      align = c('l', 'c', 'c', 'l', 'c', 'c', 'c'),
      digits = 3
    ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "center"
    ) %>%
    kableExtra::column_spec(1, bold = TRUE, color = "darkblue") %>%
    kableExtra::column_spec(2:3, background = "#f7f7f7")

  # Add cell merging for repeated values using collapse_rows
  formatted_table <- formatted_table %>%
    kableExtra::collapse_rows(columns = 1:3, valign = "middle")

  return(formatted_table)
}
