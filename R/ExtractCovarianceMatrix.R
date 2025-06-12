#' Extract Covariance Matrix from Single Parametric Survival Model
#'
#' This function extracts the variance-covariance matrix from a single fitted
#' parametric survival model with explicit control over shape parameter inclusion.
#'
#' @param dataset A data frame created by DataParametricSurv() function containing
#'   survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT
#' @param distribution Character string specifying which distribution to fit.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param parameter_name Character string for the parameter name in output table.
#'   Default is "Overall Survival"
#' @param stratify_reference_level Character string specifying which level of STRATIFY
#'   should be used as reference (coded as 1). If NULL (default), uses the first level
#' @param include_shape Logical indicating whether to include shape parameters
#'   in the covariance matrix. Only applicable for distributions with shape parameters
#'   (gamma, gengamma, gompertz). Default is FALSE
#' @param format_output Logical indicating whether to return formatted output.
#'   If TRUE (default), returns formatted data frame. If FALSE, returns raw matrix
#' @param create_kable Logical indicating whether to create formatted kable output.
#'   If TRUE, returns kableExtra formatted table with merged cells. If FALSE (default), returns data frame
#'
#' @return A data frame or matrix containing covariance matrix information
#'
#' @importFrom survival survreg Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom dplyr mutate arrange
#' @importFrom stats vcov
#' @importFrom kableExtra kable kable_styling column_spec collapse_rows
#' @importFrom magrittr %>%
#' @export
#'
#' @details
#' The function handles different distribution types:
#' \itemize{
#'   \item \strong{No Shape Parameters}: exp, weibull, lnorm, llogis
#'   \item \strong{Has Shape Parameters}: gamma (shape), gengamma (Q), gompertz (shape)
#' }
#'
#' For distributions with shape parameters, use \code{include_shape} to control
#' whether the shape parameter is included in the covariance matrix calculation.
#'
#' When \code{create_kable = TRUE}, the output includes merged cells for
#' Distribution and Parameter columns, similar to FitSurvMods formatting.
#'
#' @examples
#' \dontrun{
#' # Exponential distribution (no shape parameter)
#' cov_exp <- ExtractCovarianceMatrix(
#'   dataset = surv_data,
#'   distribution = "exp"
#' )
#'
#' # Weibull distribution with formatted kable output
#' cov_weibull_kable <- ExtractCovarianceMatrix(
#'   dataset = surv_data,
#'   distribution = "weibull",
#'   create_kable = TRUE
#' )
#'
#' # Generalized Gamma without shape in covariance matrix
#' cov_gengamma_no_shape <- ExtractCovarianceMatrix(
#'   dataset = surv_data,
#'   distribution = "gengamma",
#'   include_shape = FALSE
#' )
#'
#' # Generalized Gamma with shape in covariance matrix
#' cov_gengamma_with_shape <- ExtractCovarianceMatrix(
#'   dataset = surv_data,
#'   distribution = "gengamma",
#'   include_shape = TRUE
#' )
#'
#' # Custom stratification reference level with kable
#' cov_custom_ref <- ExtractCovarianceMatrix(
#'   dataset = surv_data,
#'   distribution = "weibull",
#'   stratify_reference_level = "Male",  # For SEX stratification
#'   create_kable = TRUE
#' )
#' }
ExtractCovarianceMatrix <- function(dataset,
                                    distribution,
                                    parameter_name = "Overall Survival",
                                    stratify_reference_level = NULL,
                                    include_shape = FALSE,
                                    format_output = TRUE,
                                    create_kable = FALSE) {

  # Validate input dataset
  required_cols <- c("SUBJID", "ARM", "STRATIFY", "SURVTIME", "CNSR", "EVENT")
  missing_cols <- setdiff(required_cols, names(dataset))
  if (length(missing_cols) > 0) {
    stop("Dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate distribution (single distribution only)
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (length(distribution) != 1 || !distribution %in% valid_distributions) {
    stop("Distribution must be exactly one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Define distribution characteristics
  distributions_with_shape <- c("gamma", "gengamma", "gompertz")
  distributions_without_shape <- c("exp", "weibull", "lnorm", "llogis")

  has_shape_param <- distribution %in% distributions_with_shape

  # Warn if include_shape is TRUE for distributions without shape parameters
  if (include_shape && !has_shape_param) {
    warning("include_shape = TRUE specified for distribution '", distribution,
            "' which has no shape parameter. Setting include_shape = FALSE.")
    include_shape <- FALSE
  }

  # Determine reference level for STRATIFY
  stratify_levels <- levels(dataset$STRATIFY)
  if (is.null(stratify_reference_level)) {
    reference_level <- stratify_levels[1]
    message("Using default reference level for STRATIFY: '", reference_level, "'")
  } else {
    if (!stratify_reference_level %in% stratify_levels) {
      stop("stratify_reference_level '", stratify_reference_level,
           "' not found in STRATIFY levels: ", paste(stratify_levels, collapse = ", "))
    }
    reference_level <- stratify_reference_level
    message("Using specified reference level for STRATIFY: '", reference_level, "'")
  }

  # Create indicator variable for stratification
  dataset$STRATIFY_INDICATOR <- ifelse(dataset$STRATIFY == reference_level, 1, 0)
  message("Stratification encoding: ", reference_level, " = 1, ",
          paste(setdiff(stratify_levels, reference_level), collapse = ", "), " = 0")

  # Check ARM levels
  arm_levels <- levels(dataset$ARM)
  if (length(arm_levels) < 2) {
    stop("ARM variable must have at least 2 levels for covariate analysis")
  }

  # Distribution name mapping
  dist_name_map <- c(
    "exp" = "Exponential",
    "weibull" = "Weibull",
    "lnorm" = "Lognormal",
    "llogis" = "Loglogistic",
    "gompertz" = "Gompertz",
    "gengamma" = "Generalized Gamma",
    "gamma" = "Gamma"
  )

  # Fit model based on distribution type
  tryCatch({

    if (distribution %in% c("exp", "weibull", "lnorm", "llogis")) {
      # Use survival::survreg for standard distributions
      survreg_dist <- switch(distribution,
                             "exp" = "exponential",
                             "weibull" = "weibull",
                             "lnorm" = "lognormal",
                             "llogis" = "loglogistic"
      )

      fit_model <- survival::survreg(
        survival::Surv(SURVTIME, EVENT) ~ ARM + STRATIFY_INDICATOR,
        data = dataset,
        dist = survreg_dist
      )

    } else {
      # Use flexsurv for advanced distributions
      fit_model <- flexsurv::flexsurvreg(
        survival::Surv(SURVTIME, EVENT) ~ ARM + STRATIFY_INDICATOR,
        data = dataset,
        dist = distribution
      )
    }

    # Extract full covariance matrix
    full_cov <- stats::vcov(fit_model)
    param_names <- rownames(full_cov)

    message("Available parameters for ", distribution, ": ", paste(param_names, collapse = ", "))

    # Find parameter indices
    intercept_idx <- grep("^\\(Intercept\\)$|^mu$", param_names)
    arm_idx <- grep("^ARM", param_names)  # Treatment effect
    stratify_idx <- grep("^STRATIFY_INDICATOR$", param_names)

    # Distribution-specific scale and shape parameters
    if (distribution == "exp") {
      scale_idx <- c()  # Exponential has no separate scale in survreg
      shape_idx <- c()
    } else if (distribution == "weibull") {
      scale_idx <- grep("^Log\\(scale\\)$", param_names)
      shape_idx <- c()  # Shape handled differently in survreg
    } else if (distribution %in% c("lnorm", "llogis")) {
      scale_idx <- grep("^Log\\(scale\\)$", param_names)
      shape_idx <- c()
    } else if (distribution == "gamma") {
      scale_idx <- grep("^rate$", param_names)
      shape_idx <- grep("^shape$", param_names)
      # Gamma has no separate intercept in flexsurv - rate serves as baseline
      if (length(intercept_idx) == 0) {
        message("Note: Gamma distribution uses 'rate' as baseline parameter (no separate intercept)")
      }
    } else if (distribution == "gengamma") {
      scale_idx <- grep("^sigma$", param_names)
      shape_idx <- grep("^Q$", param_names)
    } else if (distribution == "gompertz") {
      scale_idx <- grep("^rate$", param_names)
      shape_idx <- grep("^shape$", param_names)
      # Gompertz has no separate intercept in flexsurv - rate serves as baseline
      if (length(intercept_idx) == 0) {
        message("Note: Gompertz distribution uses 'rate' as baseline parameter (no separate intercept)")
      }
    }

    # Build parameter selection
    keep_indices <- c()
    param_labels <- c()

    # Handle distributions without separate intercept (gamma, gompertz)
    has_separate_intercept <- length(intercept_idx) > 0

    if (has_separate_intercept) {
      # Normal case: has intercept
      keep_indices <- c(keep_indices, intercept_idx[1])
      param_labels <- c(param_labels, "Intercept")

      # Add stratification effect
      if (length(stratify_idx) > 0) {
        keep_indices <- c(keep_indices, stratify_idx[1])
        param_labels <- c(param_labels, reference_level)
      }

      # Add scale parameter if available
      if (length(scale_idx) > 0) {
        keep_indices <- c(keep_indices, scale_idx[1])
        param_labels <- c(param_labels, "Scale")
      }

    } else {
      # Special case: gamma/gompertz without separate intercept
      # Rate parameter serves multiple roles
      if (length(scale_idx) > 0) {
        # First, add rate as baseline (intercept-like)
        keep_indices <- c(keep_indices, scale_idx[1])
        param_labels <- c(param_labels, "Intercept")  # Label as Intercept for consistency

        # Add stratification effect
        if (length(stratify_idx) > 0) {
          keep_indices <- c(keep_indices, stratify_idx[1])
          param_labels <- c(param_labels, reference_level)
        }

        # Add the same rate parameter again as Scale for consistency with other distributions
        keep_indices <- c(keep_indices, scale_idx[1])
        param_labels <- c(param_labels, "Scale")

        message("Note: For ", distribution, ", 'rate' parameter appears as both Intercept and Scale")
      }
    }

    # Include shape parameter if requested and available
    if (include_shape && has_shape_param && length(shape_idx) > 0) {
      keep_indices <- c(keep_indices, shape_idx[1])
      param_labels <- c(param_labels, "Shape")
      message("Including shape parameter in covariance matrix")
    } else if (has_shape_param) {
      message("Shape parameter available but excluded (include_shape = FALSE)")
    }

    # Ensure we have enough parameters
    if (length(keep_indices) < 2) {
      stop("Insufficient parameters found for covariance matrix (need at least 2)")
    }

    # For distributions where same parameter appears twice, handle covariance matrix carefully
    if (!has_separate_intercept && distribution %in% c("gamma", "gompertz")) {
      # Extract covariance submatrix with unique indices first
      unique_indices <- unique(keep_indices)
      sub_cov_unique <- full_cov[unique_indices, unique_indices, drop = FALSE]

      # Create expanded matrix to match param_labels structure
      n_params <- length(param_labels)
      expanded_cov <- matrix(0, nrow = n_params, ncol = n_params)

      # Map to expanded matrix
      param_to_unique <- match(keep_indices, unique_indices)
      for (i in 1:n_params) {
        for (j in 1:n_params) {
          expanded_cov[i, j] <- sub_cov_unique[param_to_unique[i], param_to_unique[j]]
        }
      }

      sub_cov <- expanded_cov

    } else {
      # Normal case: extract covariance submatrix
      sub_cov <- full_cov[keep_indices, keep_indices, drop = FALSE]
    }

    # Set proper row and column names
    rownames(sub_cov) <- param_labels
    colnames(sub_cov) <- param_labels

    # Return raw matrix if requested
    if (!format_output) {
      return(sub_cov)
    }

    # Convert to data frame format for output
    cov_df <- data.frame(sub_cov, check.names = FALSE)
    cov_df$Variable <- rownames(sub_cov)

    # Add metadata
    cov_df$Distribution <- dist_name_map[distribution]
    cov_df$Parameter <- parameter_name

    # Reorder columns: metadata first, then covariance values
    meta_cols <- c("Distribution", "Parameter", "Variable")
    cov_cols <- setdiff(names(cov_df), meta_cols)
    cov_df <- cov_df[, c(meta_cols, cov_cols)]

    # Round numeric values
    cov_df[cov_cols] <- lapply(cov_df[cov_cols], function(x) round(x, 10))

    # Create kable if requested
    if (create_kable) {
      formatted_table <- cov_df %>%
        kableExtra::kable(
          caption = paste("Covariance Matrix -", dist_name_map[distribution], "-", parameter_name),
          align = c('l', 'l', 'l', rep('c', length(cov_cols))),
          digits = 6,
          row.names = FALSE  # Prevent automatic row names column
        ) %>%
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed", "responsive"),
          full_width = FALSE,
          position = "center"
        ) %>%
        kableExtra::column_spec(1, bold = TRUE, color = "darkblue") %>%
        kableExtra::column_spec(2, bold = TRUE, color = "darkgreen") %>%
        kableExtra::collapse_rows(columns = c(1, 2), valign = "middle")  # Merge Distribution and Parameter columns

      return(formatted_table)
    }

    return(cov_df)

  }, error = function(e) {
    stop("Failed to process distribution '", distribution, "': ", e$message)
  })
}
