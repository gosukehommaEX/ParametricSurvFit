#' Create Kaplan-Meier Curves with Overlaid Parametric Distribution Fits
#'
#' This function creates Kaplan-Meier survival curves by treatment arm and overlays
#' fitted parametric distribution curves. Separate panels are created for each
#' stratification factor level, with parameter estimates displayed on each panel.
#' Uses survminer package for professional survival curve plotting.
#'
#' @param dataset A data frame created by DataParametricSurv() function containing
#'   survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT
#' @param distribution Character string specifying which distribution to fit and overlay.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param figure_caption Character string for figure caption. Default is
#'   "Kaplan-Meier Curves with Parametric Distribution Overlay"
#' @param stratify_name Character string specifying the original stratification
#'   variable name for panel titles. If NULL (default), uses "STRATIFY"
#' @param control_arm Character string specifying which ARM should be treated as control.
#'   Control arm will be plotted in red, treatment arm in blue. If NULL (default),
#'   uses first ARM level as control
#' @param time_scale Character string specifying the time scale for x-axis label.
#'   Default is "Months". Common options: "Months", "Years", "Weeks", "Days"
#' @param time_max Numeric value for maximum time on x-axis. If NULL (default),
#'   uses maximum observed survival time
#' @param return_plots Logical indicating whether to return plot objects.
#'   If TRUE (default), returns list of ggplot objects. If FALSE, displays plots
#'
#' @return A list of ggsurvplot objects (one for each stratification level) or displays plots
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Creates Kaplan-Meier curves for each ARM using survminer::ggsurvplot()
#'   \item Fits specified parametric distribution using flexsurv::flexsurvreg()
#'   \item Overlays parametric curves as dashed lines matching ARM colors
#'   \item Creates separate panels for each stratification factor level
#'   \item Displays parameter estimates as text annotations on each panel
#'   \item Includes risk table with proper alignment to x-axis
#' }
#'
#' Colors used:
#' \itemize{
#'   \item Control arm = red, Treatment arm = blue
#'   \item Parametric distribution lines use matching colors but with dashed line type
#' }
#'
#' @importFrom dplyr group_by filter mutate arrange summarise n
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 geom_line annotate
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Create survival dataset first
#' surv_data <- DataParametricSurv(
#'   adsl_path = "path/to/adsl.sas7bdat",
#'   adtte_path = "path/to/adtte.sas7bdat",
#'   population = "ITTFL",
#'   variable = "OS",
#'   stratify_by = "SEX"
#' )
#'
#' # Create KM plots with Weibull overlay, Placebo as control, custom time scale
#' plots_weibull <- ParametricSurvKM(
#'   dataset = surv_data,
#'   distribution = "weibull",
#'   figure_caption = "Overall Survival Analysis",
#'   stratify_name = "SEX",
#'   control_arm = "Placebo",
#'   time_scale = "Years"
#' )
#'
#' # Display plots
#' print(plots_weibull)
#'
#' # Create plots with exponential overlay and custom time scale
#' plots_exp <- ParametricSurvKM(
#'   dataset = surv_data,
#'   distribution = "exp",
#'   figure_caption = "Exponential Model Fit",
#'   stratify_name = "SEX",
#'   time_scale = "Weeks",
#'   time_max = 100
#' )
#'
#' # Save individual plots
#' ggsave("OS_SEX_M_weibull.png", plots_weibull[[1]]$plot, width = 12, height = 8)
#' ggsave("OS_SEX_F_weibull.png", plots_weibull[[2]]$plot, width = 12, height = 8)
#' }
ParametricSurvKM <- function(dataset,
                             distribution,
                             figure_caption = "Kaplan-Meier Curves with Parametric Distribution Overlay",
                             stratify_name = NULL,
                             control_arm = NULL,
                             time_scale = "Months",
                             time_max = NULL,
                             return_plots = TRUE) {

  # Define variables for NSE to avoid R CMD check notes
  ARM <- STRATIFY <- SURVTIME <- EVENT <- CNSR <- NULL
  time <- surv <- NULL

  # Load required packages
  required_packages <- c("ggplot2", "survival", "flexsurv", "survminer", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.")
    }
  }

  # Validate input dataset
  required_cols <- c("SUBJID", "ARM", "STRATIFY", "SURVTIME", "CNSR", "EVENT")
  missing_cols <- setdiff(required_cols, names(dataset))
  if (length(missing_cols) > 0) {
    stop("Dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate distribution
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (length(distribution) != 1 || !distribution %in% valid_distributions) {
    stop("Distribution must be exactly one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Set stratify_name if not provided
  if (is.null(stratify_name)) {
    stratify_name <- "STRATIFY"
  }

  # Get unique stratification levels
  stratify_levels <- unique(dataset$STRATIFY)
  arm_levels <- unique(dataset$ARM)

  # Check if we have stratification (more than one STRATIFY level)
  has_stratification <- length(stratify_levels) > 1

  # Validate ARM levels
  if (length(arm_levels) != 2) {
    stop("This function requires exactly 2 ARM levels, but found: ", length(arm_levels))
  }

  # Set control arm if not specified
  if (is.null(control_arm)) {
    control_arm <- arm_levels[1]
    message("Using default control arm: ", control_arm)
  } else if (!control_arm %in% arm_levels) {
    stop("Specified control_arm '", control_arm, "' not found in ARM levels: ",
         paste(arm_levels, collapse = ", "))
  }

  # Set colors: control arm = red, treatment arm = blue
  treatment_arm <- setdiff(arm_levels, control_arm)
  arm_colors <- c("red", "blue")
  names(arm_colors) <- c(control_arm, treatment_arm)

  # Set time_max if not provided
  if (is.null(time_max)) {
    time_max <- max(dataset$SURVTIME, na.rm = TRUE) * 1.1
  }

  # Distribution name mapping for display
  dist_name_map <- c(
    "exp" = "Exponential",
    "weibull" = "Weibull",
    "lnorm" = "Lognormal",
    "llogis" = "Loglogistic",
    "gompertz" = "Gompertz",
    "gengamma" = "Generalized Gamma",
    "gamma" = "Gamma"
  )

  # Initialize list to store plots
  plot_list <- list()

  # Create plot for each stratification level
  for (strat_level in stratify_levels) {

    # Filter data for current stratification level
    current_data <- dataset %>%
      dplyr::filter(STRATIFY == strat_level)

    # Fit Kaplan-Meier
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ ARM, data = current_data)

    # Get the actual strata names from km_fit and determine correct color order
    strata_names <- names(km_fit$strata)

    # Create color vector in the correct order based on strata names
    color_vector <- c()
    for (strata_name in strata_names) {
      # Extract ARM name from strata (format: "ARM=ArmName")
      arm_name <- gsub("ARM=", "", strata_name)
      if (arm_name == control_arm) {
        color_vector <- c(color_vector, "red")
      } else {
        color_vector <- c(color_vector, "blue")
      }
    }

    # Create title and subtitle based on stratification
    if (has_stratification) {
      plot_title <- paste0(figure_caption, " - (", stratify_name, " = ", strat_level, ") -")
    } else {
      plot_title <- paste0(figure_caption, " - Overall Population -")
    }

    # Create base ggsurvplot with minimal settings to avoid warnings
    base_plot <- suppressWarnings(survminer::ggsurvplot(
      km_fit,
      data = current_data,
      palette = color_vector,  # Use correctly ordered color vector
      conf.int = FALSE,
      pval = FALSE,
      risk.table = TRUE,
      risk.table.col = "black",  # Use black instead of strata to avoid color conflicts
      risk.table.height = 0.3,
      ncensor.plot = FALSE,
      ggtheme = ggplot2::theme_minimal(),
      tables.theme = survminer::theme_cleantable(),
      title = plot_title,
      subtitle = paste("Distribution:", dist_name_map[distribution]),
      xlab = paste0("Time (", time_scale, ")"),  # Dynamic x-axis label based on time_scale
      ylab = "Probability of Survival",
      xlim = c(0, time_max),
      ylim = c(0, 1),
      break.x.by = max(5, round(time_max/12)),
      break.y.by = 0.25,
      legend = "none",  # Remove the default ARM legend completely
      font.main = c(14, "bold", "black"),
      font.subtitle = c(12, "plain", "black"),
      font.x = c(12, "plain", "black"),
      font.y = c(12, "plain", "black"),
      risk.table.title = "Number at risk",
      risk.table.fontsize = 3.5,
      tables.y.text = FALSE
    ))

    # Fit parametric models and prepare data for overlay
    parametric_fits <- list()
    param_estimates <- list()
    arm_summary <- list()

    for (arm in arm_levels) {
      arm_data <- current_data %>%
        dplyr::filter(ARM == arm)

      # Calculate N and Events
      n_subjects <- nrow(arm_data)
      n_events <- sum(arm_data$EVENT)
      arm_summary[[as.character(arm)]] <- list(N = n_subjects, Events = n_events)

      # Parametric fit - use the appropriate distribution type
      tryCatch({
        if (distribution %in% c("exp", "weibull", "lnorm", "llogis")) {
          # Use flexsurv for consistency
          parametric_fits[[as.character(arm)]] <- flexsurv::flexsurvreg(
            survival::Surv(SURVTIME, EVENT) ~ 1,
            data = arm_data,
            dist = distribution
          )
        } else {
          # Use flexsurv for advanced distributions (gompertz, gengamma, gamma)
          parametric_fits[[as.character(arm)]] <- flexsurv::flexsurvreg(
            survival::Surv(SURVTIME, EVENT) ~ 1,
            data = arm_data,
            dist = distribution
          )
        }
      }, error = function(e) {
        warning("Failed to fit ", distribution, " distribution for ARM ", arm, ": ", e$message)
        parametric_fits[[as.character(arm)]] <- NULL
      })

      # Extract parameters for display if fit was successful
      if (!is.null(parametric_fits[[as.character(arm)]])) {
        tryCatch({
          param_info <- ExtractParams(parametric_fits[[as.character(arm)]], distribution)
          param_text <- paste0(
            param_info$Parameters, ": ",
            round(param_info$Estimates, 3),
            collapse = ", "
          )
          param_estimates[[as.character(arm)]] <- param_text
        }, error = function(e) {
          warning("Failed to extract parameters for ARM ", arm, ": ", e$message)
          param_estimates[[as.character(arm)]] <- "Parameter extraction failed"
        })
      } else {
        param_estimates[[as.character(arm)]] <- "Model fit failed"
      }
    }

    # Create time sequence for parametric curves
    time_seq <- seq(0, time_max, length.out = 200)

    # Add parametric curves to the plot (using correct ARM colors)
    for (arm in arm_levels) {
      param_fit <- parametric_fits[[as.character(arm)]]

      if (!is.null(param_fit)) {
        tryCatch({
          # Get survival predictions
          param_surv <- summary(param_fit, t = time_seq, type = "survival")[[1]]$est

          # Determine color based on control_arm
          line_color <- if (arm == control_arm) "red" else "blue"

          # Add parametric curve as dashed line
          base_plot$plot <- base_plot$plot +
            ggplot2::geom_line(
              data = data.frame(time = time_seq, surv = param_surv),
              ggplot2::aes(x = time, y = surv),
              color = line_color,
              linetype = "dashed",
              linewidth = 1,
              inherit.aes = FALSE
            )
        }, error = function(e) {
          warning("Failed to add parametric curve for ARM ", arm, ": ", e$message)
        })
      }
    }

    # Add legend for line types in bottom left (moved up slightly and arranged horizontally)
    # Using Unicode escapes for special characters to avoid non-ASCII warnings
    legend_text <- paste0(
      "\u2014\u2014 Kaplan-Meier     - - - Parametric     + Censored"
    )

    base_plot$plot <- base_plot$plot +
      ggplot2::annotate(
        "text",
        x = time_max * 0.02,
        y = 0.28,  # Moved up slightly from 0.35 to 0.28
        label = legend_text,  # Horizontal arrangement with spacing
        hjust = 0,
        vjust = 1,
        size = 3,
        color = "black"
      )

    # Add comprehensive legend with colored text matching plot lines
    # Treatment arm (blue) text - positioned above control arm
    treatment_legend <- paste0(treatment_arm, ": Events/N ", arm_summary[[treatment_arm]]$Events, "/",
                               arm_summary[[treatment_arm]]$N, ", ", param_estimates[[treatment_arm]])

    base_plot$plot <- base_plot$plot +
      ggplot2::annotate(
        "text",
        x = time_max * 0.02,
        y = 0.20,  # Treatment arm at higher position
        label = treatment_legend,
        hjust = 0,
        vjust = 1,
        size = 3,
        color = "blue",  # Blue color to match treatment arm
        fontface = "bold"
      )

    # Control arm (red) text - positioned below treatment arm
    control_legend <- paste0(control_arm, ": Events/N ", arm_summary[[control_arm]]$Events, "/",
                             arm_summary[[control_arm]]$N, ", ", param_estimates[[control_arm]])

    base_plot$plot <- base_plot$plot +
      ggplot2::annotate(
        "text",
        x = time_max * 0.02,
        y = 0.12,  # Control arm at lower position
        label = control_legend,
        hjust = 0,
        vjust = 1,
        size = 3,
        color = "red",  # Red color to match control arm
        fontface = "bold"
      )

    # Store plot with appropriate naming
    if (has_stratification) {
      plot_list[[paste0(stratify_name, "_", strat_level)]] <- base_plot
    } else {
      plot_list[["Overall"]] <- base_plot
    }

    # Display plot if return_plots is FALSE
    if (!return_plots) {
      print(base_plot)
    }
  }

  # Return plots if requested
  if (return_plots) {
    return(plot_list)
  }
}
