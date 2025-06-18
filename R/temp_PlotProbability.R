#' Create Probability Plot with 95% Confidence Limits for ARM*STRATIFY Groups
#'
#' This function creates probability plots with confidence intervals for
#' parametric survival models fitted using flexsurv. It creates separate plots
#' for each ARM*STRATIFY combination and integrates with the existing
#' DataParametricSurv and ExtractParams functions.
#'
#' @param dataset A data frame created by DataParametricSurv() function containing
#'   survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT
#' @param distribution Character string specifying the distribution type.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param time_scale Character string specifying the time scale for x-axis label.
#'   Default is "Months". Common options: "Months", "Years", "Weeks", "Days"
#' @param confidence_level Numeric value for confidence level (default: 0.95)
#' @param point_size Numeric value for point size (default: 2)
#' @param line_size Numeric value for line size (default: 1)
#' @param alpha_ribbon Numeric value for ribbon transparency (default: 0.3)
#' @param return_plots Logical indicating whether to return plot objects.
#'   If TRUE (default), returns list of ggplot objects. If FALSE, displays plots
#'
#' @return A list of ggplot2 objects containing probability plots for each ARM*STRATIFY combination
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_segment
#'   labs theme_minimal theme element_text scale_x_log10 scale_y_continuous
#'   annotation_custom unit
#' @importFrom survival survfit Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom dplyr filter mutate arrange group_by
#' @importFrom stats qnorm pexp pweibull plnorm pgamma
#' @importFrom grid textGrob gpar
#' @importFrom scales pretty_breaks
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
#' #' # Create probability plots with Weibull distribution
#' prob_plots <- PlotProbability(
#'   dataset = surv_data,
#'   distribution = "weibull",
#'   time_scale = "Months"
#' )
#'
#' # Display all plots
#' print(prob_plots)
#'
#' # Access individual plots
#' treatment_male_plot <- prob_plots[["Treatment_M"]]
#' placebo_female_plot <- prob_plots[["Placebo_F"]]
#' }
PlotProbability <- function(dataset,
                            distribution,
                            time_scale = "Months",
                            confidence_level = 0.95,
                            point_size = 2,
                            line_size = 1,
                            alpha_ribbon = 0.3,
                            return_plots = TRUE) {

  # Validate distribution
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!distribution %in% valid_distributions) {
    stop("Distribution must be one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Validate dataset columns
  required_cols <- c("SUBJID", "ARM", "STRATIFY", "SURVTIME", "CNSR", "EVENT")
  missing_cols <- setdiff(required_cols, names(dataset))
  if (length(missing_cols) > 0) {
    stop("Dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Define variables for NSE to avoid R CMD check notes
  ARM <- STRATIFY <- SURVTIME <- EVENT <- CNSR <- NULL
  time <- surv <- n_event <- cdf <- NULL

  # Get unique ARM*STRATIFY combinations
  combinations <- dataset %>%
    dplyr::group_by(ARM, STRATIFY) %>%
    dplyr::summarise(.groups = "drop") %>%
    dplyr::mutate(group_name = paste(ARM, STRATIFY, sep = "_"))

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

  # Initialize list to store plots
  plot_list <- list()

  # Create plot for each ARM*STRATIFY combination
  for (i in 1:nrow(combinations)) {
    arm_val <- combinations$ARM[i]
    stratify_val <- combinations$STRATIFY[i]
    group_name <- combinations$group_name[i]

    # Filter data for current ARM*STRATIFY combination
    current_data <- dataset %>%
      dplyr::filter(ARM == arm_val & STRATIFY == stratify_val) %>%
      dplyr::filter(!is.na(SURVTIME) & !is.na(EVENT) & SURVTIME > 0)

    if (nrow(current_data) == 0) {
      warning("No valid data for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    # Fit the specified model
    fit_model <- tryCatch({
      flexsurv::flexsurvreg(
        survival::Surv(SURVTIME, EVENT) ~ 1,
        data = current_data,
        dist = distribution
      )
    }, error = function(e) {
      warning("Failed to fit model for ARM: ", arm_val, ", STRATIFY: ", stratify_val, " - ", e$message)
      return(NULL)
    })

    if (is.null(fit_model)) {
      next
    }

    # Calculate Kaplan-Meier estimates
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = current_data)

    # Extract empirical data (only for actual events)
    km_summary <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      n_event = km_fit$n.event
    )

    # Filter to only include time points where actual events occurred
    empirical_data <- km_summary %>%
      dplyr::filter(n_event > 0) %>%
      dplyr::mutate(cdf = 1 - surv) %>%
      dplyr::filter(cdf > 0.001 & cdf < 0.999) %>%
      dplyr::arrange(time)

    if (nrow(empirical_data) == 0) {
      warning("No valid empirical data points for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    # Create time sequence for model predictions
    time_min <- min(current_data$SURVTIME) * 0.5
    time_max <- max(current_data$SURVTIME) * 1.2
    time_seq <- seq(time_min, time_max, length.out = 200)

    # Get model predictions with confidence intervals
    model_pred <- tryCatch({
      summary(fit_model, t = time_seq, type = "cumhaz", ci = TRUE)
    }, error = function(e) {
      # Fallback to basic prediction without CI
      basic_pred <- summary(fit_model, t = time_seq, type = "survival")[[1]]
      data.frame(
        time = time_seq,
        est = 1 - basic_pred$est,
        lcl = NA,
        ucl = NA
      )
    })

    # Handle model prediction format
    if ("[[" %in% class(model_pred) || is.list(model_pred)) {
      if (length(model_pred) > 0 && is.data.frame(model_pred[[1]])) {
        pred_df <- model_pred[[1]]
        theoretical_data <- data.frame(
          time = time_seq,
          fitted_cdf = 1 - exp(-pred_df$est),  # Convert cumhaz to cdf
          lower_cdf = if (!is.null(pred_df$lcl)) 1 - exp(-pred_df$ucl) else NA,
          upper_cdf = if (!is.null(pred_df$ucl)) 1 - exp(-pred_df$lcl) else NA
        )
      } else {
        # Fallback calculation
        surv_pred <- summary(fit_model, t = time_seq, type = "survival")[[1]]
        theoretical_data <- data.frame(
          time = time_seq,
          fitted_cdf = 1 - surv_pred$est,
          lower_cdf = NA,
          upper_cdf = NA
        )
      }
    } else {
      theoretical_data <- data.frame(
        time = time_seq,
        fitted_cdf = model_pred$est,
        lower_cdf = if (!is.null(model_pred$lcl)) model_pred$lcl else NA,
        upper_cdf = if (!is.null(model_pred$ucl)) model_pred$ucl else NA
      )
    }

    # Filter theoretical data
    theoretical_data <- theoretical_data %>%
      dplyr::filter(fitted_cdf > 0.001 & fitted_cdf < 0.999)

    # Extract model parameters for legend
    model_params <- fit_model$res[, "est"]
    names(model_params) <- rownames(fit_model$res)

    # Create model statistics for legend
    n_total <- nrow(current_data)
    n_events <- sum(current_data$EVENT)
    n_censored <- n_total - n_events

    # Extract parameters using ExtractParams function for consistency
    extracted_params <- ExtractParams(fit_model, distribution)

    # Create parameter text using ExtractParams results with right-aligned values (v30 style)
    param_lines <- c()
    for (i in seq_along(extracted_params$Parameters)) {
      param_name <- extracted_params$Parameters[i]
      param_value <- sprintf("%.3f", extracted_params$Estimates[i])

      # Use fixed spacing like v30 - match the exact format of "Distribution    "
      if (param_name == "Intercept") {
        param_line <- paste0("Intercept       ", sprintf("%12s", param_value))
      } else if (param_name == "Scale") {
        param_line <- paste0("Scale           ", sprintf("%12s", param_value))
      } else if (param_name == "Shape") {
        param_line <- paste0("Shape           ", sprintf("%12s", param_value))
      } else if (param_name == "Rate") {
        param_line <- paste0("Rate            ", sprintf("%12s", param_value))
      } else if (param_name == "Meanlog") {
        param_line <- paste0("Meanlog         ", sprintf("%12s", param_value))
      } else if (param_name == "Sdlog") {
        param_line <- paste0("Sdlog           ", sprintf("%12s", param_value))
      } else {
        # Fallback for any other parameter names
        spacing_needed <- max(0, 16 - nchar(param_name))
        param_line <- paste0(param_name, paste(rep(" ", spacing_needed), collapse = ""), sprintf("%12s", param_value))
      }

      param_lines <- c(param_lines, param_line)
    }

    # Combine all parameter lines
    param_text <- paste(param_lines, collapse = "\n")

    # Create complete legend text with right-aligned values
    legend_text <- paste0(
      "Distribution    ", sprintf("%12s", dist_name_map[distribution]), "\n",
      "Observations    ", sprintf("%12d", n_total), "\n",
      "Event           ", sprintf("%12d", n_events), "\n",
      "Censored        ", sprintf("%12d", n_censored), "\n",
      param_text
    )

    # Create plot title
    plot_title <- paste0(
      dist_name_map[distribution], " Probability Plot With ", confidence_level * 100, "% Confidence Limits\n",
      "ARM: ", arm_val, ", STRATIFY: ", stratify_val
    )

    # Prepare data for legend display
    # Create event data for points
    event_data <- empirical_data %>%
      dplyr::mutate(type = "Events")

    # Create censored data for rug plot
    censored_data <- current_data %>%
      dplyr::filter(EVENT == 0) %>%
      dplyr::mutate(type = "Right Censored",
                    cdf = 0.01)  # Small y value for rug plot

    # Create the plot with proper legend
    p <- ggplot2::ggplot() +
      # Add confidence ribbon (if available)
      {if (!all(is.na(theoretical_data$lower_cdf))) {
        ggplot2::geom_ribbon(
          data = theoretical_data,
          ggplot2::aes(x = time, ymin = lower_cdf, ymax = upper_cdf),
          alpha = alpha_ribbon,
          fill = "lightblue"
        )
      }} +
      # Add fitted line
      ggplot2::geom_line(
        data = theoretical_data,
        ggplot2::aes(x = time, y = fitted_cdf),
        color = "blue",
        size = line_size
      ) +
      # Add empirical points (events) with legend
      ggplot2::geom_point(
        data = event_data,
        ggplot2::aes(x = time, y = cdf, shape = type),
        size = point_size,
        color = "black",
        fill = "white"
      ) +
      # Add rug plot for censored observations with legend
      {if (nrow(censored_data) > 0) {
        ggplot2::geom_point(
          data = censored_data,
          ggplot2::aes(x = SURVTIME, y = cdf, shape = type),
          size = 3,
          color = "black"
        )
      }} +
      # Set shapes manually
      ggplot2::scale_shape_manual(
        name = "",
        values = c("Events" = 21, "Right Censored" = 124),  # 21 = circle, 124 = |
        guide = ggplot2::guide_legend(override.aes = list(size = c(3, 4)))
      ) +
      # Use log scale for x-axis
      ggplot2::scale_x_log10() +
      # Set y-axis as probability scale starting from 0%
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = paste0(seq(0, 100, 20), "%"),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      # Customize axes and labels
      ggplot2::labs(
        title = plot_title,
        x = paste0("Time to Event (", time_scale, ")"),
        y = "Cumulative Probability"
      ) +
      # Add theme with legend at bottom
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 9),
        panel.grid.minor = ggplot2::element_line(color = "grey90", size = 0.3),
        panel.grid.major = ggplot2::element_line(color = "grey80", size = 0.5),
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.box = "horizontal",
        legend.margin = ggplot2::margin(t = 10),
        legend.box.background = ggplot2::element_rect(fill = "white", color = "black", size = 0.5),
        legend.box.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10)
      )

    # Add statistics box (no need for separate legend creation)
    p <- p +
      ggplot2::annotation_custom(
        grob = grid::textGrob(
          legend_text,
          x = 0.02, y = 0.98,
          just = c("left", "top"),
          gp = grid::gpar(fontsize = 8, fontfamily = "mono")
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      )

    # Store plot
    plot_list[[group_name]] <- p

    # Display plot if return_plots is FALSE
    if (!return_plots) {
      print(p)
    }
  }

  # Return plots if requested
  if (return_plots) {
    return(plot_list)
  }
}
