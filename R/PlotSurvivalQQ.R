#' Create Survival Function Q-Q Plot with log(-log(S)) transformation
#'
#' This function creates Q-Q plots using S(t) directly with log(-log(S)) transformation
#' for parametric survival models fitted using flexsurv. It creates separate plots
#' for each ARM*STRATIFY combination with native flexsurv confidence intervals.
#'
#' @param dataset A data frame created by DataParametricSurv() function containing
#'   survival data with columns: SUBJID, ARM, STRATIFY, SURVTIME, CNSR, EVENT
#' @param distribution Character string specifying the distribution type.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param time_scale Character string specifying the time scale for x-axis label.
#'   Default is "Months". Common options: "Months", "Years", "Weeks", "Days"
#' @param stratify_name Character string specifying the original stratification
#'   variable name for plot titles. If NULL (default), uses "STRATIFY"
#' @param confidence_level Numeric value for confidence level (default: 0.95)
#' @param point_size Numeric value for point size (default: 2)
#' @param line_size Numeric value for line size (default: 1)
#' @param alpha_ribbon Numeric value for ribbon transparency (default: 0.3)
#' @param return_plots Logical indicating whether to return plot objects.
#'   If TRUE (default), returns list of ggplot objects. If FALSE, displays plots
#'
#' @return A list of ggplot2 objects containing S(t)-based Q-Q plots for each ARM*STRATIFY combination
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon
#'   labs theme_minimal theme element_text scale_x_log10 scale_y_continuous
#'   annotation_custom unit element_rect margin
#' @importFrom survival survfit Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom dplyr filter mutate arrange group_by summarise
#' @importFrom stats qnorm
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
#' # Create S(t)-based Q-Q plots with Weibull distribution
#' qq_plots <- PlotSurvivalQQ(
#'   dataset = surv_data,
#'   distribution = "weibull",
#'   time_scale = "Months"
#' )
#'
#' # Display all plots
#' print(qq_plots)
#'
#' # Access individual plots
#' treatment_male_plot <- qq_plots[["Treatment_M"]]
#' placebo_female_plot <- qq_plots[["Placebo_F"]]
#' }
PlotSurvivalQQ <- function(dataset,
                           distribution,
                           time_scale = "Months",
                           stratify_name = NULL,
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
  time <- surv <- n_event <- NULL

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

  # Set stratify_name if not provided
  if (is.null(stratify_name)) {
    stratify_name <- "STRATIFY"
  }

  # Calculate S(t) confidence intervals using flexsurv
  calculate_survival_confidence <- function(fit_model, time_seq, confidence_level) {

    # Use flexsurv's built-in confidence interval calculation for S(t)
    # Remove B=0 to enable confidence intervals
    surv_summary <- summary(fit_model, t = time_seq, type = "survival", ci = TRUE)[[1]]

    # Extract survival estimates and confidence bounds
    surv_est <- surv_summary$est
    surv_lower <- surv_summary$lcl
    surv_upper <- surv_summary$ucl

    # Check if confidence intervals are available
    if (is.null(surv_lower) || is.null(surv_upper)) {
      warning("Confidence intervals not available from flexsurv")
      surv_lower <- surv_est
      surv_upper <- surv_est
    }

    return(data.frame(
      time = time_seq,
      fitted_surv = surv_est,
      lower_surv = surv_lower,
      upper_surv = surv_upper
    ))
  }

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
      dplyr::mutate(
        log_neg_log_surv = log(-log(surv))  # y-axis transformation: log(-log(S))
      ) %>%
      dplyr::filter(surv > 0.001 & surv < 0.999, is.finite(log_neg_log_surv)) %>%
      dplyr::arrange(time)

    if (nrow(empirical_data) == 0) {
      warning("No valid empirical data points for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    # Create time sequence for model predictions - use data-driven range
    time_min <- min(current_data$SURVTIME[current_data$SURVTIME > 0]) * 0.5
    time_max <- max(current_data$SURVTIME) * 1.5
    time_seq <- exp(seq(log(time_min), log(time_max), length.out = 300))

    # Get model predictions with native flexsurv confidence intervals
    theoretical_data <- calculate_survival_confidence(
      fit_model = fit_model,
      time_seq = time_seq,
      confidence_level = confidence_level
    )

    # Apply transformations - ensure survival values are within valid range
    theoretical_data <- theoretical_data %>%
      dplyr::mutate(
        fitted_surv_safe = pmax(0.001, pmin(0.999, fitted_surv)),
        lower_surv_safe = pmax(0.001, pmin(0.999, lower_surv)),
        upper_surv_safe = pmax(0.001, pmin(0.999, upper_surv)),
        # Apply log(-log(S)) transformation
        log_neg_log_fitted = log(-log(fitted_surv_safe)),
        log_neg_log_lower = log(-log(upper_surv_safe)),  # Note: CI bounds flip
        log_neg_log_upper = log(-log(lower_surv_safe))
      ) %>%
      dplyr::filter(is.finite(log_neg_log_fitted))

    # Calculate y-axis range based on both empirical data AND theoretical confidence intervals
    empirical_y_values <- empirical_data$log_neg_log_surv[is.finite(empirical_data$log_neg_log_surv)]
    if (length(empirical_y_values) == 0) {
      warning("No finite y values for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    # Get theoretical Y range from confidence intervals
    theoretical_y_values <- c(
      theoretical_data$log_neg_log_fitted[is.finite(theoretical_data$log_neg_log_fitted)],
      theoretical_data$log_neg_log_lower[is.finite(theoretical_data$log_neg_log_lower)],
      theoretical_data$log_neg_log_upper[is.finite(theoretical_data$log_neg_log_upper)]
    )

    # Combine empirical and theoretical ranges
    all_y_values <- c(empirical_y_values, theoretical_y_values)
    y_min <- min(all_y_values, na.rm = TRUE) - 0.2  # Add padding
    y_max <- max(all_y_values, na.rm = TRUE) + 0.2

    # Calculate breaks and labels for y-axis with expanded range
    desired_percentiles <- c(99.9, 99, 98, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 2, 1, 0.1)
    desired_survival <- desired_percentiles / 100
    desired_survival[desired_survival <= 0] <- 0.001
    desired_survival[desired_survival >= 1] <- 0.999
    desired_breaks <- log(-log(desired_survival))

    # Filter to include breaks within the expanded range
    valid_indices <- which(desired_breaks >= y_min & desired_breaks <= y_max)

    # Ensure we have at least some key percentiles
    if (length(valid_indices) < 3) {
      key_indices <- which(desired_percentiles %in% c(99, 90, 50, 10, 1))
      valid_indices <- sort(unique(c(valid_indices, key_indices)))
      valid_indices <- valid_indices[desired_breaks[valid_indices] >= y_min & desired_breaks[valid_indices] <= y_max]
    }

    final_breaks <- desired_breaks[valid_indices]
    final_labels <- as.character(desired_percentiles[valid_indices])

    # Extract model parameters for legend
    n_total <- nrow(current_data)
    n_events <- sum(current_data$EVENT)
    n_censored <- n_total - n_events

    # Extract parameters using ExtractParams function
    extracted_params <- ExtractParams(fit_model, distribution)

    # Create parameter text
    param_lines <- c()
    for (j in seq_along(extracted_params$Parameters)) {
      param_name <- extracted_params$Parameters[j]
      param_value <- sprintf("%.3f", extracted_params$Estimates[j])

      if (param_name == "Intercept") {
        param_line <- paste0("Intercept       ", sprintf("%12s", param_value))
      } else if (param_name == "Scale") {
        param_line <- paste0("Scale           ", sprintf("%12s", param_value))
      } else if (param_name == "Shape") {
        param_line <- paste0("Shape           ", sprintf("%12s", param_value))
      } else if (param_name == "Rate") {
        param_line <- paste0("Rate            ", sprintf("%12s", param_value))
      } else {
        spacing_needed <- max(0, 16 - nchar(param_name))
        param_line <- paste0(param_name, paste(rep(" ", spacing_needed), collapse = ""), sprintf("%12s", param_value))
      }

      param_lines <- c(param_lines, param_line)
    }

    param_text <- paste(param_lines, collapse = "\n")

    # Create complete legend text
    legend_text <- paste0(
      "Distribution    ", sprintf("%12s", dist_name_map[distribution]), "\n",
      "Observations    ", sprintf("%12d", n_total), "\n",
      "Event           ", sprintf("%12d", n_events), "\n",
      "Censored        ", sprintf("%12d", n_censored), "\n",
      param_text
    )

    # Create plot title
    plot_title <- paste0(
      dist_name_map[distribution], " Survival Q-Q Plot With ", confidence_level * 100, "% Confidence Limits\n",
      "ARM: ", arm_val, ", ", stratify_name, ": ", stratify_val
    )

    # Create event data for points
    event_data <- empirical_data %>%
      dplyr::mutate(type = "Events")

    # Create censored data positioned within the plot area (not at minimum)
    censored_data <- current_data %>%
      dplyr::filter(EVENT == 0) %>%
      dplyr::filter(SURVTIME >= time_min & SURVTIME <= time_max) %>%
      dplyr::mutate(
        type = "Right Censored",
        log_neg_log_surv = y_min + 0.1  # Slightly above minimum for visibility
      )

    # Determine x-axis limits and breaks
    x_min <- min(time_min, min(empirical_data$time) * 0.8)
    x_max <- max(time_max, max(empirical_data$time) * 1.2)

    x_breaks <- scales::pretty_breaks(n = 6)(c(x_min, x_max))
    x_breaks <- x_breaks[x_breaks > 0]

    # Create the plot with comprehensive warning suppression
    p <- ggplot2::ggplot() +
      # Add confidence ribbon with warning suppression
      suppressWarnings(ggplot2::geom_ribbon(
        data = theoretical_data,
        ggplot2::aes(x = time, ymin = log_neg_log_lower, ymax = log_neg_log_upper),
        alpha = alpha_ribbon,
        fill = "lightblue"
      )) +
      # Add fitted line with warning suppression
      suppressWarnings(ggplot2::geom_line(
        data = theoretical_data,
        ggplot2::aes(x = time, y = log_neg_log_fitted),
        color = "blue",
        size = line_size
      )) +
      # Add empirical points (events)
      suppressWarnings(ggplot2::geom_point(
        data = event_data,
        ggplot2::aes(x = time, y = log_neg_log_surv, shape = type),
        size = point_size + 1,
        color = "black",
        fill = "white"
      )) +
      # Add censored observations
      {if (nrow(censored_data) > 0) {
        suppressWarnings(ggplot2::geom_point(
          data = censored_data,
          ggplot2::aes(x = SURVTIME, y = log_neg_log_surv, shape = type),
          size = 4,
          color = "black"
        ))
      }} +
      # Set shapes manually for legend
      ggplot2::scale_shape_manual(
        name = "",
        values = c("Events" = 21, "Right Censored" = 124),
        guide = ggplot2::guide_legend(override.aes = list(size = c(4, 5)))
      ) +
      # Use log scale for x-axis with data-driven range
      ggplot2::scale_x_log10(
        limits = c(x_min, x_max),
        breaks = x_breaks,
        labels = as.character(x_breaks)
      ) +
      # Set y-axis with survival probability labels
      ggplot2::scale_y_continuous(
        name = "Survival Probability (%)",
        limits = c(y_min, y_max),
        expand = c(0, 0),
        breaks = final_breaks,
        labels = final_labels
      ) +
      # Customize axes and labels
      ggplot2::labs(
        title = plot_title,
        x = paste0("Time to Event (", time_scale, ")")
      ) +
      # Add theme
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 9),
        panel.grid.minor = ggplot2::element_line(color = "grey90", size = 0.3),
        panel.grid.major = ggplot2::element_line(color = "grey80", size = 0.5),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1),
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.box = "horizontal",
        legend.margin = ggplot2::margin(t = 10),
        legend.box.background = ggplot2::element_rect(fill = "white", color = "black", size = 0.5),
        legend.box.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10)
      )

    # Calculate dynamic box size based on text content
    text_lines <- strsplit(legend_text, "\n")[[1]]
    max_chars <- max(nchar(text_lines))  # Longest line character count
    num_lines <- length(text_lines)      # Number of lines

    # Calculate box dimensions (with padding)
    char_width <- 0.07  # Approximate width per character in inches for mono font
    line_height <- 0.12 # Approximate height per line in inches
    padding_width <- 0.3  # Left/right padding
    padding_height <- 0.2 # Top/bottom padding

    box_width <- max_chars * char_width + padding_width
    box_height <- num_lines * line_height + padding_height

    # Add statistics box with dynamic white background and black border
    p <- p +
      ggplot2::annotation_custom(
        grob = grid::rectGrob(
          x = 0.02, y = 0.98,
          width = grid::unit(box_width, "inches"),
          height = grid::unit(box_height, "inches"),
          just = c("left", "top"),
          gp = grid::gpar(fill = "white", col = "black", lwd = 1)
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
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
