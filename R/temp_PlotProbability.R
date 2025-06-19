#' Create Probability Plot with log(-log(S)) transformation for ARM*STRATIFY Groups
#'
#' This function creates probability plots with log(-log(S)) transformation for
#' parametric survival models fitted using flexsurv. It creates separate plots
#' for each ARM*STRATIFY combination following SAS LIFEREG PROBPLOT procedure.
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
#' @return A list of ggplot2 objects containing probability plots for each ARM*STRATIFY combination
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_segment
#'   labs theme_minimal theme element_text scale_x_log10 scale_y_continuous
#'   annotation_custom unit
#' @importFrom survival survfit Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom dplyr filter mutate arrange group_by
#' @importFrom stats qnorm pexp pweibull plnorm pgamma vcov
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
#' # Create probability plots with Weibull distribution
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

  # Set stratify_name if not provided
  if (is.null(stratify_name)) {
    stratify_name <- "STRATIFY"
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
        cdf = 1 - surv,
        log_neg_log_surv = log(-log(surv))  # y-axis transformation: log(-log(S))
      ) %>%
      dplyr::filter(surv > 0.001 & surv < 0.999, is.finite(log_neg_log_surv)) %>%
      dplyr::arrange(time)

    if (nrow(empirical_data) == 0) {
      warning("No valid empirical data points for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    # Create time sequence for model predictions - extend range to cover full plot area
    time_seq <- seq(0.1, 100, length.out = 300)  # Cover full x-axis range

    # Get model predictions with confidence intervals
    model_pred <- tryCatch({
      summary(fit_model, t = time_seq, type = "survival", ci = TRUE)
    }, error = function(e) {
      # Fallback to basic prediction without CI
      basic_pred <- summary(fit_model, t = time_seq, type = "survival")
      list(basic_pred)
    })

    # Handle model prediction format
    if (is.list(model_pred) && length(model_pred) > 0) {
      pred_df <- model_pred[[1]]

      # Check if confidence intervals are available
      if (all(c("est", "lcl", "ucl") %in% names(pred_df))) {
        theoretical_data <- data.frame(
          time = time_seq,
          fitted_surv = pred_df$est,
          lower_surv = pred_df$lcl,
          upper_surv = pred_df$ucl
        )
        has_ci <- TRUE
      } else {
        theoretical_data <- data.frame(
          time = time_seq,
          fitted_surv = pred_df$est
        )
        has_ci <- FALSE
      }
    } else {
      theoretical_data <- data.frame(
        time = time_seq,
        fitted_surv = model_pred$est
      )
      has_ci <- FALSE
    }

    # Apply transformations (keep all data - don't filter by survival bounds)
    theoretical_data <- theoretical_data %>%
      dplyr::mutate(log_neg_log_fitted = log(-log(pmax(0.001, pmin(0.999, fitted_surv)))))

    # Add CI transformations if available
    if (has_ci && all(c("lower_surv", "upper_surv") %in% names(theoretical_data))) {
      theoretical_data <- theoretical_data %>%
        dplyr::mutate(
          # Constrain survival values to avoid log(0) issues but keep full range
          lower_surv_safe = pmax(0.001, pmin(0.999, lower_surv)),
          upper_surv_safe = pmax(0.001, pmin(0.999, upper_surv)),
          log_neg_log_lower = log(-log(upper_surv_safe)),  # Note: CI bounds flip
          log_neg_log_upper = log(-log(lower_surv_safe))
        ) %>%
        dplyr::filter(
          is.finite(log_neg_log_fitted),
          is.finite(log_neg_log_lower),
          is.finite(log_neg_log_upper)
        ) %>%
        dplyr::select(-lower_surv_safe, -upper_surv_safe)  # Remove temporary columns
    } else {
      theoretical_data <- theoretical_data %>%
        dplyr::filter(is.finite(log_neg_log_fitted))
      has_ci <- FALSE
    }

    # Calculate y-axis range based on empirical data only
    empirical_y_values <- empirical_data$log_neg_log_surv[is.finite(empirical_data$log_neg_log_surv)]
    if (length(empirical_y_values) == 0) {
      warning("No finite y values for ARM: ", arm_val, ", STRATIFY: ", stratify_val)
      next
    }

    first_event_y <- min(empirical_y_values)
    y_min <- first_event_y  # No padding below first event
    y_max <- max(empirical_y_values) + 0.5

    # Calculate breaks and labels beforehand
    desired_percentiles <- c(1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99)
    desired_survival <- 1 - desired_percentiles / 100
    desired_survival[desired_survival <= 0] <- 0.001
    desired_survival[desired_survival >= 1] <- 0.999
    desired_breaks <- log(-log(desired_survival))

    # Filter to include breaks within range, but always include 1, 2, 5
    valid_indices <- which(desired_breaks >= y_min & desired_breaks <= y_max)
    first_three_indices <- 1:3  # Always include 1, 2, 5
    all_indices <- sort(unique(c(first_three_indices, valid_indices)))

    final_breaks <- desired_breaks[all_indices]
    final_labels <- as.character(desired_percentiles[all_indices])

    # Extract model parameters for legend
    n_total <- nrow(current_data)
    n_events <- sum(current_data$EVENT)
    n_censored <- n_total - n_events

    # Extract parameters using ExtractParams function for consistency
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
    if (has_ci) {
      plot_title <- paste0(
        dist_name_map[distribution], " Probability Plot With ", confidence_level * 100, "% Confidence Limits\n",
        "ARM: ", arm_val, ", ", stratify_name, ": ", stratify_val
      )
    } else {
      plot_title <- paste0(
        dist_name_map[distribution], " Probability Plot\n",
        "ARM: ", arm_val, ", ", stratify_name, ": ", stratify_val
      )
    }

    # Create event data for points
    event_data <- empirical_data %>%
      dplyr::mutate(type = "Events")

    # Create censored data positioned exactly at the y-axis minimum and within x-axis range
    censored_data <- current_data %>%
      dplyr::filter(EVENT == 0) %>%
      dplyr::filter(SURVTIME >= 0.1 & SURVTIME <= 100) %>%  # Filter to x-axis range
      dplyr::mutate(
        type = "Right Censored",
        log_neg_log_surv = y_min  # Position exactly at y-axis minimum
      )

    # Create the plot
    p <- suppressWarnings(ggplot2::ggplot() +
                            # Add confidence ribbon if available (suppress warnings for out-of-range values)
                            {if (has_ci) {
                              ggplot2::geom_ribbon(
                                data = theoretical_data,
                                ggplot2::aes(x = time, ymin = log_neg_log_lower, ymax = log_neg_log_upper),
                                alpha = alpha_ribbon,
                                fill = "lightblue"
                              )
                            }} +
                            # Add fitted line (suppress warnings for out-of-range values)
                            ggplot2::geom_line(
                              data = theoretical_data,
                              ggplot2::aes(x = time, y = log_neg_log_fitted),
                              color = "blue",
                              size = line_size
                            ) +
                            # Add empirical points (events) - larger size
                            ggplot2::geom_point(
                              data = event_data,
                              ggplot2::aes(x = time, y = log_neg_log_surv, shape = type),
                              size = point_size + 1,  # Increase size by 1
                              color = "black",
                              fill = "white"
                            ) +
                            # Add censored observations as vertical lines - larger size
                            {if (nrow(censored_data) > 0) {
                              ggplot2::geom_point(
                                data = censored_data,
                                ggplot2::aes(x = SURVTIME, y = log_neg_log_surv, shape = type),
                                size = 4,  # Increase from 3 to 4
                                color = "black"
                              )
                            }} +
                            # Set shapes manually for legend - adjust legend sizes too
                            ggplot2::scale_shape_manual(
                              name = "",
                              values = c("Events" = 21, "Right Censored" = 124),
                              guide = ggplot2::guide_legend(override.aes = list(size = c(4, 5)))  # Increase legend sizes
                            ) +
                            # Use log scale for x-axis with fixed range and clean labels
                            ggplot2::scale_x_log10(
                              limits = c(0.1, 100),  # Fixed range from 0.1 to 100
                              breaks = c(0.1, 1, 10, 100),
                              labels = c("0.1", "1", "10", "100")  # Remove .0 from integer values
                            ) +
                            # Set y-axis with pre-calculated breaks and labels and no clipping
                            ggplot2::scale_y_continuous(
                              name = "Percent",
                              limits = c(y_min, y_max),
                              expand = c(0, 0),
                              breaks = final_breaks,
                              labels = final_labels,
                              oob = scales::oob_keep  # Keep out-of-bounds values instead of removing them
                            ) +
                            # Customize axes and labels
                            ggplot2::labs(
                              title = plot_title,
                              x = paste0("Time to Event (", time_scale, ")")
                            ) +
                            # Add theme with black border
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
                            ))

    # Add statistics box
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
