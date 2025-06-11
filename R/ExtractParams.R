#' Extract Parameters from Fitted Parametric Survival Models
#'
#' This function extracts parameter estimates, standard errors, and 95% confidence
#' intervals from fitted parametric survival models using flexsurv package.
#'
#' @param fit_model A fitted survival model object from flexsurv::flexsurvreg()
#' @param distribution Character string specifying the distribution type.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#'
#' @return A list containing parameter estimates and confidence intervals with components Parameters, Estimates, SEs, LCLs, and UCLs
#'
#' @details
#' The function handles different parameterizations for each distribution:
#' \itemize{
#'   \item \strong{Exponential}: Rate and Intercept (log scale)
#'   \item \strong{Weibull}: Shape, Scale, and Intercept (log scale)
#'   \item \strong{Log-normal}: Intercept (meanlog) and Scale (sdlog)
#'   \item \strong{Log-logistic}: Intercept (log scale) and Scale
#'   \item \strong{Gompertz}: Intercept (log scale), Rate, and Shape
#'   \item \strong{Generalized Gamma}: Intercept (mu), Scale (sigma), and Shape (Q)
#'   \item \strong{Gamma}: Shape, Rate, and Intercept (log scale)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(flexsurv)
#' library(survival)
#'
#' # Fit a Weibull model
#' fit_weibull <- flexsurvreg(Surv(time, status) ~ 1,
#'                           data = lung,
#'                           dist = "weibull")
#'
#' # Extract parameters
#' params_weibull <- ExtractParams(fit_weibull, "weibull")
#' print(params_weibull)
#'
#' # Fit an exponential model
#' fit_exp <- flexsurvreg(Surv(time, status) ~ 1,
#'                       data = lung,
#'                       dist = "exp")
#'
#' # Extract parameters
#' params_exp <- ExtractParams(fit_exp, "exp")
#' print(params_exp)
#'
#' # Fit a gamma model
#' fit_gamma <- flexsurvreg(Surv(time, status) ~ 1,
#'                         data = lung,
#'                         dist = "gamma")
#'
#' # Extract parameters
#' params_gamma <- ExtractParams(fit_gamma, "gamma")
#' print(params_gamma)
#' }
ExtractParams <- function(fit_model, distribution) {

  # Validate distribution parameter
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!distribution %in% valid_distributions) {
    stop("Distribution must be one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Get basic model results
  res <- fit_model[['res']]

  if (distribution == 'exp') {
    # Exponential: Scale, Rate and Intercept
    Scale.est <- res[1,  'est']
    Scale.SE  <- res[1,   'se']
    Scale.LCL <- res[1, 'L95%']
    Scale.UCL <- res[1, 'U95%']

    Rate.est <- 1 / Scale.est
    Rate.SE  <- Scale.SE / (Scale.est ^ 2)
    Rate.LCL <- 1 / Scale.UCL
    Rate.UCL <- 1 / Scale.LCL

    # Intercept (log scale)
    Intercept.est <- log(Rate.est)
    Intercept.SE  <- Rate.SE / Rate.est
    Intercept.LCL <- log(Rate.LCL)
    Intercept.UCL <- log(Rate.UCL)

    return(
      list(
        Parameters = c('Intercept', 'Rate'),
        Estimates  = c(Intercept.est, Rate.est),
        SEs  = c(Intercept.SE, Rate.SE),
        LCLs = c(Intercept.LCL, Rate.LCL),
        UCLs = c(Intercept.UCL, Rate.UCL)
      )
    )

  } else if (distribution == 'weibull') {
    # Weibull: Shape, Scale and Intercept
    Shape.est <- res['shape', 'est']
    Shape.SE  <- res['shape', 'se']
    Shape.LCL <- res['shape', 'L95%']
    Shape.UCL <- res['shape', 'U95%']

    Scale.est <- res['scale', 'est']
    Scale.SE  <- res['scale', 'se']
    Scale.LCL <- res['scale', 'L95%']
    Scale.UCL <- res['scale', 'U95%']

    # Intercept (log scale)
    Intercept.est <- log(Scale.est)
    Intercept.SE  <- Scale.SE / Scale.est
    Intercept.LCL <- log(Scale.LCL)
    Intercept.UCL <- log(Scale.UCL)

    return(
      list(
        Parameters = c('Intercept', 'Scale', 'Shape'),
        Estimates = c(Intercept.est, Scale.est, Shape.est),
        SEs  = c(Intercept.SE, Scale.SE, Shape.SE),
        LCLs = c(Intercept.LCL, Scale.LCL, Shape.LCL),
        UCLs = c(Intercept.UCL, Scale.UCL, Shape.UCL)
      )
    )

  } else if (distribution == 'lnorm') {
    # Log-normal: Intercept (meanlog) and Scale (sdlog)
    Meanlog.est <- res['meanlog', 'est']
    Meanlog.SE  <- res['meanlog', 'se']
    Meanlog.LCL <- res['meanlog', 'L95%']
    Meanlog.UCL <- res['meanlog', 'U95%']

    sdlog.est <- res['sdlog', 'est']
    sdlog.SE  <- res['sdlog', 'se']
    sdlog.LCL <- res['sdlog', 'L95%']
    sdlog.UCL <- res['sdlog', 'U95%']

    return(
      list(
        Parameters = c('Intercept', 'Scale'),
        Estimates = c(Meanlog.est, sdlog.est),
        SEs  = c(Meanlog.SE, sdlog.SE),
        LCLs = c(Meanlog.LCL, sdlog.LCL),
        UCLs = c(Meanlog.UCL, sdlog.UCL)
      )
    )

  } else if (distribution == 'llogis') {
    # Log-logistic: Shape, Scale and Intercept
    Shape.est <- res['shape', 'est']
    Shape.SE  <- res['shape', 'se']
    Shape.LCL <- res['shape', 'L95%']
    Shape.UCL <- res['shape', 'U95%']

    Scale.est <- 1 / Shape.est
    Scale.SE  <- Shape.SE / (Shape.est ^ 2)
    Scale.LCL <- 1 / Shape.UCL
    Scale.UCL <- 1 / Shape.LCL

    alpha.est <- res['scale', 'est']
    alpha.SE  <- res['scale', 'se']
    alpha.LCL <- res['scale', 'L95%']
    alpha.UCL <- res['scale', 'U95%']

    Intercept.est <- log(alpha.est)
    Intercept.SE  <- alpha.SE / alpha.est
    Intercept.LCL <- log(alpha.LCL)
    Intercept.UCL <- log(alpha.UCL)

    return(
      list(
        Parameters = c('Intercept', 'Scale'),
        Estimates = c(Intercept.est, Scale.est),
        SEs = c(Intercept.SE, Scale.SE),
        LCLs = c(Intercept.LCL, Scale.LCL),
        UCLs = c(Intercept.UCL, Scale.UCL)
      )
    )

  } else if (distribution == 'gompertz') {
    # Gompertz: Shape, Rate and Intercept (log scale)
    Shape.est <- res['shape', 'est']
    Shape.SE  <- res['shape', 'se']
    Shape.LCL <- res['shape', 'L95%']
    Shape.UCL <- res['shape', 'U95%']

    Rate.est <- res['rate', 'est']
    Rate.SE  <- res['rate', 'se']
    Rate.LCL <- res['rate', 'L95%']
    Rate.UCL <- res['rate', 'U95%']

    # Intercept (log scale)
    Intercept.est <- log(Rate.est)
    Intercept.SE  <- Rate.SE / Rate.est
    Intercept.LCL <- log(Rate.LCL)
    Intercept.UCL <- log(Rate.UCL)

    return(
      list(
        Parameters = c('Intercept', 'Rate', 'Shape'),
        Estimates = c(Intercept.est, Rate.est, Shape.est),
        SEs = c(Intercept.SE, Rate.SE, Shape.SE),
        LCLs = c(Intercept.LCL, Rate.LCL, Shape.LCL),
        UCLs = c(Intercept.UCL, Rate.UCL, Shape.UCL)
      )
    )

  } else if (distribution == 'gengamma') {
    # Generalized Gamma: Scale (sigma), Shape (Q) and Intercept (mu)
    Scale.est <- res['sigma', 'est']
    Scale.SE  <- res['sigma', 'se']
    Scale.LCL <- res['sigma', 'L95%']
    Scale.UCL <- res['sigma', 'U95%']

    Shape.est <- res['Q', 'est']
    Shape.SE  <- res['Q', 'se']
    Shape.LCL <- res['Q', 'L95%']
    Shape.UCL <- res['Q', 'U95%']

    Intercept.est <- res['mu', 'est']
    Intercept.SE  <- res['mu', 'se']
    Intercept.LCL <- res['mu', 'L95%']
    Intercept.UCL <- res['mu', 'U95%']

    return(
      list(
        Parameters = c('Intercept', 'Scale', 'Shape'),
        Estimates = c(Intercept.est, Scale.est, Shape.est),
        SEs = c(Intercept.SE, Scale.SE, Shape.SE),
        LCLs = c(Intercept.LCL, Scale.LCL, Shape.LCL),
        UCLs = c(Intercept.UCL, Scale.UCL, Shape.UCL)
      )
    )

  } else if (distribution == 'gamma') {
    # Gamma: Shape, Rate and Intercept (log scale)
    Shape.est <- res['shape', 'est']
    Shape.SE  <- res['shape', 'se']
    Shape.LCL <- res['shape', 'L95%']
    Shape.UCL <- res['shape', 'U95%']

    Rate.est <- res['rate', 'est']
    Rate.SE  <- res['rate', 'se']
    Rate.LCL <- res['rate', 'L95%']
    Rate.UCL <- res['rate', 'U95%']

    # Intercept (log scale)
    Intercept.est <- log(Rate.est)
    Intercept.SE  <- Rate.SE / Rate.est
    Intercept.LCL <- log(Rate.LCL)
    Intercept.UCL <- log(Rate.UCL)

    return(
      list(
        Parameters = c('Intercept', 'Shape', 'Rate'),
        Estimates = c(Intercept.est, Shape.est, Rate.est),
        SEs = c(Intercept.SE, Shape.SE, Rate.SE),
        LCLs = c(Intercept.LCL, Shape.LCL, Rate.LCL),
        UCLs = c(Intercept.UCL, Shape.UCL, Rate.UCL)
      )
    )
  }
}
