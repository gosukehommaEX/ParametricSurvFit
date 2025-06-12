# ParametricSurvFit

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/ParametricSurvFit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/ParametricSurvFit/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/ParametricSurvFit/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/ParametricSurvFit?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/ParametricSurvFit)](https://CRAN.R-project.org/package=ParametricSurvFit)
<!-- badges: end -->

## Overview

`ParametricSurvFit` is a comprehensive R package designed for fitting multiple parametric survival distributions to clinical trial data. It provides streamlined functions for analyzing ADaM datasets commonly used in pharmaceutical research, with support for stratified analysis, covariance matrix extraction, and formatted output tables.

## Features

- **Multiple Distribution Support**: Fit exponential, Weibull, log-normal, log-logistic, Gompertz, generalized gamma, and gamma distributions
- **ADaM Dataset Integration**: Direct support for ADSL and ADTTE datasets
- **Stratified Analysis**: Flexible stratification by any categorical variable
- **Covariance Matrix Extraction**: Extract variance-covariance matrices from fitted models
- **Formatted Output**: Professional tables with kableExtra formatting
- **Parameter Extraction**: Comprehensive extraction of distribution parameters with confidence intervals
- **Population Filtering**: Support for ITT, Safety, and Randomized populations

## Installation

You can install the development version of ParametricSurvFit from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ParametricSurvFit
devtools::install_github("gosukehommaEX/ParametricSurvFit")
```

## Dependencies

The package requires the following R packages:
- `dplyr` (>= 1.0.0)
- `purrr` (>= 0.3.0)
- `tidyr` (>= 1.1.0)
- `flexsurv` (>= 2.0)
- `survival` (>= 3.2.0)
- `haven` (>= 2.4.0)
- `kableExtra` (>= 1.3.0)
- `rlang` (>= 0.4.0)
- `stats` (>= 3.5.0)

## Quick Start

```r
library(ParametricSurvFit)

# 1. Create survival dataset from ADaM files
surv_data <- DataParametricSurv(
  adsl_path = "path/to/adsl.sas7bdat",
  adtte_path = "path/to/adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = "SEX"
)

# 2. Fit multiple parametric distributions
results <- FitSurvMods(
  dataset = surv_data,
  distributions = c("exp", "weibull", "lnorm", "llogis"),
  table_caption = "Parametric Models for Overall Survival"
)

# 3. Extract covariance matrix for a specific distribution
cov_matrix <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "weibull",
  stratify_reference_level = "M",  # Use Male as reference
  create_kable = TRUE
)

# 4. View formatted results
print(results)
print(cov_matrix)
```

## Main Functions

### `DataParametricSurv()`
Creates a standardized dataset for parametric survival analysis from ADaM datasets.

### `FitSurvMods()`
Fits multiple parametric survival distributions and creates formatted output tables.

### `ExtractParams()`
Extracts parameters, standard errors, and confidence intervals from fitted models.

### `ExtractCovarianceMatrix()` 🆕
Extracts variance-covariance matrices from fitted parametric survival models.

**Key Parameters:**
- `dataset`: Output from `DataParametricSurv()`
- `distribution`: Single distribution name to fit
- `include_shape`: Include shape parameters in covariance matrix (for applicable distributions)
- `create_kable`: Return formatted kable output
- `stratify_reference_level`: Specify reference level for stratification (e.g., "M" or "F" for SEX)

## Supported Distributions

| Distribution | Parameters | Description |
|--------------|------------|-------------|
| `exp` | Rate, Intercept | Exponential distribution |
| `weibull` | Shape, Scale, Intercept | Weibull distribution |
| `lnorm` | Intercept (meanlog), Scale (sdlog) | Log-normal distribution |
| `llogis` | Intercept, Scale | Log-logistic distribution |
| `gompertz` | Intercept, Rate, Shape | Gompertz distribution |
| `gengamma` | Intercept (mu), Scale (sigma), Shape (Q) | Generalized gamma distribution |
| `gamma` | Intercept, Shape, Rate | Gamma distribution |

## Example Workflows

### Basic Analysis
```r
# Simple overall survival analysis without stratification
surv_data <- DataParametricSurv(
  adsl_path = "adsl.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = NULL  # No stratification
)

results <- FitSurvMods(surv_data, distributions = c("exp", "weibull"))
```

### Stratified Analysis
```r
# Analysis stratified by sex
surv_data_sex <- DataParametricSurv(
  adsl_path = "adsl.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = "SEX"
)

results_sex <- FitSurvMods(
  dataset = surv_data_sex,
  distributions = c("exp", "weibull", "gamma"),
  table_caption = "OS Analysis by Sex",
  stratify_name = "SEX"
)
```

### Covariance Matrix Analysis
```r
# Extract covariance matrix for Weibull distribution
cov_weibull <- ExtractCovarianceMatrix(
  dataset = surv_data_sex,
  distribution = "weibull",
  parameter_name = "Overall Survival",
  stratify_reference_level = "M",  # Use Male as reference
  create_kable = TRUE
)

# Include shape parameter for gamma distribution
cov_gamma_with_shape <- ExtractCovarianceMatrix(
  dataset = surv_data_sex,
  distribution = "gamma",
  include_shape = TRUE,
  stratify_reference_level = "F"  # Use Female as reference
)
```

## Covariance Matrix Features

The `ExtractCovarianceMatrix()` function provides:
- **Distribution-specific handling**: Adapts to each distribution's parameter structure
- **Shape parameter control**: Optional inclusion of shape parameters
- **Reference level specification**: Control stratification encoding (e.g., "M" vs "F" for SEX)
- **Formatted output**: Professional kable tables with merged cells
- **Raw matrix option**: For advanced statistical analysis

### Distributions with Shape Parameters
- **Gamma**: `include_shape = TRUE` includes the shape parameter
- **Generalized Gamma**: `include_shape = TRUE` includes the Q parameter
- **Gompertz**: `include_shape = TRUE` includes the shape parameter

### Common Stratification Reference Levels
- **SEX**: Use `"M"` (Male) or `"F"` (Female)
- **REGION**: Use specific region codes like `"US"`, `"EU"`, etc.
- **AGE_GROUP**: Use `"<65"`, `">=65"`, etc.

## Output Format

The package provides professionally formatted tables with:
- Grouped headers by distribution
- Parameter estimates with standard errors
- 95% confidence intervals
- Stratified results when applicable
- Interactive hover effects
- Exportable formats (HTML, PDF, etc.)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Homma, G. (2025). ParametricSurvFit: Parametric Survival Model Fitting and Analysis. 
R package version 0.1.0. https://github.com/gosukehommaEX/ParametricSurvFit
```

## Acknowledgments

- Built using the excellent [flexsurv](https://github.com/chjackson/flexsurv) package
- Table formatting powered by [kableExtra](https://github.com/haozhu233/kableExtra)
- Data manipulation using the [tidyverse](https://www.tidyverse.org/) ecosystem
