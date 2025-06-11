# ParametricSurvFit

<!-- badges: start -->
[![R-CMD-check](https://github.com/yourusername/ParametricSurvFit/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/ParametricSurvFit/actions)
[![Codecov test coverage](https://codecov.io/gh/yourusername/ParametricSurvFit/branch/main/graph/badge.svg)](https://codecov.io/gh/yourusername/ParametricSurvFit?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/ParametricSurvFit)](https://CRAN.R-project.org/package=ParametricSurvFit)
<!-- badges: end -->

## Overview

`ParametricSurvFit` is a comprehensive R package designed for fitting multiple parametric survival distributions to clinical trial data. It provides streamlined functions for analyzing ADaM datasets commonly used in pharmaceutical research, with support for stratified analysis and formatted output tables.

## Features

- **Multiple Distribution Support**: Fit exponential, Weibull, log-normal, log-logistic, Gompertz, generalized gamma, and gamma distributions
- **ADaM Dataset Integration**: Direct support for ADSLPH and ADTTE datasets
- **Stratified Analysis**: Flexible stratification by any categorical variable
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
devtools::install_github("yourusername/ParametricSurvFit")
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

## Quick Start

```r
library(ParametricSurvFit)

# 1. Create survival dataset from ADaM files
surv_data <- DataParametricSurv(
  adslph_path = "path/to/adslph.sas7bdat",
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

# 3. View formatted results
print(results)
```

## Main Functions

### `DataParametricSurv()`
Creates a standardized dataset for parametric survival analysis from ADaM datasets.

**Key Parameters:**
- `adslph_path`: Path to ADSLPH dataset
- `adtte_path`: Path to ADTTE dataset  
- `population`: Analysis population ("ITTFL", "SAFFL", or "RANDFL")
- `variable`: Survival endpoint ("OS", "PFS", etc.)
- `stratify_by`: Stratification variable name or NULL

### `FitSurvMods()`
Fits multiple parametric survival distributions and creates formatted output tables.

**Key Parameters:**
- `dataset`: Output from `DataParametricSurv()`
- `distributions`: Vector of distribution names to fit
- `format_output`: Whether to return formatted table (TRUE) or raw data (FALSE)
- `table_caption`: Caption for the output table

### `ExtractParams()`
Extracts parameters, standard errors, and confidence intervals from fitted models.

**Key Parameters:**
- `fit_model`: Fitted flexsurv model object
- `distribution`: Distribution name

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
  adslph_path = "adslph.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = NULL  # No stratification
)

results <- FitSurvMods(surv_data, distributions = c("exp", "weibull"))
```

### Stratified Analysis
```r
# Analysis stratified by region
surv_data_region <- DataParametricSurv(
  adslph_path = "adslph.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "PFS",
  stratify_by = "REGION"
)

results_region <- FitSurvMods(
  dataset = surv_data_region,
  distributions = c("exp", "weibull", "gamma"),
  table_caption = "PFS Analysis by Region",
  stratify_name = "REGION"
)
```

### Custom Distribution Selection
```r
# Fit only specific distributions
results_custom <- FitSurvMods(
  dataset = surv_data,
  distributions = c("weibull", "lnorm", "gengamma"),
  format_output = FALSE  # Return raw data frame
)
```

## Output Format

The package provides professionally formatted tables with:
- Grouped headers by distribution
- Parameter estimates with standard errors
- 95% confidence intervals
- Stratified results when applicable
- Interactive hover effects
- Exportable formats (HTML, PDF, etc.)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Homma, G. (2025). ParametricSurvFit: Parametric Survival Model Fitting and Analysis. 
R package version 0.1.0. https://github.com/yourusername/ParametricSurvFit
```

## Contact

Gosuke Homma - my.name.is.gosuke@gmail.com

Project Link: [https://github.com/yourusername/ParametricSurvFit](https://github.com/yourusername/ParametricSurvFit)

## Acknowledgments

- Built using the excellent [flexsurv](https://github.com/chjackson/flexsurv) package
- Table formatting powered by [kableExtra](https://github.com/haozhu233/kableExtra)
- Data manipulation using the [tidyverse](https://www.tidyverse.org/) ecosystem
