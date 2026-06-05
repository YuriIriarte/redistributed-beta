# RB: An R package for modeling unit data using the Redistributed-Beta distribution

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20550857.svg)](https://doi.org/10.5281/zenodo.20550857)

R package for modeling data on the unit interval $(0,1)$ using the Redistributed-Beta (RB) distribution.

## Overview

The **RB** package implements the Redistributed-Beta distribution, a flexible model for bounded data on the unit interval $(0,1)$.

The package provides tools for:

- Density, distribution function, random generation, and raw moment computation
- Moment-based and maximum likelihood estimation
- Descriptive measures, including mean, variance, coefficient of variation, skewness, and kurtosis
- Model comparison with beta and Kumaraswamy distributions
- Bootstrap goodness-of-fit procedures
- Reproducibility scripts for simulations, figures, and real-data applications

All scripts used in the associated manuscript are included to support computational reproducibility.

## Installation

You can install the stable version used in the paper directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("YuriIriarte/redistributed-beta")
```

## Basic usage

```markdown
## Basic Usage

```r
library(RB)

# Generate data from the Redistributed-Beta distribution
set.seed(123)
x <- rRB(n = 100, alpha = 2, beta = 3)

# Evaluate the density and distribution function
dRB(0.4, alpha = 2, beta = 3)
pRB(0.4, alpha = 2, beta = 3)

# Compute raw moments and summary measures
mRB(1:4, alpha = 2, beta = 3)
meanRB(alpha = 2, beta = 3)
varRB(alpha = 2, beta = 3)
skewRB(alpha = 2, beta = 3)
kurtRB(alpha = 2, beta = 3)

# Fit the model by maximum likelihood
fit <- fitRB_mle(x)

# Summarize the fitted model
summary(fit)
```

## Reproducibility Resources

All scripts used in the simulation studies, figures, and real-data
applications reported in the paper *A Support-Preserving Redistribution Mechanism for Constructing 
Flexible Bounded Distributions* are distributed within the package.

After installing and loading the package, the available scripts can
be listed using:

```r
list.files(system.file("paper", package = "RB"))
```

For example, the Monte Carlo simulation script can be accessed
directly through:

```r
file.edit(system.file("paper/simulations.R", package = "RB"))
```

Additional scripts for figures and applications can be accessed in
the same way.

Documentation for all implemented functions is available through the
standard R help system. For example:

```r
?dCB
?fitCB_mle
?gofCB_boot
```

## Citation

If you use this package, please cite:

```r
@Manual{iriarte2026,
  title = {{RB}: R package for modeling unit data with the Redistributed-Beta distribution},
  author = {Yuri A. Iriarte and Juan M. Astorga},
  year = {2026},
  note = {R package version 1.0.0},
  doi = {10.5281/zenodo.20550857},
  url = {https://doi.org/10.5281/zenodo.20550857}
}
```

Iriarte, Y. A., and Astorga, J. M. (2026). RB: R package for modeling unit data with the Redistributed-Beta distribution. 
R package version 1.0.0. DOI: https://doi.org/10.5281/zenodo.20550857.
