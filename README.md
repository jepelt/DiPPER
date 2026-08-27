# DiPPER

This is an R package implementing DiPPER (**Di**fferential **P**revalence via
**P**robabilistic **E**stimation in **R**).

DiPPER is a Bayesian hierarchical model-based approach designed for
*differential prevalence analysis*, that is, for analyzing
associations between the *presence/absence status* of taxonomic features 
(e.g. microbial species or genera) and an external variable (e.g. host 
disease status).

Alongside the variable of interest (e.g. disease status), DiPPER can adjust for
*covariates* (e.g. age, sex, BMI) and it automatically accounts for sequencing
depth by default. DiPPER can also be used in *longitudinal* (repeated measures)
designs.

## Advantages of DiPPER

* **Robust in boundary cases:** DiPPER produces finite differential prevalence
  estimates and uncertainty intervals even in boundary cases (e.g., when a
  taxon is completely absent or fully present in one study group) where some
  standard approaches fail.
* **Longitudinal data support:** DiPPER can handle longitudinal (or
  repeated measures) data by including random intercepts in the model.
* **Inherent multiplicity adjustment:** DiPPER provides differential
  prevalence estimates and uncertainty intervals that are naturally and
  automatically adjusted for multiple testing through hierarchical shrinkage.
  Differentially prevalent taxa can thus be detected without any post-hoc
  multiplicity adjustments.
* **High sensitivity:** When benchmarked against alternative methods, DiPPER
  detected potentially differentially prevalent taxa with high sensitivity
  while maintaining the error rate at the nominal level
  (see [the pre-print](https://arxiv.org/abs/2602.05938) for details).

## Installation

**1. Install the CmdStan backend**

DiPPER relies on CmdStan for Bayesian inference, so you first need to install
the `cmdstanr` package and configure the C++ toolchain.

```r
# Install cmdstanr from the official stan-dev repository
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           getOption("repos")))

# Set up the C++ toolchain (Windows users may be prompted to install Rtools)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Install the actual CmdStan backend (only needs to be done once)
cmdstanr::install_cmdstan()
```

**2. Install DiPPER**

Once CmdStan is ready, install DiPPER. Using `BiocManager` automatically
installs all required Bioconductor dependencies as well.

```r
# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install the development version from GitHub
BiocManager::install("jepelt/DiPPER")
```

DiPPER compiles its Stan models during installation, which takes a few minutes
and produces verbose output. This is expected, and it is done only once.

## Example

Below is a minimal workflow using the built-in `tse_hintikka` dataset, which
compares gut microbiomes of rats on a high-fat versus low-fat diet
(N = 20 + 20) while controlling for xylo-oligosaccharide supplementation. The
dataset originates from
[Hintikka et al. (2021)](https://www.mdpi.com/1660-4601/18/8/4049).

```r
library(DiPPER)

# Load the dataset (a TreeSummarizedExperiment agglomerated to genus level)
data("tse_hintikka")

# Run DiPPER.
# The first term in the formula (here: Fat) is automatically used as the
# variable of interest. XOS (xylo-oligosaccharide supplementation) is included
# as a covariate. The name of the assay (here: 'counts') must also be provided.
fit <- dipper(
    tse = tse_hintikka,
    formula = ~ Fat + XOS,
    assay.type = "counts",
    niter = 400 # This needs to be increased in practice!
)

# Extract summarized results as a data.frame
res <- summary(fit)

# Create a forest plot (showing only 'significant' taxa).
plot(fit)
```

<img src="man/figures/plot_hintikka.png" width="600" alt="Forest plot of the Hintikka dataset">

## Documentation

The package vignette covers the above example in more detail, along with an
example using longitudinal data and a description of the summary output. After
installing DiPPER, open it with:

```r
vignette("DiPPER")
```
