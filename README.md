# DiPPER

This is an R package implementing DiPPER (**Di**fferential **P**revalence via
**P**robabilistic **E**stimation in **R**).

DiPPER is a Bayesian hierarchical model-based approach designed for differential
prevalence analysis, particularly in microbiome studies.

A pre-print of the paper introducing DiPPER can be found
[here](https://arxiv.org/abs/2602.05938).

## Advantages of DiPPER

* **Inherent multiplicity adjustment:** DiPPER provides differential
    prevalence estimates and uncertainty intervals that are *naturally* and
    automatically adjusted for multiple testing through hierarchical shrinkage.
* **Robust in boundary cases:** DiPPER produces finite differential prevalence
    estimates and uncertainty intervals even in boundary cases (e.g. when a taxon
    is completely absent or fully present in one study group) where some
    standard approaches fail.

## Installation

You can install the development version of DiPPER from GitHub using:

```r
# install.packages("remotes")
remotes::install_github("jepelt/DiPPER")
```

DiPPER also requires the cmdstanr package and a working Stan toolchain.
You can install them as follows:

```r
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev/",
                           getOption("repos")))

# Set up the C++ toolchain (Windows users may need Rtools)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Install the Stan backend (only needs to be done once)
cmdstanr::install_cmdstan()
```

## Example Usage

Below is a simple example workflow using the example data included in the
package.

```r
library(DiPPER)

# Load example data (TreeSummarizedExperiment object)
# This dataset compares (N = 20 + 20) rats on a High/Low fat diet.
data("tse_hintikka")

# Run DiPPER. 
# The first term in the formula (here: Fat) is automatically 
# used as the variable of interest. XOS (xylo-oligosaccharide supplementation)
# is included as a covariate to adjust for.
# Note: When using DiPPER for the first time, it may take around two minutes
# to run the function due to the compilation of the Stan model.

fit <- DiPPER(
  tse = tse_hintikka,
  formula = ~ Fat + XOS,
  tax_rank = "Genus",
  seed = 1
)

# Extract summarized results as a data.frame
res <- summary(fit)

# Create a forest plot (showing only 'significant' taxa)
plot(fit, show_taxa = "significant")
```
