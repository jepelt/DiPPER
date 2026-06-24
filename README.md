# DiPPER

This is an R package implementing DiPPER (**Di**fferential **P**revalence via
**P**robabilistic **E**stimation in **R**).

DiPPER is a Bayesian hierarchical model-based approach designed for differential
prevalence analysis, particularly in microbiome studies.

A pre-print of the paper introducing DiPPER can be found
[here](https://arxiv.org/abs/2602.05938).

## Advantages of DiPPER

* **Inherent multiplicity adjustment:** DiPPER provides differential
  prevalence estimates and uncertainty intervals that are naturally and
  automatically adjusted for multiple testing through hierarchical shrinkage.
  Differentially prevalent taxa can thus be detected without any post-hoc
  multiplicity adjustments.
* **Robust in boundary cases:** DiPPER produces finite differential prevalence
  estimates and uncertainty intervals even in boundary cases (e.g., when a taxon
  is completely absent or fully present in one study group) where some
  standard approaches fail.
* **Longitudinal data support:** DiPPER can handle longitudinal (or
  repeated measures) data by including random intercepts in the model.


## Installation

DiPPER relies on packages from Bioconductor (for handling microbiome data) and
Stan (for Bayesian inference). Please ensure these are installed first:

```r
# 1. Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("SummarizedExperiment", "TreeSummarizedExperiment"))

# 2. Install cmdstanr
install.packages("cmdstanr",
                 repos = c("[https://stan-dev.r-universe.dev/](https://stan-dev.r-universe.dev/)",
                           getOption("repos")))

# Set up the C++ toolchain (Windows users may need Rtools)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Install the Stan backend (only needs to be done once)
cmdstanr::install_cmdstan()
```

Once the dependencies are ready, you can install the development version of
DiPPER from GitHub using:

```r
# install.packages("remotes")
remotes::install_github("jepelt/DiPPER")
```

## Example Usage

### Example 1 (Basic cross-sectional setting)

Below is a simple example workflow using the built-in dataset, which compares
gut microbiomes of rats on a high-fat versus low-fat diet (N = 20 + 20) while
controlling for xylo-oligosaccharide supplementation. The dataset originates 
from [Hintikka et al. (2021)](https://www.mdpi.com/1660-4601/18/8/4049).

```r
library(DiPPER)

# Load dataset, provided as a SummarizedExperiment object, and agglomerated to 
# genus level.
data("tse_hintikka")

# Run DiPPER.
# The first term in the formula (here: Fat) is automatically used as the
# variable of interest. XOS (xylo-oligosaccharide supplementation) is included
# as a covariate. The name of the assay (here: 'counts') must also be provided.
#
# Note: When using DiPPER for the first time, compiling the Stan model may
# take a few minutes and produce verbose C++ output. This is expected.
fit <- dipper(
    tse = tse_hintikka,
    formula = ~ Fat + XOS,
    assay.type = "counts"
)

# Extract summarized results as a data.frame
res <- summary(fit)

# Create a forest plot (showing only 'significant' taxa)
# It can be seen how, for instance, Sellimonas genus has higher prevalence in
# High fat group (vs. Low fat group), while Faecalibaculum genus has higher 
# prevalence in Low fat group.
plot(fit, show.taxa = "significant")
```

<img src="man/figures/plot_hintikka.png" width="600" alt="Forest plot of the Hintikka dataset">

### Example 2 (Longitudinal setting)

This example demonstrates how to run DiPPER on repeated measures data by
including a random intercept for the subject. We use a built-in subset of
the [Vatanen et al. (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4950857/)
dataset. The dataset includes 1-3 gut microbiome samples from 40 infants 
during their first year of life.

The abundance data are provided as relative abundances (proportions), and
thus we cannot control for sequencing depth now.

```r
library(DiPPER)

# Load dataset, provided as a SummarizedExperiment object, and agglomerated to 
# genus level.
data("tse_vatanen")

# The metadata includes variables subject_id (indicating infant), age (in
# months), antibiotics (use of antibiotics: yes/no) and gender (male/female).
tse_vatanen |>
    SummarizedExperiment::colData() |>
    as.data.frame() |>
    head()

# Run DiPPER.
# By putting age as the first variable of the formula it is automatically set as
# the variable of interest. Gender and antibiotics use are controlled for.
# Lastly, the longitudinality of the data is accounted for by specifying
# subject_id as a random intercept term by writing (1 | subject_id).

# Note that because the abundance data are proportions, we must specify
# data.type = "relabundance" in addition to providing the assay name
# assay.type = "relative_abundance". Moreover, as we cannot now control for the
# sequencing depth, we must set read.depth = FALSE.
fit2 <- dipper(
    tse = tse_vatanen,
    formula = ~ age + antibiotics + gender + (1 | subject_id),
    assay.type = "relative_abundance",
    data.type = "relabundance",
    read.depth = FALSE,
    iter.sampling = 500 # Set to a lower value to speed up the example
)

# Extract the results for age 
res2 <- summary(fit2)

# Create a forest plot (showing 'significant' taxa only)
# As expected, the prevalence of several taxa increases with infant age.
plot(fit2, show.taxa = "significant")
```

<img src="man/figures/plot_vatanen.png" width="600" alt="Forest plot of the Vatanen dataset">
