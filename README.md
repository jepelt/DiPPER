# DiPPER

This is an R package implementing DiPPER (**Di**fferential **P**revalence via 
**P**robabilistic **E**stimation in **R**).

DiPPER is a Bayesian hierarchical model designed for differential prevalence 
analysis, particularly in microbiome studies. Unlike standard frequentist 
approaches (e.g., Wald test) which may fail or yield infinite estimates in 
boundary cases (such as when a taxon is completely absent in one group), 
DiPPER produces robust, finite estimates through Bayesian regularization. 
Furthermore, the model provides differential prevalence estimates and 
uncertainty intervals that are inherently adjusted for multiplicity.

A pre-print of the paper introducing DiPPER can be found 
[here](https://arxiv.org/abs/2602.05938).

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
# This dataset compares (N = 20 + 20) rats on High/Low fat diet.
data("tse_hintikka")

# Run DiPPER. 
# The first term in the formula (here: Fat) is automatically 
# used as the variable of interest. XOS (xylo-oligosaccaride supplementation)
# is included as a covariate to adjust for.
# Note: This may take around two minutes to run for the first time due to
# compilation.

fit <- DiPPER(
  tse = tse_hintikka,
  formula = ~ Fat + XOS,
  tax_rank = 'Genus',
  seed = 1
)

# Extract summarized results as a data.frame
res <- summary(fit)

# Create a forest plot (showing only 'significant' taxa)
plot(fit, show_taxa = "significant")
```
