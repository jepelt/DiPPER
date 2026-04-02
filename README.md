# DiPPER

This is an R package implementing DiPPER (Differential Prevalence via
Probabilistic Estimation in R).

DiPPER is a Bayesian hierarchical model designed
for differential prevalence analysis (especially) in microbiome studies.
Unlike standard frequentist approaches (based on e.g., the Wald test),
which may fail or yield infinite estimates in boundary cases (i.e., when a
taxon is completely absent in one group), DiPPER produces robust, finite
estimates through Bayesian regularization. Furthermore, the model provides
differential prevalence estimates and uncertainty intervals that are inherently
adjusted for multiplicity.

A pre-print of the paper introducing DiPPER can be found
[here](https://arxiv.org/abs/2602.05938).


## Installation

You can install the development version of DiPPER from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("jepelt/DiPPER")
```

Note: DiPPER requires cmdstanr and a working Stan backend to compile
and run the Bayesian models. You can install and configure it with:
R

```r
install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)

cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Check system setting
cmdstanr::install_cmdstan()
```

## Example Usage

Below is a simple example workflow using the example data included in the
package.

```r
library(DiPPER)

# Load example data (TreeSummarizedExperiment object)
data("tse_thomas")

# Run DiPPER
# The first term in the formula (here: study_condition) is automatically
# used as the main variable of interest.
fit <- DiPPER(
  tse = tse_thomas,
  formula = ~ study_condition + age + bmi + sex
)

# Extract the results
res <- summary_dipper(fit)

# Create a forest plot of the results (show only 'significant' results)
plot_dipper(fit, show_taxa = 'significant')
```
