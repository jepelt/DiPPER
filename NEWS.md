# DiPPER 0.99.1

* `data.type` is now a required argument with no default. A new value,
  `"abundance"`, covers non-count, non-proportional abundance data
  such as CPM or TPM values.
* `read.depth` must now be chosen explicitly for count data. For other data
  types, it defaults to `FALSE`. Moreover, read depths supplied as a
  metadata column are now log10-transformed and added to the model
  automatically.
* `print()` now reports also the estimated values of the hyperparameters of the
  hierarchical Laplace prior. For this to work `keep.pars` of `dipper()` now
  defaults to `c("beta", "tau", "nu")`.
* `verbose` argument is added to `dipper()`. Setting `verbose = FALSE`
  silences all messages, including those from Stan.
* `threshold` is now always interpreted on the scale of the supplied data.
  For instance, for `data.type = "counts"`, `threshold = 2` means that count
  values of 2 or below are interpreted as absence.
* The vignette and the details section of `?dipper` have been updated regarding
  the type of the input data, controlling for sequencing depth, and the
  choice between the symmetric and the asymmetric hierarchical prior.
* Bugs in parsing model formula have been fixed.

# DiPPER 0.99.0

* Initial Bioconductor submission.
* Provides an implementation of DiPPER (Differential Prevalence via
Probabilistic Estimation in R), a Bayesian hierarchical modeling approach for
differential prevalence analysis, especially in microbiome studies.
* Supports both cross-sectional and longitudinal (repeated measures) designs,
  with covariate adjustment and automatic control for sequencing depth.
* Posterior inference is performed with CmdStan, which must be installed
  separately before installing DiPPER.
