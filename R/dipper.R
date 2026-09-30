#' DiPPER (Differential Prevalence via Probabilistic Estimation in R)
#'
#' This is the main wrapper function for running DiPPER. It first prepares the
#' input data and then fits the model using CmdStanR.
#'
#' @inheritParams prep_dipper_data
#' @inheritParams run_dipper
#' @param ... Additional arguments passed directly to \code{run_dipper}
#'   (e.g., cmdstanr specific arguments).
#'
#' @details
#' DiPPER is a Bayesian hierarchical model-based approach designed for
#' differential prevalence analysis, especially for microbiome studies.
#' It is designed for study designs where there is one
#' \strong{variable of interest} (e.g. treatment group, disease status, etc.)
#' and potentially other covariates (e.g. age, sex, etc.) that may confound
#' the associations between microbe prevalences and the variable of interest.
#' DiPPER can also be applied to longitudinal or repeated measures data.
#'
#' Technically, DiPPER models the presence/absence of taxonomic features (e.g.
#' genera or species) using multiple logistic regression models. The models are
#' connected by a common (Asymmetric) Laplace prior for the variable of
#' interest. The hierarchical structure guarantees robust (and always finite)
#' differential prevalence estimates and uncertainty intervals that are also
#' effectively multiplicity-adjusted.
#'
#' \subsection{Input data and preprocessing}{
#' The \code{dipper} function takes either a \code{(Tree)SummarizedExperiment}
#' object or an abundance matrix and a metadata \code{data.frame} as input.
#' The abundance matrix can be sequencing counts, relative abundances, other
#' non-negative abundances (e.g. pathway abundances), or presence/absence
#' data. The \code{data.type} argument is required, as it determines how
#' presence is defined and whether read depth can be derived from the data.
#'
#' Before model fitting, the input data undergoes automated preprocessing.
#' Continuous abundance data are converted to presence/absence format based on
#' the specified \code{threshold}. By default, a feature is absent only if its
#' value is exactly zero. Furthermore, prevalence filtering is applied: taxa
#' that are not present in at least \code{min.present} samples, or absent in at
#' least \code{min.absent} samples, are excluded from the analysis.
#' }
#'
#' \subsection{Model formula and the variable of interest}{
#' The \code{dipper} function uses the standard \code{formula} argument to
#' specify the model formula. A possible random intercept is included in a
#' standard lme4 style as the last term in the formula (e.g.
#' \code{+ (1 | subject_id)}).
#'
#' By default, the first term in the \code{formula} is treated as the variable
#' of interest, and the results are provided only for this variable. Any other
#' term of the formula can be selected as the variable of interest with
#' \code{var.of.interest}.
#'
#' The variable of interest must be either binary (a factor or character
#' variable with exactly two levels) or continuous. The other variables in the
#' model can, however, be factors with more than two levels.
#' }
#'
#' \subsection{Controlling for read depth (\code{read.depth})}{
#' As the observed presence/absence status of microbes may depend on the
#' sequencing (read) depth, DiPPER can control for it by adding log10 read
#' depth as a covariate. Because whether this is appropriate depends on
#' preprocessing that DiPPER cannot detect, \code{read.depth} must be chosen
#' explicitly for count data:
#' \itemize{
#'   \item Use \code{TRUE} for non-rarefied sequencing counts.
#'   \item Give a metadata column name when the read depths were computed
#'     manually, e.g. before heavy filtering of the count data.
#'   \item Use \code{FALSE} when read depth is already accounted for, e.g. by
#'     rarefying the data.
#' }
#' For other data types, read depth cannot be derived from the abundances
#' themselves, so it is controlled for only if a metadata column of
#' pre-computed read depths is supplied.
#' }
#'
#' \subsection{Symmetric or asymmetric prior (\code{symmetric})}{
#' By default, \code{dipper} uses the Asymmetric Laplace prior as the
#' hierarchical prior for the differential prevalence parameters of interest.
#' This means that the (a)symmetry of the Laplace prior is estimated from the
#' data. However, this choice may emphasize some technical biases in the data.
#' This can occur, for instance, when read depth correlates with the variable
#' of interest but one is not able to control for the read depth in the
#' analysis. Therefore, if such technical biases are possible or likely, the
#' user can choose to use a symmetric Laplace prior instead by setting
#' \code{symmetric = TRUE}.
#'
#' The two versions often give rather similar results, but the symmetric
#' version generally leads to slightly more robust and/or conservative results
#' while the asymmetric version is more sensitive to detect small effects.
#'
#' The parameters of the hierarchical prior are reported by
#' \code{\link{print.dipper_fit}}:
#' \itemize{
#'   \item \code{nu} indicates the asymmetry. It is constrained between 0 and
#'   1, with 0.5 meaning symmetry, values below 0.5 indicating positive
#'   skewness and values above 0.5 negative skewness.
#'   \item \code{tau} is the prior scale, that is, small values mean that the
#'   differential prevalence parameters are shrunk strongly towards zero/each
#'   other.
#' }
#'
#' If \code{nu} is clearly away from 0.5 (e.g. below 0.30 or above 0.70) and
#' \code{tau} is at the same time small (e.g. < 0.10), the asymmetry of the
#' prior may have a somewhat strong effect on the results. This may lead to a
#' large number of \sQuote{significant} findings in the same direction. If one
#' suspects that these findings may be due to some technical biases in the
#' data, it may be advisable to use the symmetric version of DiPPER.
#' }
#'
#' \subsection{MCMC sampling and prior settings}{
#' The posterior distribution computation for DiPPER is performed using the
#' Hamiltonian Monte Carlo algorithm via CmdStanR. The total number of
#' posterior samples can be controlled via arguments \code{niter} and
#' \code{chains}. Higher numbers lead to higher accuracy of the posterior
#' statistics but increase the computation time. Nevertheless, the default
#' values \code{niter = 2000} and \code{chains = 4} should be sufficient for
#' most cases. Higher values are generally recommended only if indicated by
#' the automatic MCMC diagnostics (e.g., too high R-hat values, or too low
#' effective sample sizes).
#'
#' Lastly, adjusting the default prior distribution settings
#' (\code{prior.alpha.sd}, \code{prior.tau.sd}, etc.) is generally not
#' recommended unless the user is highly experienced with Bayesian modeling
#' and the details of DiPPER.
#' }
#'
#' @references
#' Pelto, J., et al. (2026). DiPPER: A Bayesian approach to differential
#' prevalence analysis with applications in microbiome studies.
#' arXiv preprint. \url{https://arxiv.org/abs/2602.05938}
#'
#' @return A list object of class \code{dipper_fit} containing:
#' \describe{
#'   \item{draws}{A matrix of posterior draws, with one row per MCMC draw and
#'   one column per retained parameter. By default the draws of \code{beta}
#'   (the parameters of interest) and of the hierarchical prior parameters
#'   \code{tau} and \code{nu} are retained; see \code{keep.pars}.}
#'   \item{dipper_data}{A list containing the prepared data passed to Stan.}
#'   \item{symmetric}{Logical indicating if a symmetric Laplace prior for
#'   differential prevalence parameters was used.}
#'   \item{diagnostics}{A list of MCMC diagnostics. Always
#'   contains \code{num_divergent}, \code{num_max_treedepth},
#'   \code{time_total}, \code{n_chains} and \code{n_iter_sampling}. If
#'   \code{run.diagnostics = TRUE}, it additionally contains
#'   \code{max_rhat}, \code{min_ess_bulk}, \code{min_ess_tail} and
#'   \code{diagnostics_level}.}
#'   \item{stanfit}{The \code{CmdStanMCMC} object returned by CmdStanR.
#'   Present only if \code{keep.stanfit = TRUE}. Note that this is an R6
#'   object that depends on the cmdstanr namespace, so a fit containing it
#'   cannot be reloaded on a machine where cmdstanr is unavailable.}
#' }
#'
#' @export
#'
#' @examples
#' data("tse_hintikka")
#'
#' # Run DiPPER
#' # Note: niter = 400, chains = 2 and cores = 2 are used here for speed, so
#' # convergence warnings are expected. In real applications, use higher
#' # values (e.g. the default niter = 2000, chains = 4, and cores = 4).
#' if (instantiate::stan_cmdstan_exists()) {
#'     fit <- dipper(
#'         tse = tse_hintikka,
#'         assay.type = "counts",
#'         data.type = "counts",
#'         formula = ~ Fat + XOS,
#'         read.depth = TRUE,
#'         niter = 400,
#'         chains = 2,
#'         cores = 2
#'     )
#'
#'     print(fit)
#'
#'     res <- summary(fit)
#'     head(res)
#' }
dipper <- function(tse = NULL,
                   assay.type = NULL,
                   assay = NULL,
                   meta = NULL,
                   data.type = NULL,
                   formula,
                   var.of.interest = NULL,
                   read.depth = NULL,
                   symmetric = FALSE,
                   threshold = 0,
                   min.present = 5,
                   min.absent = min.present,
                   niter = 2000,
                   niter.warmup = floor(niter / 2),
                   chains = 4,
                   cores = 4,
                   adapt.delta = 0.8,
                   max.treedepth = 10,
                   run.diagnostics = TRUE,
                   diagnostics.level = c("basic", "full"),
                   keep.pars = c("beta", "tau", "nu"),
                   keep.stanfit = FALSE,
                   seed = 1,
                   print.progress = 200,
                   verbose = TRUE,
                   prior.alpha.sd = 4.0,
                   prior.tau.sd = 1.0,
                   prior.nu.sd = 0.05,
                   prior.cov.sd = 1.0,
                   prior.reads.mean = 2.0,
                   prior.reads.sd = 2.0,
                   prior.sigma.subj = 1.0,
                   ...) {

    # 0. Validate and match arguments ------------------------------------------

    # data.type and read.depth are validated by prep_dipper_data(), so that
    # both entry points give the same errors.
    diagnostics.level <- match.arg(diagnostics.level)


    # 1. Prepare the data ------------------------------------------------------
    prepared_data <- prep_dipper_data(
        tse = tse,
        assay.type = assay.type,
        assay = assay,
        meta = meta,
        data.type = data.type,
        formula = formula,
        var.of.interest = var.of.interest,
        read.depth = read.depth,
        threshold = threshold,
        min.present = min.present,
        min.absent = min.absent,
        verbose = verbose
    )

    # 2. Run the core model ----------------------------------------------------
    result <- run_dipper(
        prep.data = prepared_data,
        symmetric = symmetric,
        niter = niter,
        niter.warmup = niter.warmup,
        chains = chains,
        cores = cores,
        adapt.delta = adapt.delta,
        max.treedepth = max.treedepth,
        run.diagnostics = run.diagnostics,
        diagnostics.level = diagnostics.level,
        keep.pars = keep.pars,
        keep.stanfit = keep.stanfit,
        seed = seed,
        print.progress = print.progress,
        verbose = verbose,
        prior.alpha.sd = prior.alpha.sd,
        prior.tau.sd = prior.tau.sd,
        prior.nu.sd = prior.nu.sd,
        prior.cov.sd = prior.cov.sd,
        prior.reads.mean = prior.reads.mean,
        prior.reads.sd = prior.reads.sd,
        prior.sigma.subj = prior.sigma.subj,
        ...
    )

    return(result)
}
