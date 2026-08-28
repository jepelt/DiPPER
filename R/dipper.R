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
#' **variable of interest**
#' (e.g. treatment group, disease status, etc.) and potentially other
#' covariates (e.g. age, sex, etc.) that may confound the associations between
#' microbe prevalences and the variable of interest. DiPPER can also be applied
#' to longitudinal or repeated measures data.
#'
#' Technically, DiPPER models the presence/absence of taxonomic features (e.g.
#' genera or species) using multiple logistic regression models. The models are
#' connected by a common (Asymmetric) Laplace prior for the variable of
#' interest. The hierarchical structure guarantees robust (and always finite)
#' differential prevalence estimates and uncertainty intervals that are also
#' effectively multiplicity-adjusted.
#'
#' The \code{dipper} function takes either a \code{(Tree)SummarizedExperiment}
#' object or an abundance matrix and a metadata \code{data.frame} as input. The
#' abundance matrix can be either counts (default), relative abundances, or
#' presence/absence data. This is indicated by the \code{data.type} argument.
#'
#' Before model fitting, the input data undergoes automated preprocessing.
#' Continuous abundance data (counts or relative abundances) are converted to
#' presence/absence format based on the specified \code{threshold}. Furthermore,
#' prevalence filtering is applied: taxa that are not present in at least
#' \code{min.present} samples, or absent in at least \code{min.absent} samples,
#' are excluded from the analysis.
#'
#' The \code{dipper} function uses the standard \code{formula} argument to
#' specify the model formula. A possible random intercept is included in a
#' standard lme4 style as the last term in the formula (e.g.
#' \code{+ (1 | subject_id)}).
#'
#' By default, the first term in the \code{formula} is treated as the variable
#' of interest. It is important to note that the hierarchical (Asymmetric)
#' Laplace prior, and thus the automatic multiplicity adjustment, is applied
#' exclusively to this variable.
#'
#' As the observed presence/absence status of microbes may depend on the
#' sequencing (read) depth, DiPPER by default controls for this. This is done
#' automatically if a counts matrix is provided and \code{read.depth = TRUE}.
#' If only relative abundances (proportions) are available, one must set
#' \code{read.depth = FALSE} (or provide a separate sequencing depth variable).
#'
#' By default, \code{dipper} uses the Asymmetric Laplace prior as the
#' hierarchical prior. This encodes an empirical observation that within a given
#' microbiome study, most non-zero prevalence differences tend to share the same
#' direction. However, this choice may sometimes lead to an undesirably strong
#' bias in the estimates. Therefore, the user can choose to use a symmetric
#' Laplace prior instead by setting \code{symmetric = TRUE}. This option gives
#' generally more robust and/or conservative results that may be more in line
#' with the results of standard (frequentist) approaches.
#'
#' The posterior distribution computation for DiPPER is performed using the
#' Hamiltonian Monte Carlo algorithm via CmdStanR. The total number of posterior
#' samples can be controlled via arguments \code{niter} and \code{chains}.
#' Higher numbers lead to higher accuracy of the posterior statistics but
#' increase the computation time. Nevertheless, the default values
#' \code{niter = 2000} and \code{chains = 4} should be sufficient for most
#' cases. Higher values are generally recommended only if indicated by the
#' automatic MCMC diagnostics (e.g., too high R-hat values, or too low
#' effective sample sizes).
#'
#' Lastly, adjusting the default prior distribution settings
#' (\code{prior.alpha.sd}, \code{prior.tau.sd}, etc.) is generally not
#' recommended unless the user is highly experienced with Bayesian modeling
#' and the details of DiPPER.
#'
#' @references
#' Pelto, J., et al. (2026). DiPPER: A Bayesian approach to differential
#' prevalence analysis with applications in microbiome studies.
#' arXiv preprint. \url{https://arxiv.org/abs/2602.05938}
#'
#' @return A list object of class \code{dipper_fit} containing:
#' \describe{
#'   \item{draws}{A matrix of posterior draws, with one row per MCMC draw and
#'   one column per retained parameter. By default only the draws of
#'   \code{beta} (the parameters of interest) are retained; see
#'   \code{keep.pars}.}
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
#' # Note: niter = 400, chains = 1 and cores = 1 are used here for speed.
#' # In real applications, use higher values (e.g. the default niter = 2000,
#' # chains = 4, and cores = 4).
#' if (instantiate::stan_cmdstan_exists()) {
#'     fit <- dipper(
#'         tse = tse_hintikka,
#'         formula = ~ Fat + XOS,
#'         assay.type = "counts",
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
                   formula,
                   assay = NULL,
                   meta = NULL,
                   assay.type = NULL,
                   data.type = c("counts", "relabundance", "pa"),
                   var.of.interest = NULL,
                   read.depth = TRUE,
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
                   keep.pars = "beta",
                   keep.stanfit = FALSE,
                   seed = 1,
                   print.progress = 200,
                   prior.alpha.sd = 4.0,
                   prior.tau.sd = 1.0,
                   prior.nu.sd = 0.05,
                   prior.cov.sd = 1.0,
                   prior.reads.mean = 2.0,
                   prior.reads.sd = 2.0,
                   prior.sigma.subj = 1.0,
                   ...) {

    # 0. Validate and match arguments ------------------------------------------
    diagnostics.level <- match.arg(diagnostics.level)
    data.type <- match.arg(data.type)


    # 1. Prepare the data ------------------------------------------------------
    prepared_data <- prep_dipper_data(
        tse = tse,
        formula = formula,
        assay = assay,
        meta = meta,
        assay.type = assay.type,
        data.type = data.type,
        var.of.interest = var.of.interest,
        read.depth = read.depth,
        threshold = threshold,
        min.present = min.present,
        min.absent = min.absent
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
