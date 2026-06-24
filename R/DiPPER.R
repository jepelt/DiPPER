#' Run DiPPER
#'
#' This is a wrapper function that sequentially runs \code{prep_dipper_data}
#' and \code{run_dipper}.
#'
#' @inheritParams prep_dipper_data
#' @inheritParams run_dipper
#' @param ... Additional arguments passed directly to \code{run_dipper}
#'   (e.g., cmdstanr specific arguments).
#'
#' @return A list object of class "dipper_fit".
#' @export
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
                   warmup = floor(niter / 2),
                   chains = 4,
                   cores = 4,
                   adapt.delta = 0.8,
                   max.treedepth = 10,
                   run.diagnostics = TRUE,
                   diagnostics.level = c("basic", "full"),
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
        warmup = warmup,
        chains = chains,
        cores = cores,
        adapt.delta = adapt.delta,
        max.treedepth = max.treedepth,
        run.diagnostics = run.diagnostics,
        diagnostics.level = diagnostics.level,
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
