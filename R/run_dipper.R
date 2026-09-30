#' Validate the arguments of run_dipper()
#'
#' @return \code{TRUE} invisibly. Called for its errors.
#' @keywords internal
#' @noRd
.dipper_validate_run_args <- function(prep.data, keep.pars, keep.stanfit,
                                      niter, niter.warmup) {
    req_names <- c("y", "X", "N", "K", "P", "design_matrix_cols")
    if (!all(req_names %in% names(prep.data))) {
        stop("Invalid prep.data. Generate using prep_dipper_data().",
             call. = FALSE)
    }

    if (!is.null(keep.pars) && !is.character(keep.pars)) {
        stop("'keep.pars' must be a character vector or NULL.",
             call. = FALSE)
    }

    if (!is.logical(keep.stanfit) || length(keep.stanfit) != 1) {
        stop("'keep.stanfit' must be TRUE or FALSE.", call. = FALSE)
    }

    if (niter - niter.warmup <= 0) {
        stop("'niter' must be strictly greater than 'niter.warmup'.",
             call. = FALSE)
    }

    invisible(TRUE)
}


#' Prior means and standard deviations for the covariates
#'
#' The variable of interest (first design matrix column) is excluded. The read
#' depth covariate, if present, gets its own prior.
#'
#' @return A list with numeric vectors \code{mean} and \code{sd}, each of
#'   length \code{P - 1}.
#' @keywords internal
#' @noRd
.dipper_cov_priors <- function(prep.data, prior.cov.sd,
                               prior.reads.mean, prior.reads.sd) {
    P <- prep.data$P

    if (P <= 1) {
        return(list(mean = numeric(0), sd = numeric(0)))
    }

    cov_names <- prep.data$design_matrix_cols[-1]
    p_mean <- rep(0, P - 1)
    p_sd <- rep(prior.cov.sd, P - 1)

    if (!is.null(prep.data$read.depth.var)) {
        read_idx <- which(cov_names == prep.data$read.depth.var)
        if (length(read_idx) > 0) {
            p_mean[read_idx] <- prior.reads.mean
            p_sd[read_idx] <- prior.reads.sd
        }
    }

    list(mean = p_mean, sd = p_sd)
}


#' Build the Stan data list and select the Stan model
#'
#' @return A list with elements \code{stan_data} (list passed to Stan) and
#'   \code{model_name} (name of the Stan model without the .stan extension).
#' @keywords internal
#' @noRd
.dipper_stan_setup <- function(prep.data, symmetric,
                               prior.alpha.sd, prior.tau.sd, prior.nu.sd,
                               prior.cov.sd, prior.reads.mean, prior.reads.sd,
                               prior.sigma.subj) {
    cov_priors <- .dipper_cov_priors(
        prep.data = prep.data,
        prior.cov.sd = prior.cov.sd,
        prior.reads.mean = prior.reads.mean,
        prior.reads.sd = prior.reads.sd
    )

    stan_data <- list(
        N = prep.data$N,
        K = prep.data$K,
        P = prep.data$P,
        y = prep.data$y,
        X = prep.data$X,
        prior_alpha_mean = 0.0,
        prior_alpha_sd = prior.alpha.sd,
        prior_tau_sd = prior.tau.sd,
        prior_cov_mean = as.array(cov_priors$mean),
        prior_cov_sd = as.array(cov_priors$sd)
    )

    is_longitudinal <- isTRUE(prep.data$is_longitudinal)

    if (is_longitudinal) {
        stan_data$S <- prep.data$S
        stan_data$subj <- prep.data$subj
        stan_data$prior_sigma_subj <- prior.sigma.subj
    }

    # The asymmetry parameter nu exists only in the asymmetric models
    if (!symmetric) {
        stan_data$prior_nu_sd <- prior.nu.sd
    }

    model_name <- paste0(
        "dipper_dp_",
        if (is_longitudinal) "long_" else "",
        if (symmetric) "sym" else "asym"
    )

    list(stan_data = stan_data, model_name = model_name)
}


#' Translate print.progress and verbose into CmdStanR settings
#'
#' @return A list with elements \code{refresh} (integer) and
#'   \code{show_messages} (logical).
#' @keywords internal
#' @noRd
.dipper_progress_settings <- function(print.progress, verbose) {
    # verbose = FALSE always silences the sampler's progress output as well.
    if (!verbose) {
        print.progress <- FALSE
    }

    if (is.logical(print.progress)) {
        refresh <- if (print.progress) 200L else 0L
    } else if (is.numeric(print.progress) && print.progress > 0) {
        refresh <- as.integer(print.progress)
    } else {
        refresh <- 0L
    }

    list(refresh = refresh, show_messages = refresh > 0)
}


#' Check that a DiPPER Stan model has been compiled
#'
#' The Stan models are compiled when DiPPER is installed. If CmdStan was
#' installed only afterwards, the executables are missing and sampling would
#' fail with an uninformative cmdstanr error.
#'
#' @param mod A CmdStanModel object.
#' @return \code{TRUE} invisibly. Called for its error.
#' @keywords internal
#' @noRd
.dipper_check_compiled <- function(mod) {
    exe <- tryCatch(mod$exe_file(), error = function(e) character(0))

    if (length(exe) == 0 || !nzchar(exe[1]) || !file.exists(exe[1])) {
        stop(
            "The DiPPER Stan models have not been compiled. This happens if ",
            "DiPPER was installed before CmdStan. Reinstall DiPPER, then ",
            "restart R.",
            call. = FALSE
        )
    }

    invisible(TRUE)
}


#' Parameters checked by the convergence diagnostics
#'
#' @return Character vector of Stan parameter names, or \code{NULL} to check
#'   all parameters, including latent and auxiliary ones.
#' @keywords internal
#' @noRd
.dipper_diagnostic_vars <- function(diagnostics.level, symmetric, P) {
    if (diagnostics.level == "full") {
        return(NULL)
    }

    target_vars <- c("alpha", "beta", "tau")

    if (!symmetric) {
        target_vars <- c(target_vars, "nu")
    }
    if (P > 1) {
        target_vars <- c(target_vars, "beta_cov")
    }

    target_vars
}


#' Maximum returning NA when x contains only NAs
#'
#' @keywords internal
#' @noRd
.dipper_max <- function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}


#' Minimum returning NA when x contains only NAs
#'
#' @keywords internal
#' @noRd
.dipper_min <- function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}


#' Convergence warning messages from MCMC diagnostics
#'
#' @param diagnostics List with elements \code{num_divergent},
#'   \code{max_rhat}, \code{min_ess_bulk} and \code{min_ess_tail}.
#' @param niter Total number of iterations per chain, used for the
#'   recommendation.
#' @return Character vector of messages (empty if no issues were found).
#' @keywords internal
#' @noRd
.dipper_convergence_warnings <- function(diagnostics, niter) {
    warn_msg <- character(0)
    rec_iter <- niter * 2

    divs <- diagnostics$num_divergent
    max_rhat <- diagnostics$max_rhat
    min_ess_bulk <- diagnostics$min_ess_bulk
    min_ess_tail <- diagnostics$min_ess_tail

    if (divs > 0) {
        warn_msg <- c(warn_msg, sprintf(
            "%d divergent transitions. Try adapt.delta = 0.95.", divs
        ))
    }
    if (is.finite(max_rhat) && max_rhat >= 1.01) {
        warn_msg <- c(warn_msg, sprintf(
            "Max R-hat is %.3f. Try niter = %d.", max_rhat, rec_iter
        ))
    }
    if (is.finite(min_ess_bulk) && min_ess_bulk < 400) {
        warn_msg <- c(warn_msg, sprintf(
            "Min bulk ESS is %.1f. Try niter = %d.", min_ess_bulk, rec_iter
        ))
    }
    if (is.finite(min_ess_tail) && min_ess_tail < 400) {
        warn_msg <- c(warn_msg, sprintf(
            "Min tail ESS is %.1f. Try niter = %d.", min_ess_tail, rec_iter
        ))
    }

    warn_msg
}


#' Fit the DiPPER model using CmdStan
#'
#' @param prep.data The list returned by \code{prep_dipper_data}.
#' @param symmetric Logical. If TRUE, a symmetric Laplace prior (nu = 0.5) is
#'   used for the differential prevalence parameters of interest. If FALSE
#'   (default), the asymmetry parameter nu is estimated from the data. See the
#'   Details section of \code{\link{dipper}} for guidance on choosing between
#'   the two.
#' @param niter Total number of MCMC iterations per chain (including warmup).
#'   Default is 2000.
#' @param niter.warmup Number of warmup MCMC iterations per chain.
#'   Defaults to floor(niter / 2).
#' @param chains Number of MCMC chains. Default is 4.
#' @param cores Number of CPU cores to use. Default is 4.
#' @param seed Random seed for reproducibility. Default is 1.
#' @param adapt.delta Target average acceptance probability in MCMC sampling.
#'   Default is 0.8.
#' @param max.treedepth Maximum depth of the trees in MCMC sampling. Default is
#'   10.
#' @param run.diagnostics Logical. Whether to run convergence diagnostics
#'   after MCMC sampling. Default is TRUE.
#' @param diagnostics.level Character. Level of MCMC diagnostics to run.
#'   "basic" (default) checks main parameters (alpha, beta, covariates,
#'   hyperparameters).
#'   "full" checks all parameters, including latent and auxiliary variables.
#' @param keep.pars Character vector of Stan parameter names whose posterior
#'   draws are retained in the returned object, or NULL to retain all of them.
#'   The default, \code{c("beta", "tau", "nu")}, keeps the differential
#'   prevalence parameters of the variable of interest together with the two
#'   hyperparameters of the hierarchical Laplace prior. \code{nu} is dropped
#'   automatically when \code{symmetric = TRUE}.
#' @param keep.stanfit Logical. Whether to include the underlying CmdStanR
#'   fit R6 object in the returned object (e.g. for manual model fit checking).
#'   Default is FALSE.
#' @param print.progress How often to print MCMC progress. Set to 0 or
#'   FALSE to disable. Default is 200. Ignored if verbose = FALSE.
#' @param verbose Logical. Whether to print any progress messages issued by
#'   DiPPER functions. Default is TRUE. Setting verbose = FALSE overrides
#'   \code{print.progress}.
#' @param prior.alpha.sd Prior standard deviation for alpha (fixed intercept
#'   for centered variables). Default is 4.0.
#' @param prior.tau.sd Prior standard deviation for tau (global scale of the
#'   Symmetric/Asymmetric Laplace prior). Default is 1.0.
#' @param prior.nu.sd Prior standard deviation for nu (asymmetry parameter).
#'   Ignored if symmetric = TRUE. Default is 0.05.
#' @param prior.cov.sd Default prior standard deviation for covariates.
#'   Default is 1.0.
#' @param prior.reads.mean Prior mean for the read depth covariate.
#'   Default is 2.0.
#' @param prior.reads.sd Prior standard deviation for the read depth covariate.
#'   Default is 2.0.
#' @param prior.sigma.subj Prior standard deviation for the half-normal prior
#'   distribution of the taxon-specific random intercept standard deviations.
#'   Default is 1.0.
#' @param ... Additional arguments passed to the \code{$sample()} method of the
#'   CmdStanR model object.
#' @return A list of class \code{dipper_fit} with elements \code{draws} (a
#'   matrix of posterior draws, with one column per retained parameter. See
#'   \code{keep.pars}), \code{dipper_data}, \code{symmetric} and
#'   \code{diagnostics}. If \code{keep.stanfit = TRUE}, the CmdStanR fit R6
#'   object is included as \code{stanfit}.
#' @keywords internal
run_dipper <- function(prep.data,
                       symmetric = FALSE,
                       niter = 2000,
                       niter.warmup = floor(niter / 2),
                       chains = 4,
                       cores = 4,
                       seed = 1,
                       adapt.delta = 0.8,
                       max.treedepth = 10,
                       run.diagnostics = TRUE,
                       diagnostics.level = c("basic", "full"),
                       keep.pars = c("beta", "tau", "nu"),
                       keep.stanfit = FALSE,
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

    diagnostics.level <- match.arg(diagnostics.level)


    # 1. Input validation ------------------------------------------------------
    .dipper_validate_run_args(
        prep.data = prep.data,
        keep.pars = keep.pars,
        keep.stanfit = keep.stanfit,
        niter = niter,
        niter.warmup = niter.warmup
    )

    if (!instantiate::stan_cmdstan_exists()) {
        stop(
            "CmdStan not found. Install it with cmdstanr::install_cmdstan(), ",
            "then reinstall DiPPER so that the Stan models get compiled.",
            call. = FALSE
        )
    }

    # The asymmetry parameter nu only exists in the asymmetric models, so
    # asking the sampler for it would fail
    if (symmetric && !is.null(keep.pars) && "nu" %in% keep.pars) {
        keep.pars <- setdiff(keep.pars, "nu")
        if (length(keep.pars) == 0) {
            keep.pars <- "beta"
        }
    }


    # 2. Stan data and model selection -----------------------------------------
    setup <- .dipper_stan_setup(
        prep.data = prep.data,
        symmetric = symmetric,
        prior.alpha.sd = prior.alpha.sd,
        prior.tau.sd = prior.tau.sd,
        prior.nu.sd = prior.nu.sd,
        prior.cov.sd = prior.cov.sd,
        prior.reads.mean = prior.reads.mean,
        prior.reads.sd = prior.reads.sd,
        prior.sigma.subj = prior.sigma.subj
    )


    # 3. Configure MCMC printing -----------------------------------------------
    progress <- .dipper_progress_settings(print.progress, verbose)


    # 4. Load the model and run sampling ---------------------------------------
    actual_sampling <- niter - niter.warmup

    if (verbose) {
        message("Preparing Stan model...")
    }

    mod <- instantiate::stan_package_model(
        name = setup$model_name,
        package = "DiPPER"
    )
    .dipper_check_compiled(mod)

    if (verbose) {
        message(sprintf(
            "Starting sampling with %d chains on %d cores...", chains, cores
        ))
    }

    if (progress$refresh == 0 && verbose) {
        message("Progress printing is disabled. Please wait...")
    }

    fit <- mod$sample(
        data = setup$stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = cores,
        iter_warmup = niter.warmup,
        iter_sampling = actual_sampling,
        adapt_delta = adapt.delta,
        max_treedepth = max.treedepth,
        refresh = progress$refresh,
        init = 0.1,
        show_messages = progress$show_messages,
        show_exceptions = FALSE,
        ...
    )


    # 5. Diagnostics -----------------------------------------------------------
    diag_sum <- fit$diagnostic_summary(quiet = TRUE)

    diagnostics <- list(
        num_divergent = sum(diag_sum$num_divergent),
        num_max_treedepth = sum(diag_sum$num_max_treedepth),
        time_total = fit$time()$total,
        n_chains = chains,
        n_iter_sampling = actual_sampling
    )

    if (run.diagnostics) {
        if (verbose) {
            message("Sampling completed. Checking diagnostics...")
        }

        target_vars <- .dipper_diagnostic_vars(
            diagnostics.level = diagnostics.level,
            symmetric = symmetric,
            P = prep.data$P
        )

        # posterior::ess_* warns when the ESS estimate is capped at the number
        # of draws. Only the minimum ESS is used below, so the cap is
        # irrelevant here and the warning would fire on every run.
        summ <- suppressWarnings(
            fit$summary(target_vars, "rhat", "ess_bulk", "ess_tail")
        )

        diagnostics$max_rhat <- .dipper_max(summ$rhat)
        diagnostics$min_ess_bulk <- .dipper_min(summ$ess_bulk)
        diagnostics$min_ess_tail <- .dipper_min(summ$ess_tail)
        diagnostics$diagnostics_level <- diagnostics.level

        warn_msg <- .dipper_convergence_warnings(diagnostics, niter)

        if (length(warn_msg) > 0) {
            warning(
                "Convergence issues detected:\n- ",
                paste(warn_msg, collapse = "\n- "),
                call. = FALSE
            )
        } else if (verbose) {
            message("All MCMC diagnostics are within acceptable limits.")
        }
    }


    # 6. Extract the posterior draws -------------------------------------------
    if (!is.null(keep.pars) && length(keep.pars) > 0) {
        draws <- fit$draws(variables = keep.pars, format = "matrix")
        if (keep.stanfit) {
            fit$.__enclos_env__$private$draws_ <-
                fit$draws(variables = keep.pars)
        }
    } else {
        draws <- fit$draws(format = "matrix")
    }

    draws <- as.matrix(draws)

    # Drop the posterior draws_matrix class and its attributes, to avoid
    # dependency on the posterior package
    attributes(draws) <- attributes(draws)[c("dim", "dimnames")]


    # 7. Return output object -------------------------------------------------
    out <- list(
        draws = draws,
        dipper_data = prep.data,
        symmetric = symmetric,
        diagnostics = diagnostics
    )

    if (keep.stanfit) {
        out$stanfit <- fit
    }

    structure(out, class = "dipper_fit")
}
