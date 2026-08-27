#' Run the DiPPER model using CmdStan
#'
#' @param prep.data The list returned by \code{prep_dipper_data}.
#' @param symmetric Logical. If TRUE, a symmetric Laplace prior (nu = 0.5)
#'   is used for the differential prevalence parameters of interest. If FALSE
#'   (default), the prior is allowed to be asymmetric.
#' @param niter Total number of MCMC iterations per chain (including warmup).
#'   Default is 2000.
#' @param niter.warmup Number of warmup MCMC iterations per chain.
#'   Defaults to floor(niter / 2).
#' @param chains Number of MCMC chains. Default is 4.
#' @param cores Number of CPU cores to use. Default is 4.
#' @param seed Random seed for reproducibility. Default is 1.
#' @param adapt.delta Target average acceptance probability. Default is 0.8.
#' @param max.treedepth Maximum depth of the trees. Default is 10.
#' @param run.diagnostics Logical. Whether to run convergence diagnostics
#'   after sampling. Default is TRUE.
#' @param diagnostics.level Character. Level of MCMC diagnostics to run.
#'   "basic" (default) checks main parameters (alpha, beta, covariates,
#'   hyperparameters).
#'   "full" checks all parameters, including latent and auxiliary variables.
#' @param keep.pars Character vector of Stan parameter names whose posterior
#'   draws are retained in the returned object, or NULL to retain all of them.
#'   The default, "beta", keeps only the differential prevalence parameters of
#'   the variable of interest, which is what \code{summary} and \code{plot}
#'   need. Note that this argument does not affect convergence diagnostics.
#' @param keep.stanfit Logical. Whether to include the underlying CmdStanR
#'   fit R6 object in the returned object (e.g. for manual model fit checking).
#'   Default is FALSE.
#' @param print.progress How often to print MCMC progress. Set to 0 or
#'   FALSE to disable. Default is 200.
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
#'   matrix of posterior draws, with one column per retained parameter),
#'   \code{dipper_data}, \code{symmetric} and \code{diagnostics}. If
#'   \code{keep.stanfit = TRUE}, the CmdStanR fit R6 object is included as
#'   \code{stanfit}.
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
                       keep.pars = "beta",
                       keep.stanfit = FALSE,
                       print.progress = 200,
                       prior.alpha.sd = 4.0,
                       prior.tau.sd = 1.0,
                       prior.nu.sd = 0.05,
                       prior.cov.sd = 1.0,
                       prior.reads.mean = 2.0,
                       prior.reads.sd = 2.0,
                       prior.sigma.subj = 1.0,
                       ...) {

    diagnostics.level <- match.arg(diagnostics.level)


    # 1. Dependency and input validation ---------------------------------------
    if (!instantiate::stan_cmdstan_exists()) {
        stop(
            "CmdStan not found. Install it with cmdstanr::install_cmdstan(), ",
            "then reinstall DiPPER so that the Stan models get compiled.",
            call. = FALSE
        )
    }

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


    # 2. Extract covariates and setup priors -----------------------------------
    P <- prep.data$P
    col_names <- prep.data$design_matrix_cols

    if (P > 1) {
        cov_names <- col_names[-1]
        p_mean <- rep(0, P - 1)
        p_sd <- rep(prior.cov.sd, P - 1)

        if (!is.null(prep.data$read.depth.var)) {
            read_idx <- grep(prep.data$read.depth.var, cov_names)
            if (length(read_idx) > 0) {
                p_mean[read_idx] <- prior.reads.mean
                p_sd[read_idx] <- prior.reads.sd
            }
        }
    } else {
        p_mean <- numeric(0)
        p_sd <- numeric(0)
    }


    # 3. Build base Stan data list ---------------------------------------------
    stan_data <- list(
        N = prep.data$N,
        K = prep.data$K,
        P = P,
        y = prep.data$y,
        X = prep.data$X,
        prior_alpha_mean = 0.0,
        prior_alpha_sd = prior.alpha.sd,
        prior_tau_sd = prior.tau.sd,
        prior_cov_mean = as.array(p_mean),
        prior_cov_sd = as.array(p_sd)
    )


    # 4. Select Stan code based on longitudinality and (a)symmetry -------------
    is_longitudinal <- isTRUE(prep.data$is_longitudinal)

    if (is_longitudinal) {
        stan_data$S <- prep.data$S
        stan_data$subj <- prep.data$subj
        stan_data$prior_sigma_subj <- prior.sigma.subj

        if (!symmetric) {
            stan_data$prior_nu_sd <- prior.nu.sd
            file_name <- "dipper_dp_long_asym.stan"
        } else {
            file_name <- "dipper_dp_long_sym.stan"
        }
    } else {
        if (!symmetric) {
            stan_data$prior_nu_sd <- prior.nu.sd
            file_name <- "dipper_dp_asym.stan"
        } else {
            file_name <- "dipper_dp_sym.stan"
        }
    }


    # 5. Configure MCMC printing -----------------------------------------------
    if (is.logical(print.progress)) {
        if (print.progress) {
            refresh_val <- 200L
            show_msgs <- TRUE
        } else {
            refresh_val <- 0L
            show_msgs <- FALSE
        }
    } else if (is.numeric(print.progress) && print.progress > 0) {
        refresh_val <- as.integer(print.progress)
        show_msgs <- TRUE
    } else {
        refresh_val <- 0L
        show_msgs <- FALSE
    }


    # 6. Compile and run sampling ----------------------------------------------

    # Calculate actual sampling iterations for the sampler
    actual_sampling <- niter - niter.warmup
    if (actual_sampling <= 0) {
        stop("'niter' must be strictly greater than 'niter.warmup'.",
             call. = FALSE)
    }

    message("Preparing Stan model...")

    model_name <- sub("\\.stan$", "", file_name)

    mod <- instantiate::stan_package_model(
        name = model_name,
        package = "DiPPER"
    )

    message(sprintf(
        "Starting sampling with %d chains on %d cores...", chains, cores
    ))

    if (refresh_val == 0) {
        message("Progress printing is disabled. Please wait...")
    }

    fit <- mod$sample(
        data = stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = cores,
        iter_warmup = niter.warmup,
        iter_sampling = actual_sampling,
        adapt_delta = adapt.delta,
        max_treedepth = max.treedepth,
        refresh = refresh_val,
        init = 0.1,
        show_messages = show_msgs,
        show_exceptions = FALSE,
        ...
    )


    # 7. Diagnostics -----------------------------------------------------------
    diag_sum <- fit$diagnostic_summary(quiet = TRUE)

    diagnostics <- list(
        num_divergent = sum(diag_sum$num_divergent),
        num_max_treedepth = sum(diag_sum$num_max_treedepth),
        time_total = fit$time()$total,
        n_chains = chains,
        n_iter_sampling = actual_sampling
    )

    if (run.diagnostics) {
        message("Sampling completed. Checking diagnostics...")

        divs <- diagnostics$num_divergent

        if (diagnostics.level == "basic") {
            target_vars <- c("alpha", "beta", "tau")

            if (!symmetric) {
                target_vars <- c(target_vars, "nu")
            }
            if (P > 1) {
                target_vars <- c(target_vars, "beta_cov")
            }
        } else {
            target_vars <- NULL
        }

        summ <- suppressWarnings(
            fit$summary(target_vars, "rhat", "ess_bulk", "ess_tail")
        )

        max_rhat <- suppressWarnings(max(summ$rhat, na.rm = TRUE))
        min_ess_bulk <- suppressWarnings(min(summ$ess_bulk, na.rm = TRUE))
        min_ess_tail <- suppressWarnings(min(summ$ess_tail, na.rm = TRUE))

        diagnostics$max_rhat <- max_rhat
        diagnostics$min_ess_bulk <- min_ess_bulk
        diagnostics$min_ess_tail <- min_ess_tail
        diagnostics$diagnostics_level <- diagnostics.level

        warn_msg <- c()
        rec_iter <- niter * 2

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
                "Min bulk ESS is %.1f. Try niter = %d.",
                min_ess_bulk, rec_iter
            ))
        }
        if (is.finite(min_ess_tail) && min_ess_tail < 400) {
            warn_msg <- c(warn_msg, sprintf(
                "Min tail ESS is %.1f. Try niter = %d.",
                min_ess_tail, rec_iter
            ))
        }

        if (length(warn_msg) > 0) {
            warning(
                "Convergence issues detected:\n- ",
                paste(warn_msg, collapse = "\n- "),
                call. = FALSE
            )
        } else {
            message("All MCMC diagnostics are within acceptable limits.")
        }
    }


    # 8. Extract the posterior draws -------------------------------------------
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


    # 9. Return output object -------------------------------------------------
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
