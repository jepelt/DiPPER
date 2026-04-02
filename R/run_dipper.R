#' Run the DiPPER model using cmdstanr
#'
#' @param prep_data The list returned by \code{prep_dipper_data}.
#' @param iter_sampling Number of post-warmup MCMC iterations per chain.
#' @param iter_warmup Number of warmup MCMC iterations per chain. Defaults
#'   to iter_sampling.
#' @param run_diagnostics Logical. Whether to run convergence diagnostics.
#'   Default is TRUE.
#' @param chains Number of MCMC chains.
#' @param cores Number of CPU cores to use.
#' @param seed Random seed for reproducibility.
#' @param adapt_delta Target average acceptance probability (default 0.8).
#' @param max_treedepth Maximum depth of the trees (default 10).
#' @param prior_nu_sd Prior standard deviation for nu.
#' @param prior_tau_sd Prior standard deviation for tau.
#' @param prior_alpha_sd Prior standard deviation for alpha (intercept).
#' @param prior_cov_sd Default prior standard dev for covariates.
#' @param prior_reads_mean Prior mean for the read depth covariate.
#' @param prior_reads_sd Prior standard dev for the read depth covariate.
#' @param print_progress How often to print MCMC progress. Set to 0 or
#'   FALSE to disable. Default is 200.
#'
#' @export
run_dipper <- function(prep_data,
                       iter_sampling = 1000,
                       iter_warmup = iter_sampling,
                       run_diagnostics = TRUE,
                       chains = 4,
                       cores = 4,
                       seed = 1,
                       adapt_delta = 0.8,
                       max_treedepth = 10,
                       prior_nu_sd = 0.05,
                       prior_tau_sd = 1.0,
                       prior_alpha_sd = 4.0,
                       prior_cov_sd = 1.0,
                       prior_reads_mean = 2.0,
                       prior_reads_sd = 2.0,
                       print_progress = 200) {

    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
        stop("cmdstanr required. Install via cmdstanr::install_cmdstanr()")
    }

    req_names <- c("y", "X", "N", "K", "P", "design_matrix_cols")
    if (!all(req_names %in% names(prep_data))) {
        stop("Invalid prep_data. Generate using prep_dipper_data().")
    }

    P <- prep_data$P
    col_names <- prep_data$design_matrix_cols

    if (P > 1) {
        cov_names <- col_names[-1]
        p_mean <- rep(0, P - 1)
        p_sd <- rep(prior_cov_sd, P - 1)

        if (!is.null(prep_data$read_depth_var)) {
            read_idx <- grep(prep_data$read_depth_var, cov_names)
            if (length(read_idx) > 0) {
                p_mean[read_idx] <- prior_reads_mean
                p_sd[read_idx] <- prior_reads_sd
            }
        }
    } else {
        p_mean <- numeric(0)
        p_sd <- numeric(0)
    }

    stan_data <- list(
        N = prep_data$N,
        K = prep_data$K,
        P = P,
        y = prep_data$y,
        X = prep_data$X,
        prior_alpha_mean = 0.0,
        prior_alpha_sd = prior_alpha_sd,
        prior_nu_sd = prior_nu_sd,
        prior_tau_sd = prior_tau_sd,
        prior_cov_mean = as.array(p_mean),
        prior_cov_sd = as.array(p_sd)
    )

    stan_file <- system.file(
        "stan", "dipper_dp_asym.stan",
        package = "DiPPER"
    )

    if (stan_file == "") {
        stop("Stan file not found in 'inst/stan/'.")
    }

    if (is.logical(print_progress) && !print_progress) {
        refresh_val <- 0
        show_msgs <- FALSE
    } else if (is.numeric(print_progress) && print_progress > 0) {
        refresh_val <- as.integer(print_progress)
        show_msgs <- TRUE
    } else {
        refresh_val <- 0
        show_msgs <- FALSE
    }

    message("Compiling or loading cached Stan model...")
    message("Note: The first run may take a while as the model is being compiled.")
    mod <- suppressWarnings(suppressMessages(
        cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
    ))

    message(sprintf(
        "Starting MCMC sampling with %d chains on %d cores...",
        chains, cores
    ))

    if (refresh_val == 0) {
        message("Progress printing is disabled. Please wait...")
    }

    fit <- mod$sample(
        data = stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = cores,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = refresh_val,
        init = 0.1,
        show_messages = show_msgs,
        show_exceptions = FALSE
    )

    if (run_diagnostics) {
        message("Sampling completed. Checking diagnostics...")

        diag_sum <- fit$diagnostic_summary()
        divs <- sum(diag_sum$num_divergent)

        summ <- fit$summary(NULL, "rhat", "ess_bulk", "ess_tail")
        max_rhat <- max(summ$rhat, na.rm = TRUE)
        min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
        min_ess_tail <- min(summ$ess_tail, na.rm = TRUE)

        warn_msg <- c()
        rec_iter <- iter_sampling * 2

        if (divs > 0) {
            warn_msg <- c(warn_msg, sprintf(
                "%d divergent transitions. Increase adapt_delta (e.g., 0.95).",
                divs
            ))
        }
        if (max_rhat >= 1.01) {
            warn_msg <- c(warn_msg, sprintf(
                "Max R-hat is %.3f. Try iter_sampling = %d.",
                max_rhat, rec_iter
            ))
        }
        if (min_ess_bulk < 400) {
            warn_msg <- c(warn_msg, sprintf(
                "Min bulk ESS is %.1f. Try iter_sampling = %d.",
                min_ess_bulk, rec_iter
            ))
        }
        if (min_ess_tail < 400) {
            warn_msg <- c(warn_msg, sprintf(
                "Min tail ESS is %.1f. Try iter_sampling = %d.",
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

    structure(
        list(
            stanfit = fit,
            dipper_data = prep_data
        ),
        class = "dipper_fit"
    )
}
