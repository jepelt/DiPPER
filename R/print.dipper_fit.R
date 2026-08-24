#' Print a summary of a DiPPER model fit
#'
#' @param x A \code{dipper_fit} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} invisibly.
#' @method print dipper_fit
#' @export
print.dipper_fit <- function(x, ...) {

    cat("DiPPER Model Fit\n")
    cat("----------------\n")

    # Formula
    form_str <- paste(deparse(x$dipper_data$formula), collapse = " ")
    form_str <- gsub("\\s+", " ", form_str)
    cat("Model formula:       ", form_str, "\n")

    # Prior
    prior_str <- ifelse(x$symmetric, "Symmetric", "Asymmetric")
    cat("Hierarchical prior:  ", prior_str, "Laplace\n")

    # Data dimensions
    cat("Features modeled:    ", x$dipper_data$K, "\n")
    cat("Samples:             ", x$dipper_data$N, "\n")

    # Posterior samples
    draws_dim <- dim(x$stanfit$draws())
    n_iters <- draws_dim[1]
    n_chains <- draws_dim[2]
    tot_draws <- n_iters * n_chains
    cat("Posterior draws:     ", tot_draws,
        sprintf(" (%d chains x %d iterations)", n_chains, n_iters), "\n")

    # Diagnostics (instantaneous from sampler metadata)
    diag_sum <- x$stanfit$diagnostic_summary()
    divs <- sum(diag_sum$num_divergent)
    treedepths <- sum(diag_sum$num_max_treedepth)

    if (divs == 0 && treedepths == 0) {
        cat("MCMC diagnostics:    ",
            "No divergent transitions or max treedepth hits\n")
    } else {
        cat("MCMC diagnostics:    ", divs, "divergent transitions,",
            treedepths, "max treedepth hits\n")
    }

    # Computation time (CmdStanR wall time for sampling)
    if (!is.null(x$stanfit$time()$total)) {
        run_time <- round(x$stanfit$time()$total, 2)
        cat("MCMC sampling time:  ", run_time, "seconds\n")
    }

    invisible(x)
}
