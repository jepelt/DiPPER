#' Print a summary of a DiPPER model fit
#'
#' @param x A \code{dipper_fit} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} invisibly.
#' @method print dipper_fit
#' @export
#'
#' @examples
#' # Load pre-run model fit for the example dataset (tse_hintikka)
#' data("fit_example")
#'
#' print(fit_example)
print.dipper_fit <- function(x, ...) {

    if (!inherits(x, "dipper_fit")) {
        stop("Input must be a 'dipper_fit' object.", call. = FALSE)
    }

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

    # Everything below is read from the stored diagnostics rather than from
    # the CmdStanMCMC object, so that it also works for saved fits.
    diag <- x$diagnostics

    if (!is.null(diag)) {

        # Posterior samples
        if (!is.null(diag$n_chains) && !is.null(diag$n_iter_sampling)) {
            tot_draws <- diag$n_chains * diag$n_iter_sampling
            cat("Posterior draws:     ", tot_draws,
                sprintf(" (%d chains x %d iterations)",
                        diag$n_chains, diag$n_iter_sampling), "\n")
        }

        if (diag$num_divergent == 0 && diag$num_max_treedepth == 0) {
            cat("MCMC diagnostics:    ",
                "No divergent transitions or max treedepth hits\n")
        } else {
            cat("MCMC diagnostics:    ", diag$num_divergent,
                "divergent transitions,",
                diag$num_max_treedepth, "max treedepth hits\n")
        }

        if (!is.null(diag$max_rhat) && is.finite(diag$max_rhat)) {
            cat("Convergence:         ",
                sprintf(
                    "max R-hat %.3f, min bulk ESS %.0f, min tail ESS %.0f",
                    diag$max_rhat, diag$min_ess_bulk, diag$min_ess_tail
                ),
                "\n")
        }

        if (!is.null(diag$time_total)) {
            cat("MCMC sampling time:  ",
                round(diag$time_total, 2), "seconds\n")
        }
    }

    invisible(x)
}
