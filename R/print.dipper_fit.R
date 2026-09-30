#' Posterior summary of a hierarchical prior hyperparameter
#'
#' Returns the posterior median and credible interval of a scalar Stan
#' parameter, formatted for printing, or \code{NULL} if the parameter was not
#' retained in the fit object (see \code{keep.pars}).
#'
#' @param draws Matrix of posterior draws.
#' @param par Name of the parameter (e.g. "tau").
#' @param prob Probability mass of the credible interval.
#' @param digits Number of decimals.
#' @return A string such as `tau 0.08 [0.03, 0.13]`, or NULL.
#'
#' @keywords internal
#' @noRd
.dipper_par_summary <- function(draws, par, prob = 0.90, digits = 2) {

    if (is.null(draws) || !(par %in% colnames(draws))) {
        return(NULL)
    }

    x <- draws[, par]
    qs <- stats::quantile(x, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))

    sprintf(
        "%s %.*f [%.*f, %.*f]",
        par, digits, stats::median(x), digits, qs[1], digits, qs[2]
    )
}


#' Print a summary of a DiPPER model fit
#'
#' @param x A \code{dipper_fit} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The hierarchical prior line indicates the type of the chosen hierarchical
#' prior (symmetric or asymmetric) and its hyperparameters with their 90 percent
#' credible intervals. (For a symmetric fit, \code{nu} is fixed at 0.5.) See the
#' \sQuote{Symmetric or asymmetric prior} section in \code{\link{dipper}} for
#' how to interpret these values and when a symmetric prior may be preferable.
#'
#' @return The object \code{x} invisibly.
#' @method print dipper_fit
#' @export
#' @importFrom stats median quantile
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

    # Prior, with the posterior of its parameters when they were retained
    prior_str <- paste(
        ifelse(x$symmetric, "Symmetric", "Asymmetric"), "Laplace"
    )

    tau_str <- .dipper_par_summary(x$draws, "tau")
    nu_str <- if (x$symmetric) {
        # Only worth stating alongside an actual estimate
        if (is.null(tau_str)) NULL else "nu fixed at 0.5"
    } else {
        .dipper_par_summary(x$draws, "nu")
    }

    par_str <- paste(c(tau_str, nu_str), collapse = ", ")
    if (nzchar(par_str)) {
        prior_str <- paste0(prior_str, " (", par_str, ")")
    }

    cat("Hierarchical prior:  ", prior_str, "\n")

    # Data dimensions. K_unfiltered is absent in fits saved before it was
    # recorded, in which case only the modelled count is shown.
    n_in <- x$dipper_data$K_unfiltered
    n_kept <- x$dipper_data$K

    if (!is.null(n_in) && n_in > n_kept) {
        cat("Features/taxa:       ",
            sprintf("%d retained, %d removed by filtering",
                    n_kept, n_in - n_kept), "\n")
    } else {
        cat("Features/taxa:       ", n_kept, "\n")
    }

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
