#' Summarize DiPPER results for the variable of interest
#'
#' This function computes and summarizes the posterior quantities (differential
#' prevalence estimate, credible interval, pseudo q-value) for the variable of
#' interest for each taxonomic feature.
#'
#' @param object A \code{dipper_fit} object returned by \code{\link{dipper}}.
#' @param prob Numeric. Probability mass for the credible interval (analogous to
#'   confidence level in frequentist statistics). Default is 0.95 (i.e.,
#'   95 percent credible intervals are produced).
#' @param scale Character string indicating the scale of the estimates. Either
#'   \code{"log_odds"} (default) or \code{"odds_ratio"}.
#' @param original.scale Logical. If \code{TRUE} (default), continuous
#'   covariates are back-transformed to their original scale. If \code{FALSE},
#'   results are presented per one standard deviation increase. Ignored for
#'   categorical variables.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The function computes the median posterior estimates and credible interval
#' boundaries (based on the \code{(1 - prob) / 2} and \code{1 - (1 - prob) / 2}
#' posterior quantiles) from the posterior draws of the differential prevalence
#' parameters (beta) of the variable of interest.
#'
#' Note that as the boundaries of the credible intervals depend on the
#' \code{prob} level, so does the \code{significant} logical flag.
#'
#' \code{max_prob} is the maximum of the posterior probabilities that the
#' parameter is strictly positive or strictly negative. Consequently, it
#' indicates the highest "significance level" at which the result would be
#' considered "statistically significant" (i.e., the credible interval
#' excluding zero). It thus bears some resemblance to a multiplicity-adjusted
#' two-sided \emph{p}-value (\emph{q}-value), but it does not have any strict
#' statistical interpretation. Due to the finite number of posterior draws,
#' zero values are replaced with \code{0.5 * (1 / n_draws)}.
#'
#' @return A \code{data.frame} containing the summarized results for each
#' taxon, sorted by \code{pseudo_q} (smallest to largest). The columns are:
#' \describe{
#'   \item{taxon}{Name of the taxonomic feature.}
#'   \item{log_or / odds_ratio}{Median posterior estimate on the chosen scale.}
#'   \item{lwr}{Lower bound of the credible interval.}
#'   \item{upr}{Upper bound of the credible interval.}
#'   \item{significant}{Logical. \code{TRUE} if the credible interval on the
#'   log-odds scale does not include zero (or on the odds ratio scale does not
#'   include one).}
#'   \item{pseudo_q}{Bayesian pseudo q-value.}
#' }
#'
#' @method summary dipper_fit
#' @export
#' @importFrom stats median quantile
#'
#' @examples
#' library(DiPPER)
#'
#' # Load pre-run model fit for the example dataset (tse_hintikka)
#' data("fit_example")
#'
#' # Summarize on log-odds scale with default 95 percent interval
#' res_log <- summary(fit_example)
#' head(res_log)
#'
#' # Summarize on odds ratio scale with 90 percent interval
#' res_or <- summary(fit_example, prob = 0.90, scale = "odds_ratio")
#' head(res_or)
summary.dipper_fit <- function(object,
                               prob = 0.95,
                               scale = c("log_odds", "odds_ratio"),
                               original.scale = TRUE,
                               ...) {

    scale <- match.arg(scale)

    if (!inherits(object, "dipper_fit")) {
        stop("Input must be a 'dipper_fit' object.")
    }

    alpha <- 1 - prob
    lower_prob <- alpha / 2
    upper_prob <- 1 - (alpha / 2)

    draws <- object$stanfit$draws("beta", format = "matrix")
    taxa_names <- object$dipper_data$taxa_names
    var_int <- object$dipper_data$var.of.interest

    if (ncol(draws) != length(taxa_names)) {
        stop("Mismatch between the number of parameters and taxa names.")
    }

    # Back-transform continuous variables to original scale if requested
    if (original.scale &&
        var_int %in% names(object$dipper_data$continuous_scales)) {
        scale_factor <- object$dipper_data$continuous_scales[[var_int]]
        draws <- draws / scale_factor
    }

    est_median <- apply(draws, 2, stats::median)
    ci_lower <- apply(draws, 2, stats::quantile, probs = lower_prob)
    ci_upper <- apply(draws, 2, stats::quantile, probs = upper_prob)

    prob_pos <- colMeans(draws > 0)
    prob_neg <- colMeans(draws < 0)
    max_prob <- pmax(prob_pos, prob_neg)

    pseudo_q <- 2 * (1 - max_prob)

    # Prevent pseudo_q from being exactly 0 to avoid log(0) in volcano plots
    n_draws <- nrow(draws)
    pseudo_q[pseudo_q == 0] <- 0.5 * (1 / n_draws)

    # Significance is determined on the log-odds scale
    significant <- (ci_lower > 0) | (ci_upper < 0)

    if (scale == "odds_ratio") {
        est_median <- exp(est_median)
        ci_lower <- exp(ci_lower)
        ci_upper <- exp(ci_upper)
    }

    res <- data.frame(
        taxon = taxa_names,
        estimate = est_median,
        lwr = ci_lower,
        upr = ci_upper,
        significant = significant,
        pseudo_q = pseudo_q,
        stringsAsFactors = FALSE
    )

    est_name <- ifelse(scale == "odds_ratio", "odds_ratio", "log_or")
    colnames(res)[2] <- est_name

    # Sort by pseudo_q (smallest pseudo q-value at the top)
    res <- res[order(res$pseudo_q), ]
    rownames(res) <- NULL

    return(res)
}
