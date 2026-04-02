#' Summarize the results from a DiPPER model fit
#'
#' @param fit A dipper_fit object returned by run_dipper.
#' @param prob Probability mass for the credible interval (default 0.95).
#' @param scale Character string, either "log_odds" or "odds_ratio".
#'
#' @return A data.frame containing the summarized results for each taxon.
#' @export
#' @importFrom stats median quantile
summary_dipper <- function(fit,
                           prob = 0.95,
                           scale = c("log_odds", "odds_ratio")) {

    scale <- match.arg(scale)

    if (!inherits(fit, "dipper_fit")) {
        stop("Input must be a 'dipper_fit' object.")
    }

    alpha <- 1 - prob
    lower_prob <- alpha / 2
    upper_prob <- 1 - (alpha / 2)

    draws <- fit$stanfit$draws("beta", format = "matrix")

    est_median <- apply(draws, 2, stats::median)
    ci_lower <- apply(draws, 2, stats::quantile, probs = lower_prob)
    ci_upper <- apply(draws, 2, stats::quantile, probs = upper_prob)

    prob_pos <- colMeans(draws > 0)
    prob_neg <- colMeans(draws < 0)
    max_prob <- pmax(prob_pos, prob_neg)

    # Multilevel shrinkage intrinsically corrects for multiple testing.
    # We output a pseudo q-value based on the directional probability.
    pseudo_q <- 2 * (1 - max_prob)

    significant <- (ci_lower > 0) | (ci_upper < 0)

    if (scale == "odds_ratio") {
        est_median <- exp(est_median)
        ci_lower <- exp(ci_lower)
        ci_upper <- exp(ci_upper)
    }

    res <- data.frame(
        taxon = fit$dipper_data$taxa_names,
        estimate = est_median,
        lwr = ci_lower,
        upr = ci_upper,
        pseudo_q = pseudo_q,
        significant = significant,
        stringsAsFactors = FALSE
    )

    est_name <- ifelse(scale == "odds_ratio", "odds_ratio", "log_or")
    colnames(res)[2] <- est_name

    # Sort by pseudo_q (smallest pseudo q-value at the top)
    res <- res[order(res$pseudo_q), ]
    rownames(res) <- NULL

    return(res)
}
