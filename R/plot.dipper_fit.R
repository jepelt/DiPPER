#' Plot method for DiPPER fit objects
#'
#' @param x A dipper_fit object returned by run_dipper.
#' @param prob Numeric. Probability mass for the credible intervals.
#'   Default is 0.95 (95% interval).
#' @param show.taxa Selection criteria for taxa: "significant" (credible
#'   intervals not overlapping OR = 1),
#'   "all", or an integer k (top k taxa ranked by pseudo_q). Default is 10.
#' @param original.scale Logical. If TRUE (default), continuous covariates
#'   are back-transformed to their original scale. If FALSE, results are
#'   presented per one standard deviation increase.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @method plot dipper_fit
#' @export
plot.dipper_fit <- function(x, prob = 0.95, show.taxa = 10,
                            original.scale = TRUE, ...) {
    if (!inherits(x, "dipper_fit")) {
        stop("Input must be a 'dipper_fit' object.")
    }

    d_data <- x$dipper_data
    taxa <- d_data$taxa_names
    var_int <- d_data$var.of.interest
    effect_name <- d_data$design_matrix_cols[1]
    var_levels <- d_data$var_levels

    alpha_level <- 1 - prob
    lower_prob <- alpha_level / 2
    upper_prob <- 1 - (alpha_level / 2)

    # 1. Extract draws for differential prevalence parameters (beta)
    draws <- x$stanfit$draws(variables = "beta", format = "matrix")

    # Back-transform continuous variables to original scale if requested
    scale_factor <- 1
    if (original.scale && var_int %in% names(d_data$continuous_scales)) {
        scale_factor <- d_data$continuous_scales[[var_int]]
        draws <- draws / scale_factor
    }

    # 2. Calculate statistics
    prob_pos <- colMeans(draws > 0)
    prob_neg <- colMeans(draws < 0)
    max_prob <- pmax(prob_pos, prob_neg)
    pseudo_q <- 2 * (1 - max_prob)

    # Exponentiate for the Odds Ratio scale
    est_median <- exp(apply(draws, 2, stats::median))
    ci_lower <- exp(apply(draws, 2, stats::quantile, probs = lower_prob))
    ci_upper <- exp(apply(draws, 2, stats::quantile, probs = upper_prob))

    summ <- data.frame(
        variable = colnames(draws),
        median = est_median,
        lower = ci_lower,
        upper = ci_upper,
        pseudo_q = pseudo_q,
        stringsAsFactors = FALSE
    )
    rownames(summ) <- NULL

    # 3. Map indices to taxa names
    regex_pat <- "^beta\\[([0-9]+)\\]$"
    idx_str <- gsub(regex_pat, "\\1", summ$variable)
    idx <- suppressWarnings(as.integer(idx_str))

    if (all(!is.na(idx)) && max(idx) <= length(taxa)) {
        summ$taxon <- taxa[idx]
    } else {
        summ$taxon <- summ$variable
    }

    # 4. Filtering logic
    plot_df <- summ
    if (is.character(show.taxa) && show.taxa == "significant") {
        plot_df <- plot_df[(plot_df$lower > 1 | plot_df$upper < 1), ]
        if (nrow(plot_df) == 0) {
            message("No significant taxa found. Returning empty plot.")
            return(invisible(NULL))
        }
    } else if (is.numeric(show.taxa)) {
        k <- as.integer(show.taxa)
        plot_df <- plot_df[order(plot_df$pseudo_q), ]
        plot_df <- utils::head(plot_df, k)
    }

    # 5. Finalizing plot data
    plot_df <- plot_df[order(plot_df$median), ]
    plot_df$taxon <- factor(plot_df$taxon, levels = plot_df$taxon)

    # 6. Formulate subtitle and title
    title_str <- paste("Effect of", var_int)
    sub_str <- NULL

    if (!is.null(var_levels) && length(var_levels) >= 2) {
        ref_lvl <- var_levels[1]

        if (effect_name != var_int && grepl(paste0("^", var_int), effect_name)){
            comp_lvl <- sub(paste0("^", var_int), "", effect_name)
        } else {
            comp_lvl <- effect_name
        }

        sub_str <- sprintf(
            "Comparing %s against %s (reference)", comp_lvl, ref_lvl
        )
    }

    # 7. Generate ggplot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = median, y = taxon)) +
        ggplot2::geom_vline(
            xintercept = 1, linetype = "dashed", color = "grey50",
            linewidth = 0.8
        ) +
        ggplot2::geom_pointrange(
            ggplot2::aes(xmin = lower, xmax = upper),
            size = 0.5
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::theme_bw() +
        ggplot2::labs(
            x = "Differential prevalence (OR)",
            y = NULL,
            title = title_str,
            subtitle = sub_str
        )

    return(p)
}
