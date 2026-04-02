#' Create a forest plot from a DiPPER fit
#'
#' @param fit A dipper_fit object returned by run_dipper.
#' @param prob Numeric. Probability mass for the credible intervals.
#'   Default is 0.95 (95% interval).
#' @param show_taxa Selection criteria for taxa: "significant" (default),
#'   "all", or an integer k (top k taxa ranked by pseudo_q).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_dipper <- function(fit, prob = 0.95, show_taxa = "significant") {
    if (!inherits(fit, "dipper_fit")) {
        stop("Input must be a 'dipper_fit' object.")
    }

    d_data <- fit$dipper_data
    taxa <- d_data$taxa_names
    var_int <- d_data$var_of_interest
    effect_name <- d_data$design_matrix_cols[1]
    var_levels <- d_data$var_levels

    alpha_level <- 1 - prob
    lower_prob <- alpha_level / 2
    upper_prob <- 1 - (alpha_level / 2)

    # 1. Extract draws and calculate statistics manually
    draws <- fit$stanfit$draws(variables = "beta", format = "matrix")

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

    # 2. Map indices to taxa names
    regex_pat <- "^beta\\[([0-9]+).*"
    idx_str <- gsub(regex_pat, "\\1", summ$variable)
    idx <- suppressWarnings(as.integer(idx_str))

    if (all(!is.na(idx)) && max(idx) <= length(taxa)) {
        summ$taxon <- taxa[idx]
    } else {
        summ$taxon <- summ$variable
    }

    # 3. Filtering logic
    plot_df <- summ
    if (show_taxa == "significant") {
        plot_df <- plot_df[(plot_df$lower > 1 | plot_df$upper < 1), ]
        if (nrow(plot_df) == 0) {
            message("No significant taxa found. Returning empty plot.")
            message("Decrease prob or set e.g. show_taxa=10 to show some results.")
            return(invisible(NULL))
        }
    } else if (is.numeric(show_taxa)) {
        k <- as.integer(show_taxa)
        plot_df <- plot_df[order(plot_df$pseudo_q), ]
        plot_df <- utils::head(plot_df, k)
    }

    # 4. Finalizing plot data
    plot_df <- plot_df[order(plot_df$median), ]
    plot_df$taxon <- factor(plot_df$taxon, levels = plot_df$taxon)

    # 5. Formulate subtitle
    if (!is.null(var_levels) && length(var_levels) >= 2) {
        ref_lvl <- var_levels[1]

        if (effect_name != var_int && grepl(paste0("^", var_int), effect_name)) {
            comp_lvl <- sub(paste0("^", var_int), "", effect_name)
        } else {
            comp_lvl <- effect_name
        }

        sub_str <- sprintf("Comparing %s against %s (reference)",
                           comp_lvl, ref_lvl)
    } else {
        sub_str <- paste("Coefficient:", effect_name)
    }

    title_str <- paste("Effect of", var_int)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = median, y = taxon)) +
        ggplot2::geom_vline(
            xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8
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
