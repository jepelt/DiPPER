# Method tests for summary, plot and print. These use the pre-computed
# fit_example object, and variants derived from it for the cases that a single
# stored fit cannot represent. None of these require CmdStan.

test_that("fit_example has the documented structure", {

    data("fit_example", package = "DiPPER")

    expect_s3_class(fit_example, "dipper_fit")
    expect_named(
        fit_example,
        c("draws", "dipper_data", "symmetric", "diagnostics")
    )

    # The fit_example should not contain the CmdStanR object
    expect_false("stanfit" %in% names(fit_example))

    # Diagnostics survive serialisation, and carry everything print() needs
    expect_type(fit_example$diagnostics, "list")
    expect_true(is.numeric(fit_example$diagnostics$num_divergent))
    expect_true(is.numeric(fit_example$diagnostics$num_max_treedepth))
    expect_true(is.numeric(fit_example$diagnostics$n_chains))
    expect_true(is.numeric(fit_example$diagnostics$n_iter_sampling))
    expect_true(is.numeric(fit_example$diagnostics$max_rhat))

    # Draws are a plain matrix, not a posterior draws_matrix, so that the
    # object can be used without the posterior package
    expect_true(is.matrix(fit_example$draws))
    expect_true(is.numeric(fit_example$draws))
    expect_identical(class(fit_example$draws), c("matrix", "array"))

    # The default keep.pars retains beta plus the two prior hyperparameters
    pars <- sub("\\[.*$", "", colnames(fit_example$draws))
    expect_setequal(pars, c("beta", "tau", "nu"))
    expect_equal(sum(pars == "beta"), fit_example$dipper_data$K)

    expect_equal(
        nrow(fit_example$draws),
        fit_example$diagnostics$n_chains *
            fit_example$diagnostics$n_iter_sampling
    )

    # The number of features before prevalence filtering is recorded
    expect_true(fit_example$dipper_data$K_unfiltered >=
                    fit_example$dipper_data$K)
})


test_that("summary.dipper_fit returns a well-formed data.frame", {

    data("fit_example", package = "DiPPER")

    # Summary object should be a data.frame with one row per taxon and the
    # columns taxon, log_or, lwr, upr, significant, and pseudo_q. The column
    # significant should be logical.
    res <- summary(fit_example)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), fit_example$dipper_data$K)
    expect_named(
        res,
        c("taxon", "log_or", "lwr", "upr", "significant", "pseudo_q")
    )
    expect_type(res$significant, "logical")

    # The taxon column comes from the prepared data
    expect_setequal(res$taxon, fit_example$dipper_data$taxa_names)

    # Each credible interval includes the point estimate
    expect_true(all(res$lwr <= res$log_or))
    expect_true(all(res$log_or <= res$upr))

    # The data.frame should be sorted by pseudo_q which should lie in the
    # interval (0, 1).
    expect_false(is.unsorted(res$pseudo_q))
    expect_true(all(res$pseudo_q > 0 & res$pseudo_q <= 1))
})


# Test for probability mass (prob) and the scale (or.scale) for results
# (log(OR) or OR)
test_that("summary.dipper_fit respects prob and or.scale", {

    data("fit_example", package = "DiPPER")

    # The results in the default log(OR) scale (or.scale = log_odds).
    res_95 <- summary(fit_example)
    res_50 <- summary(fit_example, prob = 0.50)

    # A 50% interval is never wider than a 95% one
    expect_true(
        all((res_50$upr - res_50$lwr) <= (res_95$upr - res_95$lwr) + 1e-8)
    )

    # The results expressed in the odds ratio scale.
    res_or <- summary(fit_example, or.scale = "odds_ratio")

    expect_true("odds_ratio" %in% names(res_or))
    expect_false("log_or" %in% names(res_or))

    # ORs should always be positive
    expect_true(all(res_or$odds_ratio > 0))

    # The OR is the exponentiated log(OR)
    expect_equal(res_or$odds_ratio, exp(res_95$log_or))

    # Significance should be independent of the used scale
    expect_equal(res_or$significant, res_95$significant)

    # An invalid interval level is rejected
    expect_error(summary(fit_example, prob = 1.5), "between 0 and 1")
    expect_error(summary(fit_example, prob = c(0.5, 0.9)), "single number")
})


test_that("original.scale undoes the covariate standardization", {

    data("fit_example", package = "DiPPER")

    # fit_example has a categorical variable of interest, which is not scaled.
    # Injecting a scale lets the back-transformation be tested directly.
    scaled <- fit_example
    voi <- scaled$dipper_data$var.of.interest
    scaled$dipper_data$continuous_scales <- stats::setNames(list(2), voi)

    per_sd <- summary(scaled, original.scale = FALSE)
    per_unit <- summary(scaled, original.scale = TRUE)

    # Results per one original unit are the per-SD results divided by the SD
    expect_equal(per_unit$log_or, per_sd$log_or / 2)
    expect_equal(per_unit$lwr, per_sd$lwr / 2)

    # Without a recorded scale the two agree
    res_a <- summary(fit_example, original.scale = TRUE)
    res_b <- summary(fit_example, original.scale = FALSE)
    expect_equal(res_a$log_or, res_b$log_or)
})


test_that("plot.dipper_fit returns a ggplot object", {

    data("fit_example", package = "DiPPER")

    p <- plot(fit_example, show.taxa = "all")

    expect_s3_class(p, "ggplot")

    # The data should have one row per taxon/feature, i.e. the same number of
    # rows as the filtered abundance data.
    expect_equal(nrow(p$data), fit_example$dipper_data$K)

    # Selecting Top-k taxa keeps exactly k taxa
    p_top <- plot(fit_example, show.taxa = 5)
    expect_equal(nrow(p_top$data), 5)

    # An invalid selection is rejected
    expect_error(plot(fit_example, show.taxa = "some"), "should be one of")
    expect_error(plot(fit_example, prob = 0), "between 0 and 1")
})


test_that("plot.dipper_fit returns NULL when nothing is significant", {

    data("fit_example", package = "DiPPER")

    # A very narrow interval makes every credible interval exclude OR = 1;
    # a very wide one makes none of them do so
    none <- fit_example
    none$draws[] <- 0

    expect_message(
        p <- plot(none, show.taxa = "significant"),
        "No significant taxa"
    )
    expect_null(p)

    # verbose = FALSE silences the message but still returns NULL
    expect_silent(p2 <- plot(none, show.taxa = "significant", verbose = FALSE))
    expect_null(p2)
})


test_that("print.dipper_fit reports the fit", {

    data("fit_example", package = "DiPPER")

    expect_output(print(fit_example), "DiPPER Model Fit")
    expect_output(print(fit_example), "Model formula")
    expect_output(print(fit_example), "Posterior draws")
    expect_output(print(fit_example), "MCMC diagnostics")

    # The prior line names the prior and reports its hyperparameters
    expect_output(print(fit_example), "Asymmetric Laplace")
    expect_output(print(fit_example), "tau")
    expect_output(print(fit_example), "nu")

    # The feature line reports how many were kept and how many filtered out
    expect_output(print(fit_example), "Features/taxa")
    expect_output(print(fit_example), "retained")

    # print() returns its argument invisibly
    invisible(capture.output(
        visible <- withVisible(print(fit_example))$visible
    ))
    expect_false(visible)
})


test_that("print.dipper_fit adapts to the contents of the fit", {

    data("fit_example", package = "DiPPER")

    # A symmetric fit has no nu to report
    sym <- fit_example
    sym$symmetric <- TRUE
    sym$draws <- sym$draws[, colnames(sym$draws) != "nu", drop = FALSE]
    expect_output(print(sym), "Symmetric Laplace")
    expect_output(print(sym), "nu fixed at 0.5")

    # A fit made with keep.pars = "beta" has no hyperparameters at all, as do
    # fits saved before they were retained
    beta_only <- fit_example
    beta_only$draws <- beta_only$draws[
        , grepl("^beta\\[", colnames(beta_only$draws)), drop = FALSE]
    out <- capture.output(print(beta_only))
    expect_true(any(grepl("Asymmetric Laplace", out)))
    expect_false(any(grepl("tau", out)))

    # Without K_unfiltered, only the modelled count is shown
    no_kunf <- fit_example
    no_kunf$dipper_data$K_unfiltered <- NULL
    out_k <- capture.output(print(no_kunf))
    expect_false(any(grepl("removed by filtering", out_k)))
    expect_true(any(grepl("Features/taxa", out_k)))
})


test_that("methods fail informatively without usable draws", {

    data("fit_example", package = "DiPPER")

    no_draws <- fit_example
    no_draws$draws <- NULL

    expect_error(summary(no_draws), "no posterior draws")
    expect_error(plot(no_draws), "no posterior draws")

    # A fit whose draws exclude beta cannot be summarised either
    no_beta <- fit_example
    no_beta$draws <- no_beta$draws[
        , !grepl("^beta\\[", colnames(no_beta$draws)), drop = FALSE]

    expect_error(summary(no_beta), "No draws of 'beta'")
    expect_error(plot(no_beta), "No draws of 'beta'")
})


test_that("methods reject objects of the wrong class", {

    expect_error(summary.dipper_fit(list()), "dipper_fit")
    expect_error(plot.dipper_fit(list()), "dipper_fit")
    expect_error(print.dipper_fit(list()), "dipper_fit")
})
