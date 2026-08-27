# Method tests using the pre-computed fit_example object.

test_that("summary.dipper_fit returns a well-formed data.frame", {

    data("fit_example", package = "DiPPER")

    res <- summary(fit_example)

    expect_s3_class(res, "data.frame")
    expect_named(
        res,
        c("taxon", "log_or", "lwr", "upr", "significant", "pseudo_q")
    )
    expect_equal(nrow(res), fit_example$dipper_data$K)
    expect_type(res$significant, "logical")

    # The credible interval brackets the point estimate
    expect_true(all(res$lwr <= res$log_or))
    expect_true(all(res$log_or <= res$upr))

    # Sorted by pseudo_q, which is a probability-like quantity
    expect_false(is.unsorted(res$pseudo_q))
    expect_true(all(res$pseudo_q > 0 & res$pseudo_q <= 1))

    # 'significant' agrees with the interval excluding zero
    expect_equal(res$significant, res$lwr > 0 | res$upr < 0)
})


test_that("summary.dipper_fit respects prob and or.scale", {

    data("fit_example", package = "DiPPER")

    res_95 <- summary(fit_example)
    res_50 <- summary(fit_example, prob = 0.50)

    # A 50% interval is never wider than a 95% one
    expect_true(
        all((res_50$upr - res_50$lwr) <= (res_95$upr - res_95$lwr) + 1e-8)
    )

    res_or <- summary(fit_example, or.scale = "odds_ratio")

    expect_true("odds_ratio" %in% names(res_or))
    expect_false("log_or" %in% names(res_or))
    expect_true(all(res_or$odds_ratio > 0))
    expect_equal(res_or$odds_ratio, exp(res_95$log_or))

    # Significance is determined on the log-odds scale either way
    expect_equal(res_or$significant, res_95$significant)
})


test_that("plot.dipper_fit returns a ggplot object", {

    data("fit_example", package = "DiPPER")

    p <- plot(fit_example, show.taxa = "all")

    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), fit_example$dipper_data$K)

    # Top-k selection keeps exactly k taxa
    p_top <- plot(fit_example, show.taxa = 5)
    expect_equal(nrow(p_top$data), 5)
})


test_that("print.dipper_fit works without CmdStan", {

    data("fit_example", package = "DiPPER")

    expect_output(print(fit_example), "DiPPER Model Fit")
    expect_output(print(fit_example), "Model formula")
    expect_output(print(fit_example), "Posterior draws")
    expect_output(print(fit_example), "MCMC diagnostics")

    # print() returns its argument invisibly
    expect_identical(withVisible(print(fit_example))$visible, FALSE)
})


test_that("fit_example has the documented structure", {

    data("fit_example", package = "DiPPER")

    expect_s3_class(fit_example, "dipper_fit")
    expect_named(
        fit_example,
        c("draws", "dipper_data", "symmetric", "diagnostics")
    )

    # No CmdStanR object, so the fit loads on machines without cmdstanr
    expect_false("stanfit" %in% names(fit_example))

    # Diagnostics survive serialisation, and carry everything print() needs
    expect_type(fit_example$diagnostics, "list")
    expect_true(is.numeric(fit_example$diagnostics$num_divergent))
    expect_true(is.numeric(fit_example$diagnostics$num_max_treedepth))
    expect_true(is.numeric(fit_example$diagnostics$n_chains))
    expect_true(is.numeric(fit_example$diagnostics$n_iter_sampling))

    # Draws are a plain matrix holding only beta
    expect_true(is.matrix(fit_example$draws))
    expect_true(is.numeric(fit_example$draws))
    expect_true(all(grepl("^beta\\[", colnames(fit_example$draws))))
    expect_equal(ncol(fit_example$draws), fit_example$dipper_data$K)
    expect_equal(
        nrow(fit_example$draws),
        fit_example$diagnostics$n_chains *
            fit_example$diagnostics$n_iter_sampling
    )
})


test_that("methods fail informatively without draws", {

    data("fit_example", package = "DiPPER")

    no_draws <- fit_example
    no_draws$draws <- NULL

    expect_error(summary(no_draws), "no posterior draws")
    expect_error(plot(no_draws), "no posterior draws")
})


test_that("methods reject objects of the wrong class", {

    expect_error(summary.dipper_fit(list()), "dipper_fit")
    expect_error(plot.dipper_fit(list()), "dipper_fit")
})
