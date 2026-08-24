test_that("dipper model runs and returns correct class", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("VatanenT_2016_subset", package = "DiPPER")

    # Run with absolute minimum iterations for testing purposes
    fit <- dipper(
        tse = VatanenT_2016_subset,
        formula = ~ age_point + antibiotics + gender + (1 | subject_id),
        assay.type = "relative_abundance",
        data.type = "relabundance",
        read.depth = FALSE,
        niter = 20,
        niter.warmup = 10,
        chains = 1,
        run.diagnostics = FALSE
    )

    expect_s3_class(fit, "dipper_fit")
    expect_true(!is.null(fit$stanfit))
    expect_false(fit$symmetric)
    expect_true(fit$dipper_data$is_longitudinal)
})
