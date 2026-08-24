test_that("summary and plot methods work", {
    skip_if_not(has_cmdstan(), "CmdStan not available. Skipping.")

    data("VatanenT_2016_subset", package = "DiPPER")

    # Generate a minimal fit object
    fit <- dipper(
        tse = VatanenT_2016_subset,
        formula = ~ age_point,
        assay.type = "relative_abundance",
        data.type = "relabundance",
        read.depth = FALSE,
        niter = 20,
        niter.warmup = 10,
        chains = 1,
        run.diagnostics = FALSE
    )

    # Test summary.dipper_fit
    res <- summary(fit)
    expect_s3_class(res, "data.frame")
    expect_true("taxon" %in% names(res))

    # Test plot.dipper_fit
    p <- plot(fit)
    expect_s3_class(p, "ggplot")
})
