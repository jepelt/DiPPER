# Tests that actually fit a model. These are skipped wherever CmdStan is
# unavailable.

test_that("a simple model can be fitted", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    # Note: run.diagnostics is set to FALSE for speed, and to avoid convergence
    # warnings that the (deliberately) short MCMC chains would trigger.
    fit <- dipper(
        tse = tse_hintikka,
        formula = ~ Fat + XOS,
        assay.type = "counts",
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE
    )

    # Expected output object from dipper()
    expect_s3_class(fit, "dipper_fit")
    expect_named(
        fit,
        c("draws", "dipper_data", "symmetric", "diagnostics")
    )

    # The CmdStanR object is not carried along by default
    expect_false("stanfit" %in% names(fit))

    # Asymmetric prior is used by default
    expect_false(fit$symmetric)

    # This example does not have a random intercept (i.e. not longitudinal)
    expect_false(fit$dipper_data$is_longitudinal)

    # The variable of interest is Fat
    expect_equal(fit$dipper_data$var.of.interest, "Fat")

    # (log10) Read depth is added as a covariate by default for count data
    expect_equal(fit$dipper_data$read.depth.var, "log10_read.depth")

    # Some diagnostics (divergent transitions, exceeding max treedepth,
    # number of MCMC chains, and number of retained iterations per chain) are
    # captured even when run.diagnostics = FALSE. However, max R-hat is not
    # computed in this case.
    expect_type(fit$diagnostics, "list")
    expect_true(is.numeric(fit$diagnostics$num_divergent))
    expect_true(is.numeric(fit$diagnostics$num_max_treedepth))
    expect_equal(fit$diagnostics$n_chains, 1)
    expect_equal(fit$diagnostics$n_iter_sampling, 20)
    expect_null(fit$diagnostics$max_rhat)

    # Draws are a matrix (one row per draw, one beta per taxon)
    expect_true(is.matrix(fit$draws))
    expect_true(all(grepl("^beta\\[", colnames(fit$draws))))
    expect_equal(ncol(fit$draws), fit$dipper_data$K)

    # niter - niter.warmup = 40 - 20 = 20 draws are retained
    expect_equal(nrow(fit$draws), 20)

    # print() finds the chain and iteration counts in the stored diagnostics
    expect_output(print(fit), "Posterior draws")
})


test_that("dipper runs with random intercepts", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("VatanenT_2016_subset", package = "DiPPER")

    # Note: run.diagnostics is set to FALSE for speed, and to avoid convergence
    # warnings that the (deliberately) short MCMC chains would trigger.
    fit <- dipper(
        tse = VatanenT_2016_subset,
        formula = ~ age_point + antibiotics + gender + (1 | subject_id),
        assay.type = "relative_abundance",
        data.type = "relabundance",
        read.depth = FALSE,
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE
    )

    # Expected output object from dipper() when random intercepts are included
    expect_s3_class(fit, "dipper_fit")
    expect_true(fit$dipper_data$is_longitudinal)
    expect_equal(fit$dipper_data$id_var, "subject_id")

    # Check that the number of unique subjects matches the S parameter
    expect_equal(
        fit$dipper_data$S,
        length(unique(VatanenT_2016_subset$subject_id))
    )
})


test_that("run.diagnostics = TRUE records convergence statistics", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    # Very short chains, so convergence warnings are expected here.
    # Consequently, suppressWarnings() is used to avoid printing convergence
    # warnings to the output.
    fit <- suppressWarnings(dipper(
        tse = tse_hintikka,
        formula = ~ Fat,
        assay.type = "counts",
        niter = 40,
        niter.warmup = 20,
        chains = 2,
        cores = 1,
        run.diagnostics = TRUE
    ))

    # Now the diagnostics should include max R-hat, min ESS bulk and min ESS
    # tail. By default, diagnostics level is "basic".
    expect_true(is.numeric(fit$diagnostics$max_rhat))
    expect_true(is.numeric(fit$diagnostics$min_ess_bulk))
    expect_true(is.numeric(fit$diagnostics$min_ess_tail))
    expect_equal(fit$diagnostics$diagnostics_level, "basic")
})


test_that("keep.pars = NULL retains all parameters", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    fit <- dipper(
        tse = tse_hintikka,
        formula = ~ Fat,
        assay.type = "counts",
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE,
        keep.pars = NULL
    )

    # With keep.pars = NULL the draws contain all model parameters, not just
    # the differential prevalence parameters (beta). In addition to the
    # variable(s) given in the formula (and intercept alpha), the model by
    # default automatically includes log10 read depth as covariate, so beta_cov
    # is present as well.
    vars <- colnames(fit$draws)
    expect_true(any(grepl("^beta\\[", vars)))
    expect_true(any(grepl("^alpha\\[", vars)))
    expect_true(any(grepl("^beta_cov\\[", vars)))

    # summary() should still pick out beta only (i.e. only one row per
    # taxonomic feature).
    expect_equal(nrow(summary(fit)), fit$dipper_data$K)
})


test_that("keep.stanfit = TRUE returns the CmdStanR object", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    fit <- dipper(
        tse = tse_hintikka,
        formula = ~ Fat,
        assay.type = "counts",
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE,
        keep.stanfit = TRUE
    )

    expect_true("stanfit" %in% names(fit))
    expect_true(inherits(fit$stanfit, "CmdStanMCMC"))

    # The draws matrix is present either way
    expect_true(is.matrix(fit$draws))
})


test_that("symmetric = TRUE selects the symmetric prior", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    fit <- dipper(
        tse = tse_hintikka,
        formula = ~ Fat,
        assay.type = "counts",
        symmetric = TRUE,
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE
    )

    expect_true(fit$symmetric)

    # summary() and plot() work regardless of the prior
    expect_s3_class(summary(fit), "data.frame")
})
