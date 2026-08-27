# Tests that actually fit a model. These are skipped wherever CmdStan is
# unavailable.

test_that("a simple model can be fitted", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

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

    expect_s3_class(fit, "dipper_fit")
    expect_named(
        fit,
        c("draws", "dipper_data", "symmetric", "diagnostics")
    )

    # The CmdStanR object is not carried along by default
    expect_false("stanfit" %in% names(fit))
    expect_false(fit$symmetric)
    expect_false(fit$dipper_data$is_longitudinal)
    expect_equal(fit$dipper_data$var.of.interest, "Fat")

    # Read depth is added as a covariate by default for count data
    expect_equal(fit$dipper_data$read.depth.var, "log10_read.depth")

    # Diagnostics are captured even when run.diagnostics = FALSE
    expect_type(fit$diagnostics, "list")
    expect_true(is.numeric(fit$diagnostics$num_divergent))
    expect_true(is.numeric(fit$diagnostics$num_max_treedepth))
    expect_equal(fit$diagnostics$n_chains, 1)
    expect_equal(fit$diagnostics$n_iter_sampling, 20)
    expect_null(fit$diagnostics$max_rhat)

    # Draws are a plain matrix: one row per draw, one beta per taxon
    expect_true(is.matrix(fit$draws))
    expect_equal(nrow(fit$draws), 20)
    expect_true(all(grepl("^beta\\[", colnames(fit$draws))))
    expect_equal(ncol(fit$draws), fit$dipper_data$K)

    # print() relies solely on the stored diagnostics
    expect_output(print(fit), "Posterior draws")
})


test_that("dipper runs with random intercepts", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("VatanenT_2016_subset", package = "DiPPER")

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

    expect_s3_class(fit, "dipper_fit")
    expect_true(fit$dipper_data$is_longitudinal)
    expect_equal(fit$dipper_data$id_var, "subject_id")
    expect_equal(
        fit$dipper_data$S,
        length(unique(VatanenT_2016_subset$subject_id))
    )
})


test_that("run.diagnostics = TRUE records convergence statistics", {
    skip_if_not(has_cmdstan(), "CmdStan is not available. Skipping.")

    data("tse_hintikka", package = "DiPPER")

    # Very short chains, so convergence warnings are expected here
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

    vars <- colnames(fit$draws)
    expect_true(any(grepl("^alpha\\[", vars)))
    expect_true(any(grepl("^beta\\[", vars)))

    # summary() must still pick out beta only
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
