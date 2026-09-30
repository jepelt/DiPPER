# Tests that actually fit a model. These are skipped wherever CmdStan is
# unavailable.

# Shared arguments for the tse_hintikka fits below.
hintikka_args <- function(...) {

    e <- new.env()
    utils::data("tse_hintikka", package = "DiPPER", envir = e)

    defaults <- list(
        tse = e$tse_hintikka,
        assay.type = "counts",
        data.type = "counts",
        formula = ~ Fat,
        read.depth = TRUE,
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE,
        verbose = FALSE
    )

    # Arguments to replace the defaults.
    overrides <- list(...)
    c(defaults[setdiff(names(defaults), names(overrides))], overrides)
}


test_that("a simple model can be fitted", {
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    fit <- do.call(dipper, hintikka_args(formula = ~ Fat + XOS))

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

    # log10 read depth was added as a covariate, as requested
    expect_equal(fit$dipper_data$read.depth.var, "log10_read_depth")

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

    # Draws are a plain matrix, not a posterior draws_matrix, so that the fit
    # object does not depend on the posterior package
    expect_identical(class(fit$draws), c("matrix", "array"))
    expect_setequal(names(attributes(fit$draws)), c("dim", "dimnames"))

    # By default the differential prevalence parameters and the two prior
    # hyperparameters are retained
    pars <- sub("\\[.*$", "", colnames(fit$draws))
    expect_setequal(pars, c("beta", "tau", "nu"))
    expect_equal(sum(pars == "beta"), fit$dipper_data$K)

    # niter - niter.warmup = 40 - 20 = 20 draws are retained
    expect_equal(nrow(fit$draws), 20)

    # print() finds the chain and iteration counts in the stored diagnostics
    expect_output(print(fit), "Posterior draws")
})


test_that("dipper runs with random intercepts", {
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    data("VatanenT_2016_subset", package = "DiPPER")

    fit <- dipper(
        tse = VatanenT_2016_subset,
        assay.type = "relative_abundance",
        data.type = "relabundance",
        formula = ~ age_point + antibiotics + gender + (1 | subject_id),
        read.depth = FALSE,
        niter = 40,
        niter.warmup = 20,
        chains = 1,
        cores = 1,
        run.diagnostics = FALSE,
        verbose = FALSE
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
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    # Very short chains, so convergence warnings are expected here and are
    # suppressed to keep the test output readable.
    fit <- suppressWarnings(do.call(
        dipper, hintikka_args(chains = 2, run.diagnostics = TRUE)
    ))

    # Now the diagnostics should include max R-hat, min ESS bulk and min ESS
    # tail. By default, diagnostics level is "basic".
    expect_true(is.numeric(fit$diagnostics$max_rhat))
    expect_true(is.numeric(fit$diagnostics$min_ess_bulk))
    expect_true(is.numeric(fit$diagnostics$min_ess_tail))
    expect_equal(fit$diagnostics$diagnostics_level, "basic")

    # print() reports them on the convergence line
    expect_output(print(fit), "max R-hat")
})


test_that("keep.pars controls which draws are retained", {
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    # keep.pars = NULL retains every model parameter, not just the ones the
    # methods need. Besides beta and the intercept alpha, the read depth
    # covariate makes beta_cov present as well.
    fit_all <- do.call(dipper, hintikka_args(keep.pars = NULL))
    vars <- colnames(fit_all$draws)
    expect_true(any(grepl("^beta\\[", vars)))
    expect_true(any(grepl("^alpha\\[", vars)))
    expect_true(any(grepl("^beta_cov\\[", vars)))

    # summary() should still pick out beta only (i.e. only one row per
    # taxonomic feature).
    expect_equal(nrow(summary(fit_all)), fit_all$dipper_data$K)

    # Only beta parameters retained
    fit_beta <- do.call(dipper, hintikka_args(keep.pars = "beta"))
    expect_setequal(sub("\\[.*$", "", colnames(fit_beta$draws)), "beta")
})


test_that("keep.stanfit = TRUE returns the CmdStanR object", {
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    fit <- do.call(dipper, hintikka_args(keep.stanfit = TRUE))

    expect_true("stanfit" %in% names(fit))
    expect_true(inherits(fit$stanfit, "CmdStanMCMC"))

    # The draws matrix is present either way
    expect_true(is.matrix(fit$draws))
})


test_that("symmetric = TRUE selects the symmetric prior", {
    skip_if_not(has_dipper_models(),
                "Compiled DiPPER Stan models are not available. Skipping.")

    fit <- do.call(dipper, hintikka_args(symmetric = TRUE))

    expect_true(fit$symmetric)

    # nu does not exist in the symmetric models, so it is dropped from
    # keep.pars automatically rather than being requested from the sampler
    pars <- sub("\\[.*$", "", colnames(fit$draws))
    expect_setequal(pars, c("beta", "tau"))
    expect_false("nu" %in% pars)

    # print() says so instead of reporting an estimate
    expect_output(print(fit), "nu fixed at 0.5")

    # summary() and plot() work regardless of the prior
    expect_s3_class(summary(fit), "data.frame")
})
