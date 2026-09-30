# Unit tests for the internal helpers of run_dipper().

# Minimal prep.data list with the elements that the helpers read.
make_prep <- function(P = 2, longitudinal = FALSE, read.depth.var = NULL) {
    cols <- "x_interest"
    if (P > 1) {
        cols <- c(cols, paste0("cov", seq_len(P - 1)))
        if (!is.null(read.depth.var)) {
            cols[P] <- read.depth.var
        }
    }

    prep <- list(
        y = matrix(c(0, 1, 1, 0, 1, 0), nrow = 2),
        X = matrix(0, nrow = 3, ncol = P),
        N = 3,
        K = 2,
        P = P,
        design_matrix_cols = cols,
        read.depth.var = read.depth.var,
        is_longitudinal = longitudinal,
        S = 0
    )

    if (longitudinal) {
        prep$S <- 2
        prep$subj <- c(1, 1, 2)
    }

    prep
}


test_that(".dipper_validate_run_args catches invalid arguments", {

    prep <- make_prep()

    # Valid arguments pass silently, including the default keep.pars
    expect_silent(.dipper_validate_run_args(prep, "beta", FALSE, 2000, 1000))
    expect_silent(.dipper_validate_run_args(prep, NULL, TRUE, 2000, 1000))
    expect_silent(
        .dipper_validate_run_args(prep, c("beta", "tau", "nu"), FALSE,
                                  2000, 1000)
    )

    # prep.data must contain the elements created by prep_dipper_data()
    expect_error(
        .dipper_validate_run_args(prep[c("y", "X")], "beta", FALSE, 2000, 1000),
        "Invalid prep.data"
    )

    # keep.pars must be a character vector or NULL
    expect_error(
        .dipper_validate_run_args(prep, 1, FALSE, 2000, 1000),
        "'keep.pars' must be"
    )

    # keep.stanfit must be a single logical value
    expect_error(
        .dipper_validate_run_args(prep, "beta", "yes", 2000, 1000),
        "'keep.stanfit' must be"
    )
    expect_error(
        .dipper_validate_run_args(prep, "beta", c(TRUE, FALSE), 2000, 1000),
        "'keep.stanfit' must be"
    )

    # At least one sampling iteration must remain after warmup
    expect_error(
        .dipper_validate_run_args(prep, "beta", FALSE, 1000, 1000),
        "strictly greater"
    )
})


test_that("run_dipper() validates arguments before looking for CmdStan", {

    # These errors are raised regardless of whether CmdStan is installed
    expect_error(run_dipper(list()), "Invalid prep.data")
    expect_error(
        run_dipper(make_prep(), niter = 10, niter.warmup = 10),
        "strictly greater"
    )
    expect_error(
        run_dipper(make_prep(), diagnostics.level = "extreme"),
        "should be one of"
    )
})


test_that("dipper() reports argument errors before sampling", {

    data("tse_hintikka", package = "DiPPER")

    base_args <- list(tse = tse_hintikka, assay.type = "counts",
                      data.type = "counts", formula = ~ Fat,
                      read.depth = TRUE, verbose = FALSE)

    # The data are prepared, but run_dipper() stops at argument validation,
    # so no CmdStan is needed
    expect_error(
        do.call(dipper, c(base_args, list(niter = 100, niter.warmup = 100))),
        "strictly greater"
    )
    expect_error(
        do.call(dipper, c(base_args, list(keep.pars = 1))),
        "'keep.pars' must be"
    )
    expect_error(
        do.call(dipper, c(base_args, list(keep.stanfit = "yes"))),
        "'keep.stanfit' must be"
    )

    # Arguments matched in dipper() itself
    expect_error(
        do.call(dipper, c(base_args, list(diagnostics.level = "extreme"))),
        "should be one of"
    )
})


test_that("dipper() passes the data validation errors through", {

    data("tse_hintikka", package = "DiPPER")

    # These come from prep_dipper_data(), so that both entry points report the
    # same problems in the same way
    expect_error(
        dipper(tse = tse_hintikka, assay.type = "counts",
               data.type = "counts", read.depth = TRUE, verbose = FALSE),
        "Please specify 'formula'"
    )
    expect_error(
        dipper(tse = tse_hintikka, assay.type = "counts", formula = ~ Fat,
               verbose = FALSE),
        "Please specify 'data.type'"
    )
    expect_error(
        dipper(tse = tse_hintikka, assay.type = "counts",
               data.type = "counts", formula = ~ Fat, verbose = FALSE),
        "Please specify 'read.depth' for count data"
    )
    expect_error(
        dipper(tse = tse_hintikka, assay.type = "counts",
               data.type = "clr", formula = ~ Fat, verbose = FALSE),
        "should be one of"
    )
})


test_that(".dipper_cov_priors assigns covariate and read depth priors", {

    # Only the variable of interest: no covariate priors
    none <- .dipper_cov_priors(make_prep(P = 1), 1, 2, 2)
    expect_identical(none, list(mean = numeric(0), sd = numeric(0)))

    # Covariates without read depth get the default covariate prior
    plain <- .dipper_cov_priors(make_prep(P = 3), 1.5, 2, 2)
    expect_equal(plain$mean, c(0, 0))
    expect_equal(plain$sd, c(1.5, 1.5))

    # The read depth covariate gets its own prior
    rd <- .dipper_cov_priors(
        make_prep(P = 3, read.depth.var = "log10_read_depth"),
        prior.cov.sd = 1.5, prior.reads.mean = 2.5, prior.reads.sd = 3
    )
    expect_equal(rd$mean, c(0, 2.5))
    expect_equal(rd$sd, c(1.5, 3))
})


test_that(".dipper_cov_priors matches the read depth name exactly", {

    # The read depth name is matched exactly, not as a prefix or pattern
    prep <- make_prep(P = 3, read.depth.var = "log10_read_depth")
    prep$design_matrix_cols[2] <- "log10_read_depth_raw"

    res <- .dipper_cov_priors(prep, 1, 2, 2)
    expect_equal(res$mean, c(0, 2))
})


test_that(".dipper_stan_setup selects the model and builds Stan data", {

    setup <- function(prep, symmetric) {
        .dipper_stan_setup(
            prep.data = prep, symmetric = symmetric,
            prior.alpha.sd = 4, prior.tau.sd = 1, prior.nu.sd = 0.05,
            prior.cov.sd = 1, prior.reads.mean = 2, prior.reads.sd = 2,
            prior.sigma.subj = 0.5
        )
    }

    cross <- make_prep()
    long <- make_prep(longitudinal = TRUE)

    # Model name for each combination of design and prior
    expect_equal(setup(cross, FALSE)$model_name, "dipper_dp_asym")
    expect_equal(setup(cross, TRUE)$model_name, "dipper_dp_sym")
    expect_equal(setup(long, FALSE)$model_name, "dipper_dp_long_asym")
    expect_equal(setup(long, TRUE)$model_name, "dipper_dp_long_sym")

    # Common Stan data are copied from prep.data and the prior arguments
    sd_cross <- setup(cross, FALSE)$stan_data
    expect_equal(sd_cross$N, cross$N)
    expect_equal(sd_cross$K, cross$K)
    expect_equal(sd_cross$P, cross$P)
    expect_identical(sd_cross$y, cross$y)
    expect_identical(sd_cross$X, cross$X)
    expect_equal(sd_cross$prior_alpha_mean, 0)
    expect_equal(sd_cross$prior_alpha_sd, 4)
    expect_equal(sd_cross$prior_tau_sd, 1)

    # Covariate priors are arrays, as required by Stan for vectors
    expect_true(is.array(sd_cross$prior_cov_mean))
    expect_true(is.array(sd_cross$prior_cov_sd))

    # The asymmetry prior is only included for asymmetric models
    expect_equal(sd_cross$prior_nu_sd, 0.05)
    expect_null(setup(cross, TRUE)$stan_data$prior_nu_sd)

    # Cross-sectional models have no subject structure
    expect_null(sd_cross$subj)
    expect_null(sd_cross$prior_sigma_subj)

    # Longitudinal models carry the subject indices and their prior
    sd_long <- setup(long, TRUE)$stan_data
    expect_equal(sd_long$S, long$S)
    expect_equal(sd_long$subj, long$subj)
    expect_equal(sd_long$prior_sigma_subj, 0.5)
})


test_that(".dipper_progress_settings interprets print.progress", {

    on <- function(refresh) list(refresh = refresh, show_messages = TRUE)
    off <- list(refresh = 0L, show_messages = FALSE)

    expect_equal(.dipper_progress_settings(200, TRUE), on(200L))
    expect_equal(.dipper_progress_settings(50, TRUE), on(50L))
    expect_equal(.dipper_progress_settings(TRUE, TRUE), on(200L))
    expect_equal(.dipper_progress_settings(FALSE, TRUE), off)
    expect_equal(.dipper_progress_settings(0, TRUE), off)
    expect_equal(.dipper_progress_settings(-1, TRUE), off)

    # verbose = FALSE overrides print.progress
    expect_equal(.dipper_progress_settings(50, FALSE), off)
    expect_equal(.dipper_progress_settings(TRUE, FALSE), off)
})


test_that(".dipper_check_compiled detects missing executables", {

    # A model whose executable exists passes
    exe <- tempfile()
    file.create(exe)
    on.exit(unlink(exe), add = TRUE)
    compiled <- list(exe_file = function() exe)
    expect_silent(.dipper_check_compiled(compiled))

    # Missing path, empty path, non-existent file, or an error all fail
    msg <- "have not been compiled"
    expect_error(
        .dipper_check_compiled(list(exe_file = function() character(0))),
        msg
    )
    expect_error(
        .dipper_check_compiled(list(exe_file = function() "")),
        msg
    )
    expect_error(
        .dipper_check_compiled(
            list(exe_file = function() file.path(tempdir(), "no_such_model"))
        ),
        msg
    )
    expect_error(
        .dipper_check_compiled(list(exe_file = function() stop("boom"))),
        msg
    )
})


test_that(".dipper_diagnostic_vars selects the checked parameters", {

    # Full diagnostics check all parameters
    expect_null(.dipper_diagnostic_vars("full", symmetric = FALSE, P = 3))

    # Basic diagnostics: nu only for asymmetric, beta_cov only with covariates
    expect_equal(
        .dipper_diagnostic_vars("basic", symmetric = FALSE, P = 2),
        c("alpha", "beta", "tau", "nu", "beta_cov")
    )
    expect_equal(
        .dipper_diagnostic_vars("basic", symmetric = TRUE, P = 1),
        c("alpha", "beta", "tau")
    )
})


test_that(".dipper_convergence_warnings reports each problem", {

    good <- list(num_divergent = 0L, max_rhat = 1.001,
                 min_ess_bulk = 1000, min_ess_tail = 900)
    expect_identical(.dipper_convergence_warnings(good, 2000), character(0))

    bad <- list(num_divergent = 7L, max_rhat = 1.2,
                min_ess_bulk = 50, min_ess_tail = 60)
    msgs <- .dipper_convergence_warnings(bad, niter = 1000)

    expect_length(msgs, 4)
    expect_match(msgs[1], "7 divergent transitions")
    expect_match(msgs[2], "Max R-hat is 1.200. Try niter = 2000", fixed = TRUE)
    expect_match(msgs[3], "Min bulk ESS is 50.0", fixed = TRUE)
    expect_match(msgs[4], "Min tail ESS is 60.0", fixed = TRUE)

    # Thresholds: R-hat 1.01 is already flagged, ESS 400 is not
    edge <- list(num_divergent = 0L, max_rhat = 1.01,
                 min_ess_bulk = 400, min_ess_tail = 400)
    expect_length(.dipper_convergence_warnings(edge, 2000), 1)

    # Missing diagnostics (e.g. all R-hats NA) are not reported
    missing <- list(num_divergent = 0L, max_rhat = NA_real_,
                    min_ess_bulk = NA_real_, min_ess_tail = NA_real_)
    expect_identical(.dipper_convergence_warnings(missing, 2000), character(0))
})


test_that(".dipper_max and .dipper_min handle all-NA input", {
    expect_equal(.dipper_max(c(1, NA, 3)), 3)
    expect_equal(.dipper_min(c(1, NA, 3)), 1)
    expect_identical(.dipper_max(c(NA, NA)), NA_real_)
    expect_identical(.dipper_min(c(NA, NA)), NA_real_)
})
