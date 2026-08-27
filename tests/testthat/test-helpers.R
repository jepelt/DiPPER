# Unit tests for the internal helpers of prep_dipper_data(). These need
# neither CmdStan nor the bundled datasets, so they run everywhere and are
# fast.

test_that(".dipper_validate_abundance catches invalid data", {

    counts <- matrix(c(0, 1, 2, 3), nrow = 2)

    expect_silent(.dipper_validate_abundance(counts, "counts"))

    na_mat <- counts
    na_mat[1, 1] <- NA
    expect_error(.dipper_validate_abundance(na_mat, "counts"), "missing")

    expect_error(
        .dipper_validate_abundance(counts - 1, "counts"),
        "non-negative integers"
    )
    expect_error(
        .dipper_validate_abundance(counts + 0.5, "counts"),
        "non-negative integers"
    )
    expect_error(
        .dipper_validate_abundance(counts, "relabundance"),
        "between 0 and 1"
    )
    expect_error(
        .dipper_validate_abundance(counts, "pa"),
        "0s and 1s"
    )

    rel <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
    expect_silent(.dipper_validate_abundance(rel, "relabundance"))

    pa <- matrix(c(0, 1, 1, 0), nrow = 2)
    expect_silent(.dipper_validate_abundance(pa, "pa"))
})


test_that(".dipper_to_presence_absence honours the threshold", {

    counts <- matrix(c(0, 5, 10, 0), nrow = 2)

    # Default: anything above zero is present
    expect_equal(
        .dipper_to_presence_absence(counts, "counts", 0),
        matrix(c(0, 1, 1, 0), nrow = 2)
    )

    # Whole-number threshold applies to the counts themselves
    expect_equal(
        .dipper_to_presence_absence(counts, "counts", 5),
        matrix(c(0, 0, 1, 0), nrow = 2)
    )

    # Presence/absence data is passed through untouched
    pa <- matrix(c(0, 1, 1, 0), nrow = 2)
    expect_identical(.dipper_to_presence_absence(pa, "pa", 0.5), pa)

    # A function threshold is applied per feature (row)
    fun_res <- .dipper_to_presence_absence(counts, "counts", median)
    expect_true(all(fun_res %in% c(0, 1)))

    expect_error(
        .dipper_to_presence_absence(counts, "counts", -1),
        "cannot be negative"
    )
    expect_error(
        .dipper_to_presence_absence(counts, "counts", "high"),
        "numeric value or a function"
    )
})


test_that(".dipper_check_alignment compares samples", {

    pa <- matrix(0, nrow = 2, ncol = 3,
                 dimnames = list(c("t1", "t2"), c("s1", "s2", "s3")))
    meta <- data.frame(x = 1:3, row.names = c("s1", "s2", "s3"))

    expect_silent(.dipper_check_alignment(pa, meta))

    expect_error(
        .dipper_check_alignment(pa, meta[1:2, , drop = FALSE]),
        "must match rows in metadata"
    )

    wrong <- meta
    rownames(wrong) <- c("s1", "s3", "s2")
    expect_error(.dipper_check_alignment(pa, wrong), "must match rownames")

    # Without names, a warning is issued but processing continues
    pa_nameless <- pa
    dimnames(pa_nameless) <- NULL
    expect_warning(
        .dipper_check_alignment(pa_nameless, meta),
        "identically ordered"
    )
})


test_that(".dipper_parse_random_effect splits the formula", {

    meta <- data.frame(y = 1:4, id = factor(c("a", "a", "b", "b")))

    plain <- .dipper_parse_random_effect(~ y, meta)
    expect_false(plain$is_longitudinal)
    expect_null(plain$id_var)

    long <- .dipper_parse_random_effect(~ y + (1 | id), meta)
    expect_true(long$is_longitudinal)
    expect_equal(long$id_var, "id")
    expect_equal(all.vars(long$fixed_formula), "y")

    expect_error(
        .dipper_parse_random_effect(~ y + (y | id), meta),
        "Use \\(1 \\| id\\)"
    )
    expect_error(
        .dipper_parse_random_effect(~ y + (1 | missing_id), meta),
        "not found"
    )
})


test_that(".dipper_filter_taxa filters by prevalence", {

    # Three taxa: present in 4, 2 and 0 of 4 samples
    pa <- rbind(
        c(1, 1, 1, 1),
        c(1, 1, 0, 0),
        c(0, 0, 0, 0)
    )
    rownames(pa) <- c("always", "sometimes", "never")

    # Require presence in >= 1 and absence in >= 1 sample
    kept <- .dipper_filter_taxa(pa, min.present = 1, min.absent = 1)
    expect_equal(rownames(kept), "sometimes")

    # Proportions are converted to sample counts
    kept_prop <- .dipper_filter_taxa(pa, min.present = 0.25,
                                     min.absent = 0.25)
    expect_equal(rownames(kept_prop), "sometimes")

    # No filtering at all
    expect_equal(nrow(.dipper_filter_taxa(pa, 0, 0)), 3)

    expect_error(
        .dipper_filter_taxa(pa, min.present = 4, min.absent = 4),
        "All taxa filtered out"
    )
})


test_that(".dipper_var_levels validates the variable of interest", {

    meta <- data.frame(
        binary = factor(c("Low", "High", "Low", "High"),
                        levels = c("Low", "High")),
        three = factor(c("a", "b", "c", "a")),
        constant = factor(rep("only", 4)),
        numeric = c(1.5, 2.5, 3.5, 4.5)
    )

    expect_equal(.dipper_var_levels(meta, "binary"), c("Low", "High"))
    expect_null(.dipper_var_levels(meta, "numeric"))

    expect_error(.dipper_var_levels(meta, "three"), "3 levels")
    expect_error(.dipper_var_levels(meta, "constant"), "at least two levels")
})


test_that(".dipper_scale_covariates standardizes numeric covariates", {

    meta <- data.frame(
        num = c(1, 2, 3, 4),
        depth = c(10, 20, 30, 40),
        fac = factor(c("a", "b", "a", "b"))
    )

    out <- .dipper_scale_covariates(meta, c("num", "depth", "fac"), "depth")

    # Ordinary numeric covariates are centered and scaled
    expect_equal(mean(out$meta_df$num), 0)
    expect_equal(sd(out$meta_df$num), 1)
    expect_equal(out$continuous_scales$num, sd(meta$num))

    # The read depth covariate is centered but not scaled
    expect_equal(mean(out$meta_df$depth), 0)
    expect_equal(sd(out$meta_df$depth), sd(meta$depth))
    expect_null(out$continuous_scales$depth)

    # Factors are left alone
    expect_identical(out$meta_df$fac, meta$fac)
})
