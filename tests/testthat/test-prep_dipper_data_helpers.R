# Unit tests for the internal helpers of prep_dipper_data().

test_that(".dipper_resolve_data_type requires a valid data type", {

    # All four supported types are accepted and returned unchanged
    for (dt in c("counts", "relabundance", "abundance", "pa")) {
        expect_equal(.dipper_resolve_data_type(dt), dt)
    }

    # The data type has no default, so NULL is an error that lists the options
    expect_error(.dipper_resolve_data_type(NULL), "Please specify 'data.type'")
    expect_error(.dipper_resolve_data_type(NULL), "relabundance")

    # Transformed data are not a supported type
    expect_error(.dipper_resolve_data_type("clr"), "should be one of")
})


test_that(".dipper_validate_abundance catches invalid data", {

    counts <- matrix(c(0, 1, 2, 3), nrow = 2)

    # Valid data of each type produce neither errors nor warnings
    expect_silent(.dipper_validate_abundance(counts, "counts"))
    expect_silent(
        .dipper_validate_abundance(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2),
                                   "relabundance")
    )
    expect_silent(
        .dipper_validate_abundance(matrix(c(0, 1.5, 2.7, 100), nrow = 2),
                                   "abundance")
    )
    expect_silent(
        .dipper_validate_abundance(matrix(c(0, 1, 1, 0), nrow = 2), "pa")
    )

    # Negative values are rejected for every data type, because
    # presence/absence cannot be recovered from transformed data
    for (dt in c("counts", "relabundance", "abundance", "pa")) {
        expect_error(
            .dipper_validate_abundance(counts - 1, dt),
            "negative values"
        )
    }

    # Count data must be whole numbers
    expect_error(
        .dipper_validate_abundance(counts + 0.5, "counts"),
        "must contain integers"
    )

    # Relative abundances must be proportions summing to at most 1
    expect_error(
        .dipper_validate_abundance(counts, "relabundance"),
        "between 0 and 1"
    )
    expect_error(
        .dipper_validate_abundance(matrix(c(0.7, 0.7, 0.1, 0.1), nrow = 2),
                                   "relabundance"),
        "Column sums"
    )

    # Presence/absence data must be 0s and 1s
    expect_error(.dipper_validate_abundance(counts, "pa"), "0s and 1s")

    # NAs are not allowed in abundance data
    na_mat <- counts
    na_mat[1, 1] <- NA
    expect_error(.dipper_validate_abundance(na_mat, "counts"), "missing")
})


test_that(".dipper_to_presence_absence honours presence threshold", {

    counts <- matrix(c(0, 5, 10, 0), nrow = 2)

    # By default, any non-zero value is considered presence
    expect_equal(
        .dipper_to_presence_absence(counts, "counts", 0),
        matrix(c(0, 1, 1, 0), nrow = 2)
    )

    # The threshold is always on the scale of the data, and only values
    # strictly above it count as present (here count = 5 is absent)
    expect_equal(
        .dipper_to_presence_absence(counts, "counts", 5),
        matrix(c(0, 0, 1, 0), nrow = 2)
    )

    # A threshold between 0 and 1 is also on the count scale, i.e. it is not
    # interpreted as a relative abundance
    expect_equal(
        .dipper_to_presence_absence(counts, "counts", 0.5),
        matrix(c(0, 1, 1, 0), nrow = 2)
    )

    # For relative abundances the threshold is a proportion
    rel <- matrix(c(0, 0.02, 0.30, 0), nrow = 2)
    expect_equal(
        .dipper_to_presence_absence(rel, "relabundance", 0.05),
        matrix(c(0, 0, 1, 0), nrow = 2)
    )

    # Presence/absence data is passed through untouched
    pa <- matrix(c(0, 1, 1, 0), nrow = 2)
    expect_identical(.dipper_to_presence_absence(pa, "pa", 0.5), pa)

    # A function threshold is applied per feature (row)
    fun_res <- .dipper_to_presence_absence(counts, "counts", median)
    expect_true(all(fun_res %in% c(0, 1)))

    # Invalid thresholds are rejected
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

    # If random effect term is absent, is_longitudinal should be FALSE and the
    # id_var should be NULL.
    plain <- .dipper_parse_random_effect(~ y, meta)
    expect_false(plain$is_longitudinal)
    expect_null(plain$id_var)

    # If random effect term is present, is_longitudinal should be TRUE, and
    # id_var should be read from the formula. fixed_formula should contain the
    # non-random effect terms.
    long <- .dipper_parse_random_effect(~ y + (1 | id), meta)
    expect_true(long$is_longitudinal)
    expect_equal(long$id_var, "id")
    expect_equal(all.vars(long$fixed_formula), "y")

    # The format of the random effect term should be (1 | ...)
    expect_error(
        .dipper_parse_random_effect(~ y + (y | id), meta),
        "Use \\(1 \\| id\\)"
    )

    # Check if the id variable is found in the metadata
    expect_error(
        .dipper_parse_random_effect(~ y + (1 | missing_id), meta),
        "not found"
    )
})


test_that(".dipper_resolve_var_of_interest picks and checks the variable", {

    meta <- data.frame(a = 1:4, b = 5:8)

    # Without an explicit choice, the first formula term is used and reported
    expect_message(
        voi <- .dipper_resolve_var_of_interest(~ a + b, meta, NULL,
                                               verbose = TRUE),
        "Using 'a'"
    )
    expect_equal(voi, "a")

    # verbose = FALSE silences the message
    expect_silent(
        .dipper_resolve_var_of_interest(~ a + b, meta, NULL, verbose = FALSE)
    )

    # An explicit choice is honoured
    expect_equal(
        .dipper_resolve_var_of_interest(~ a + b, meta, "b", verbose = FALSE),
        "b"
    )

    # Formula variables must exist in the metadata
    expect_error(
        .dipper_resolve_var_of_interest(~ a + missing_var, meta, NULL,
                                        verbose = FALSE),
        "missing from the metadata"
    )

    # A formula without predictors cannot define a variable of interest
    expect_error(
        .dipper_resolve_var_of_interest(~ 1, meta, NULL, verbose = FALSE),
        "no predictors"
    )
})


# Read depth handling ---------------------------------------------------------

# Minimal inputs for .dipper_add_read_depth()
rd_inputs <- function(depths = seq(10000, 29000, length.out = 10),
                      extra = NULL) {
    N <- length(depths)
    raw <- matrix(0L, nrow = 4, ncol = N,
                  dimnames = list(paste0("t", 1:4), paste0("s", seq_len(N))))
    # Distribute each sample's depth over the four features
    for (j in seq_len(N)) {
        raw[, j] <- as.integer(c(rep(depths[j] %/% 4, 3),
                                 depths[j] - 3 * (depths[j] %/% 4)))
    }
    meta <- data.frame(group = rep(c("a", "b"), length.out = N),
                       row.names = colnames(raw))
    if (!is.null(extra)) {
        meta <- cbind(meta, extra)
    }
    list(raw_mat = raw, meta_df = meta)
}

add_rd <- function(read.depth, data.type = "counts", inputs = rd_inputs(),
                   formula = ~ group) {
    .dipper_add_read_depth(read.depth, data.type, inputs$raw_mat,
                           inputs$meta_df, formula, formula)
}


test_that(".dipper_add_read_depth requires an explicit choice for counts", {

    # For count data the user must choose, because DiPPER cannot tell whether
    # the sample totals reflect the true read depth
    expect_error(add_rd(NULL), "Please specify 'read.depth' for count data")

    # For other data types the default means no read depth adjustment
    expect_null(add_rd(NULL, data.type = "relabundance")$read.depth.var)
    expect_null(add_rd(NULL, data.type = "abundance")$read.depth.var)
    expect_null(add_rd(NULL, data.type = "pa")$read.depth.var)
})


test_that(".dipper_add_read_depth validates the argument type", {

    msg <- "must be TRUE, FALSE or a single metadata column name"
    expect_error(add_rd(NA), msg)
    expect_error(add_rd(1), msg)
    expect_error(add_rd(c(TRUE, FALSE)), msg)
    expect_error(add_rd(c("a", "b")), msg)
})


test_that("read.depth = FALSE leaves the model untouched", {

    out <- add_rd(FALSE)

    expect_null(out$read.depth.var)
    expect_equal(all.vars(out$fixed_formula), "group")
    expect_false("log10_read_depth" %in% colnames(out$meta_df))
})


test_that("read.depth = TRUE derives read depth from the sample totals", {

    inputs <- rd_inputs()
    out <- add_rd(TRUE, inputs = inputs)

    expect_equal(out$read.depth.var, "log10_read_depth")
    expect_equal(out$meta_df$log10_read_depth,
                 log10(colSums(inputs$raw_mat)),
                 ignore_attr = TRUE)

    # The covariate is added to both formulas automatically
    expect_true("log10_read_depth" %in% all.vars(out$formula))
    expect_true("log10_read_depth" %in% all.vars(out$fixed_formula))

    # Sample totals are only meaningful for count data
    expect_error(
        add_rd(TRUE, data.type = "relabundance"),
        "requires data.type"
    )
})


test_that("a metadata column can be used as the read depth", {

    depths <- seq(10000, 29000, length.out = 10)
    inputs <- rd_inputs(depths, extra = data.frame(seq_depth = depths))

    out <- add_rd("seq_depth", inputs = inputs)

    # The column is log10-transformed and enters the model under the same
    # name as the automatically derived covariate
    expect_equal(out$read.depth.var, "log10_read_depth")
    expect_equal(out$meta_df$log10_read_depth, log10(depths))
    expect_true("log10_read_depth" %in% all.vars(out$fixed_formula))

    # Unknown, non-numeric and already-in-formula columns are rejected
    expect_error(add_rd("nope", inputs = inputs), "not found in metadata")
    inputs_chr <- rd_inputs(depths, extra = data.frame(seq_depth = "x"))
    expect_error(add_rd("seq_depth", inputs = inputs_chr), "must be numeric")
    expect_error(
        add_rd("seq_depth", inputs = inputs, formula = ~ group + seq_depth),
        "Do not include the read depth variable"
    )
})


test_that("read depths must be positive and non-missing", {

    depths <- seq(10000, 19000, length.out = 10)

    zero <- depths; zero[3] <- 0
    expect_error(
        add_rd("seq_depth",
               inputs = rd_inputs(depths, data.frame(seq_depth = zero))),
        "must be positive"
    )

    na_depth <- depths; na_depth[3] <- NA
    expect_error(
        add_rd("seq_depth",
               inputs = rd_inputs(depths, data.frame(seq_depth = na_depth))),
        "must be positive"
    )
})


test_that("read depth warnings fire on log-scale and near-constant input", {

    depths <- seq(10000, 29000, length.out = 10)

    # A column that has already been log10-transformed has implausibly small
    # values
    expect_warning(
        add_rd("seq_depth",
               inputs = rd_inputs(depths,
                                  data.frame(seq_depth = log10(depths)))),
        "unusually small"
    )

    # Rarefied data: read depth cannot explain detection
    expect_warning(
        add_rd(TRUE, inputs = rd_inputs(rep(10000, 10))),
        "vary very little"
    )

    # Real sequencing depths trigger neither warning
    expect_silent(add_rd(TRUE, inputs = rd_inputs(depths)))
})


test_that("a raw depth column and read.depth = TRUE agree", {

    inputs <- rd_inputs()
    totals <- colSums(inputs$raw_mat)
    inputs_col <- rd_inputs(extra = data.frame(seq_depth = totals))

    from_totals <- add_rd(TRUE, inputs = inputs)
    from_column <- add_rd("seq_depth", inputs = inputs_col)

    expect_equal(from_totals$meta_df$log10_read_depth,
                 from_column$meta_df$log10_read_depth,
                 ignore_attr = TRUE)
})


# Remaining helpers -----------------------------------------------------------

test_that(".dipper_scale_covariates standardizes numeric covariates", {

    meta <- data.frame(
        num = c(1, 2, 3, 4),
        fac = factor(c("a", "b", "a", "b"))
    )

    out <- .dipper_scale_covariates(meta, c("num", "fac"))

    # Numeric covariates are centered and scaled, and the scale is recorded so
    # that summary() and plot() can back-transform the estimates
    expect_equal(mean(out$meta_df$num), 0)
    expect_equal(sd(out$meta_df$num), 1)
    expect_equal(out$continuous_scales$num, sd(meta$num))

    # Factors are left alone
    expect_identical(out$meta_df$fac, meta$fac)

    # A constant covariate has no scale and is not rescaled
    const <- data.frame(x = rep(2, 4))
    out_const <- .dipper_scale_covariates(const, "x")
    expect_identical(out_const$meta_df$x, const$x)
    expect_length(out_const$continuous_scales, 0)
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
    kept <- .dipper_filter_taxa(pa, min.present = 1, min.absent = 1,
                                verbose = FALSE)
    expect_equal(rownames(kept), "sometimes")

    # Require presence and absence in 25 percent of samples
    kept_prop <- .dipper_filter_taxa(pa, min.present = 0.25,
                                     min.absent = 0.25, verbose = FALSE)
    expect_equal(rownames(kept_prop), "sometimes")

    # No filtering at all should retain all taxa
    expect_equal(nrow(.dipper_filter_taxa(pa, 0, 0, verbose = FALSE)), 3)

    # The number of removed and retained taxa is reported when verbose
    expect_message(
        .dipper_filter_taxa(pa, 1, 1, verbose = TRUE),
        "2 taxa removed, 1 taxa retained"
    )

    # All taxa filtered out should trigger an error
    expect_error(
        .dipper_filter_taxa(pa, min.present = 4, min.absent = 4,
                            verbose = FALSE),
        "All taxa filtered out"
    )
})


test_that(".dipper_build_design orders and centers the columns", {

    meta <- data.frame(
        group = factor(rep(c("Low", "High"), 3), levels = c("Low", "High")),
        age = c(1, 5, 2, 7, 3, 4)
    )

    X <- .dipper_build_design(~ group + age, meta, "group")

    # The intercept is dropped and the variable of interest comes first
    expect_equal(colnames(X), c("groupHigh", "age"))
    expect_true(all(abs(colMeans(X)) < 1e-12))

    # A variable that is not a formula term cannot be the variable of interest
    expect_error(
        .dipper_build_design(~ group, meta, "age"),
        "not found in design matrix"
    )
})


test_that(".dipper_build_design matches terms, not column name prefixes", {

    # 'age' must not also select the dummy columns of 'age_group', which the
    # design matrix names 'age_groupb'
    meta <- data.frame(
        age_group = factor(rep(c("a", "b"), 3)),
        age = c(1, 5, 2, 7, 3, 4)
    )

    X <- .dipper_build_design(~ age_group + age, meta, "age")

    # The column of the variable of interest is first, on its own
    expect_equal(colnames(X), c("age", "age_groupb"))
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

    # Numeric variables are continuous and have no levels
    expect_null(.dipper_var_levels(meta, "numeric"))

    expect_error(.dipper_var_levels(meta, "three"), "3 levels")
    expect_error(.dipper_var_levels(meta, "constant"), "at least two levels")

    # The error explains that the restriction applies only to the variable of
    # interest
    expect_error(.dipper_var_levels(meta, "three"), "control variables")
})
