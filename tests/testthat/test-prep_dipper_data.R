# Tests for prep_dipper_data() itself. None of these require CmdStan.

test_that("prep_dipper_data works with valid input", {

    data("tse_hintikka", package = "DiPPER")

    prep <- prep_dipper_data(
        tse = tse_hintikka,
        assay.type = "counts",
        data.type = "counts",
        formula = ~ Fat + XOS,
        read.depth = TRUE,
        verbose = FALSE
    )

    # prep_dipper_data() returns a list containing (among other objects)
    # the presence/absence matrix y and the design matrix X
    expect_type(prep, "list")
    expect_true("y" %in% names(prep))
    expect_true("X" %in% names(prep))

    # The formula ~ Fat + XOS plus automatically added log10_read_depth should
    # lead to a design matrix with 3 columns (intercept is separate).
    expect_equal(prep$P, 3)
    expect_equal(prep$read.depth.var, "log10_read_depth")
    expect_true("log10_read_depth" %in% prep$design_matrix_cols)

    # The variable of interest is the first formula term and its column is
    # moved to the front of the design matrix
    expect_equal(prep$var.of.interest, "Fat")
    expect_true(startsWith(prep$design_matrix_cols[1], "Fat"))

    # Presence/absence matrix and dimensions (number of samples N and number of
    # taxa/features K) are consistent
    expect_true(all(prep$y %in% c(0, 1)))
    expect_equal(prep$K, nrow(prep$y))
    expect_equal(prep$N, ncol(prep$y))
    expect_equal(prep$N, nrow(prep$X))
    expect_equal(length(prep$taxa_names), prep$K)

    # The number of features before prevalence filtering is recorded, so that
    # print() can report how many were removed
    expect_true(prep$K_unfiltered >= prep$K)
    expect_equal(prep$K_unfiltered, nrow(SummarizedExperiment::assay(
        tse_hintikka, "counts")))

    # Design matrix columns are centered
    expect_true(all(abs(colMeans(prep$X)) < 1e-8))

    # The example does not include random intercept -> not longitudinal.
    expect_false(prep$is_longitudinal)
    expect_null(prep$id_var)
    expect_equal(prep$S, 0)
})


test_that("prep_dipper_data handles longitudinal formulas", {

    data("VatanenT_2016_subset", package = "DiPPER")

    prep <- prep_dipper_data(
        tse = VatanenT_2016_subset,
        assay.type = "relative_abundance",
        data.type = "relabundance",
        formula = ~ age_point + antibiotics + gender + (1 | subject_id),
        read.depth = FALSE,
        verbose = FALSE
    )

    # Longitudinal example with random intercepts should set is_longitudinal =
    # TRUE and id_var = "subject_id" (from (1 | subject_id) in the formula).
    expect_true(prep$is_longitudinal)
    expect_equal(prep$id_var, "subject_id")

    # Variable of interest is the first formula term in the formula
    expect_equal(prep$var.of.interest, "age_point")

    # Number of subjects in the original TSE object matches the number of
    # subjects in the output (S).
    n_subjects <- length(unique(VatanenT_2016_subset$subject_id))
    expect_equal(prep$S, n_subjects)

    # The length of the subject ID vector (prep$subj) should match the number
    # of samples (N).
    expect_equal(length(prep$subj), prep$N)

    # The subject IDs in prep$subj should be integers between 1 and S.
    expect_true(all(prep$subj >= 1 & prep$subj <= prep$S))

    # The random effect term is dropped from the fixed formula
    expect_false(grepl("\\|", paste(deparse(prep$fixed_formula),
                                    collapse = " ")))

    # No read depth adjustment was requested
    expect_null(prep$read.depth.var)
})


test_that("prep_dipper_data requires formula, data.type and read.depth", {

    data("tse_hintikka", package = "DiPPER")

    base_args <- list(tse = tse_hintikka, assay.type = "counts",
                      verbose = FALSE)

    # The model formula has no default
    expect_error(
        do.call(prep_dipper_data, c(base_args, list(data.type = "counts",
                                                    read.depth = TRUE))),
        "Please specify 'formula'"
    )

    # The data type has no default, because it determines how presence is
    # defined and whether read depth can be derived from the data
    expect_error(
        do.call(prep_dipper_data, c(base_args, list(formula = ~ Fat))),
        "Please specify 'data.type'"
    )
    expect_error(
        do.call(prep_dipper_data, c(base_args,
                                    list(formula = ~ Fat,
                                         data.type = "clr"))),
        "should be one of"
    )

    # For count data the read depth handling must be chosen explicitly
    expect_error(
        do.call(prep_dipper_data, c(base_args,
                                    list(formula = ~ Fat,
                                         data.type = "counts"))),
        "Please specify 'read.depth' for count data"
    )
})


test_that("prep_dipper_data validates its input", {

    data("tse_hintikka", package = "DiPPER")

    # Assay name (assay.type) is missing
    expect_error(
        prep_dipper_data(tse = tse_hintikka, formula = ~ Fat,
                         data.type = "counts", read.depth = TRUE,
                         verbose = FALSE),
        "assay.type"
    )

    # assay.type that is not contained by the tse object is given
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            assay.type = "not_an_assay",
            data.type = "counts",
            formula = ~ Fat,
            read.depth = TRUE,
            verbose = FALSE
        ),
        "not found"
    )

    # Variable in the formula is not contained in (colData of) the tse object
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            assay.type = "counts",
            data.type = "counts",
            formula = ~ NotAVariable,
            read.depth = TRUE,
            verbose = FALSE
        ),
        "missing from the metadata"
    )

    # TreeSummarizedExperiment object or metadata and abundance matrix must be
    # provided
    expect_error(
        prep_dipper_data(formula = ~ Fat, data.type = "counts",
                         read.depth = TRUE, verbose = FALSE),
        "Provide either"
    )
})


test_that("prep_dipper_data rejects unsupported random effects", {

    data("VatanenT_2016_subset", package = "DiPPER")

    long_args <- list(
        tse = VatanenT_2016_subset,
        assay.type = "relative_abundance",
        data.type = "relabundance",
        read.depth = FALSE,
        verbose = FALSE
    )

    # DiPPER currently only supports random intercepts (1 | ...).
    expect_error(
        do.call(prep_dipper_data, c(
            long_args, list(formula = ~ age_point + (age_point | subject_id))
        )),
        "Use \\(1 \\| id\\)"
    )

    # ID variable given that is not found in the metadata
    expect_error(
        do.call(prep_dipper_data, c(
            long_args, list(formula = ~ age_point + (1 | not_a_column))
        )),
        "not found"
    )
})


test_that("read.depth = TRUE requires count data", {

    data("VatanenT_2016_subset", package = "DiPPER")

    expect_error(
        prep_dipper_data(
            tse = VatanenT_2016_subset,
            assay.type = "relative_abundance",
            data.type = "relabundance",
            formula = ~ age_point,
            read.depth = TRUE,
            verbose = FALSE
        ),
        "requires data.type"
    )
})


# The tests below use small synthetic matrices (see helper-synthetic.R), which
# can represent inputs that the packaged example data cannot.

test_that("matrix and metadata input is accepted for every data type", {

    counts <- make_counts()
    meta <- make_meta()

    # Counts, with and without read depth control
    with_rd <- prep_synthetic(counts, meta, data.type = "counts",
                              read.depth = TRUE)
    expect_equal(with_rd$design_matrix_cols,
                 c("groupHigh", "log10_read_depth"))

    without_rd <- prep_synthetic(counts, meta, data.type = "counts",
                                 read.depth = FALSE)
    expect_equal(without_rd$design_matrix_cols, "groupHigh")
    expect_null(without_rd$read.depth.var)

    # Other non-negative abundances, e.g. pathway abundances or CPM values.
    # Read depth cannot be derived from them, so it is simply not controlled
    # for unless a metadata column is supplied.
    abundance <- counts * 1.7
    ab <- prep_synthetic(abundance, meta, data.type = "abundance")
    expect_null(ab$read.depth.var)
    expect_equal(ab$K, with_rd$K)

    # Presence/absence input is used as is
    pa <- prep_synthetic(ifelse(counts > 0, 1, 0), meta, data.type = "pa")
    expect_true(all(pa$y %in% c(0, 1)))
    expect_equal(pa$K, with_rd$K)

    # Transformed data with negative values are rejected
    expect_error(
        prep_synthetic(log(abundance + 0.5), meta, data.type = "abundance"),
        "negative values"
    )
})


test_that("features without rownames are given default names", {

    counts <- make_counts(named = FALSE)
    meta <- make_meta()

    prep <- prep_synthetic(counts, meta, data.type = "counts",
                           read.depth = FALSE)

    expect_equal(prep$taxa_names, paste0("feature", seq_len(prep$K)))

    # summary() can label its rows, instead of failing on a NULL name vector
    expect_equal(length(prep$taxa_names), prep$K)

    # The user is told where the names come from
    # (suppressMessages() hides the unrelated variable-of-interest message)
    suppressMessages(expect_message(
        prep_dipper_data(assay = counts, meta = meta, formula = ~ group,
                         data.type = "counts", read.depth = FALSE,
                         min.present = 2, min.absent = 2),
        "no rownames"
    ))
})


test_that("prep_dipper_data warns about uninformative read depths", {

    meta <- make_meta()

    # Rarefied counts: every sample has the same total, so read depth cannot
    # explain differences in detection
    counts <- make_counts()
    rarefied <- round(sweep(counts, 2, colSums(counts), "/") * 5000)

    expect_warning(
        prep_synthetic(rarefied, meta, data.type = "counts",
                       read.depth = TRUE),
        "vary very little"
    )

    # A metadata column that has already been log10-transformed
    meta_log <- meta
    meta_log$depth <- log10(colSums(counts))
    expect_warning(
        prep_synthetic(counts, meta_log, data.type = "counts",
                       read.depth = "depth"),
        "unusually small"
    )

    # Ordinary sequencing counts trigger neither warning
    expect_no_warning(
        prep_synthetic(counts, meta, data.type = "counts", read.depth = TRUE)
    )
})


test_that("a raw read depth column and read.depth = TRUE give the same model", {

    counts <- make_counts()
    meta <- make_meta()
    meta$seq_depth <- colSums(counts)

    from_totals <- prep_synthetic(counts, meta, data.type = "counts",
                                  read.depth = TRUE)
    from_column <- prep_synthetic(counts, meta, data.type = "counts",
                                  read.depth = "seq_depth")

    expect_equal(from_totals$design_matrix_cols,
                 from_column$design_matrix_cols)
    expect_equal(from_totals$X, from_column$X)
})


test_that("continuous covariates are scaled and the scales recorded", {

    counts <- make_counts()
    meta <- make_meta()

    prep <- prep_synthetic(counts, meta, formula = ~ age + group,
                           data.type = "counts", read.depth = TRUE)

    # The variable of interest comes first, then the other covariates
    expect_equal(prep$var.of.interest, "age")
    expect_equal(prep$design_matrix_cols,
                 c("age", "groupHigh", "log10_read_depth"))

    # The scale of the continuous variable of interest is stored, so that
    # summary() and plot() can report results on the original scale
    expect_equal(prep$continuous_scales$age, sd(meta$age))

    # Read depth is centered but not scaled, so its coefficient stays on the
    # log10 scale
    expect_equal(sd(prep$X[, "log10_read_depth"]),
                 sd(log10(colSums(counts))))
})


test_that("prep_dipper_data rejects metadata with missing values", {

    counts <- make_counts()
    meta <- make_meta()
    meta$age[2] <- NA

    expect_error(
        prep_synthetic(counts, meta, formula = ~ group + age,
                       data.type = "counts", read.depth = FALSE),
        "Metadata contains NAs"
    )
})
