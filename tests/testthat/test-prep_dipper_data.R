test_that("prep_dipper_data works with valid input", {

    data("tse_hintikka", package = "DiPPER")

    prep <- prep_dipper_data(
        tse = tse_hintikka,
        formula = ~ Fat + XOS,
        assay.type = "counts"
    )

    # prep_dipper_data() returns a list containing (among other objects)
    # the presence/absence matrix y and the design matrix X
    expect_type(prep, "list")
    expect_true("y" %in% names(prep))
    expect_true("X" %in% names(prep))

    # The formula ~ Fat + XOS plus automatically added log10_read.depth should
    # lead to a design matrix with 3 columns (intercept is separate).
    expect_equal(prep$P, 3)

    # The variable of interest is the first formula term and its column is
    # moved to the front of the design matrix
    expect_equal(prep$var.of.interest, "Fat")
    expect_true(grepl("^Fat", prep$design_matrix_cols[1]))

    # Presence/absence matrix and dimensions (number of samples N and number of
    # taxa/features K) are consistent
    expect_true(all(prep$y %in% c(0, 1)))
    expect_equal(prep$K, nrow(prep$y))
    expect_equal(prep$N, ncol(prep$y))
    expect_equal(prep$N, nrow(prep$X))
    expect_equal(length(prep$taxa_names), prep$K)

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
        formula = ~ age_point + antibiotics + gender + (1 | subject_id),
        assay.type = "relative_abundance",
        data.type = "relabundance",
        read.depth = FALSE
    )

    # Longitudinal example with random intercepts should set is_longitudinal =
    # TRUE and id_var = "subject_id (from (1 | subject_id) in the formula).
    expect_true(prep$is_longitudinal)
    expect_equal(prep$id_var, "subject_id")

    # Variable of interest is the first formula term in the formula
    expect_equal(prep$var.of.interest, "age_point")

    # Number of subjects in the original TSE object matches the number of
    # subjects in the output (S).
    n_subjects <- length(unique(VatanenT_2016_subset$subject_id))
    expect_equal(prep$S, n_subjects)

    # The length of the subject ID vector (prep$subj) should match the number
    # of samples (N), i.e. the number columns in the presence/absence
    # matrix. Note that N > S here because there are multiple samples per
    # subject.
    expect_equal(length(prep$subj), prep$N)

    # The subject IDs in prep$subj should be integers between 1 and S.
    expect_true(all(prep$subj >= 1 & prep$subj <= prep$S))

    # The random effect term is dropped from the fixed formula
    expect_false(grepl("\\|", paste(deparse(prep$fixed_formula),
                                    collapse = " ")))

    # No read depth adjustment was requested
    expect_null(prep$read.depth.var)
})


test_that("prep_dipper_data validates its input", {

    data("tse_hintikka", package = "DiPPER")

    # Assay name (assay.type) is missing
    expect_error(
        prep_dipper_data(tse = tse_hintikka, formula = ~ Fat),
        "assay.type"
    )

    # assay.type that is not contained by the tse object is given
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            formula = ~ Fat,
            assay.type = "not_an_assay"
        ),
        "not found"
    )

    # Variable in the formula is not contained in (colData of) the tse object
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            formula = ~ NotAVariable,
            assay.type = "counts"
        ),
        "missing from the metadata"
    )

    # TreeSummarizedExperiment object or metadata and abundance matrix must be
    # provided
    expect_error(
        prep_dipper_data(formula = ~ Fat),
        "Provide either"
    )
})


test_that("prep_dipper_data rejects unsupported random effects", {

    data("VatanenT_2016_subset", package = "DiPPER")

    # DiPPER currently only supports random intercepts (1 | ...).
    expect_error(
        prep_dipper_data(
            tse = VatanenT_2016_subset,
            formula = ~ age_point + (age_point | subject_id),
            assay.type = "relative_abundance",
            data.type = "relabundance",
            read.depth = FALSE
        ),
        "Use \\(1 \\| id\\)"
    )

    # ID variable given that is not found in the metadata
    expect_error(
        prep_dipper_data(
            tse = VatanenT_2016_subset,
            formula = ~ age_point + (1 | not_a_column),
            assay.type = "relative_abundance",
            data.type = "relabundance",
            read.depth = FALSE
        ),
        "not found"
    )
})


test_that("read.depth = TRUE requires count data", {

    data("VatanenT_2016_subset", package = "DiPPER")

    expect_error(
        prep_dipper_data(
            tse = VatanenT_2016_subset,
            formula = ~ age_point,
            assay.type = "relative_abundance",
            data.type = "relabundance",
            read.depth = TRUE
        ),
        "requires data.type"
    )
})
