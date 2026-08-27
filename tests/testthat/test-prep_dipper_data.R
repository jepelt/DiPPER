test_that("prep_dipper_data works with valid input", {

    data("tse_hintikka", package = "DiPPER")

    prep <- prep_dipper_data(
        tse = tse_hintikka,
        formula = ~ Fat + XOS,
        assay.type = "counts"
    )

    expect_type(prep, "list")
    expect_true("y" %in% names(prep))
    expect_true("X" %in% names(prep))

    # Fat + XOS + log10_read.depth = 3 columns in the design matrix
    # (intercept is separate)
    expect_equal(prep$P, 3)

    # The variable of interest is taken from the first formula term and its
    # columns are moved to the front of the design matrix
    expect_equal(prep$var.of.interest, "Fat")
    expect_true(grepl("^Fat", prep$design_matrix_cols[1]))

    # Presence/absence matrix and dimensions are consistent
    expect_true(all(prep$y %in% c(0, 1)))
    expect_equal(prep$K, nrow(prep$y))
    expect_equal(prep$N, ncol(prep$y))
    expect_equal(prep$N, nrow(prep$X))
    expect_equal(length(prep$taxa_names), prep$K)

    # Design matrix columns are centered
    expect_true(all(abs(colMeans(prep$X)) < 1e-8))

    # Cross-sectional data
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

    expect_true(prep$is_longitudinal)
    expect_equal(prep$id_var, "subject_id")
    expect_equal(prep$var.of.interest, "age_point")

    # One random intercept per subject, indices within range
    n_subjects <- length(unique(VatanenT_2016_subset$subject_id))
    expect_equal(prep$S, n_subjects)
    expect_equal(length(prep$subj), prep$N)
    expect_true(all(prep$subj >= 1 & prep$subj <= prep$S))

    # The random effect term is dropped from the fixed formula
    expect_false(grepl("\\|", paste(deparse(prep$fixed_formula),
                                    collapse = " ")))

    # No read depth adjustment was requested
    expect_null(prep$read.depth.var)
})


test_that("prep_dipper_data validates its input", {

    data("tse_hintikka", package = "DiPPER")

    expect_error(
        prep_dipper_data(tse = tse_hintikka, formula = ~ Fat),
        "assay.type"
    )
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            formula = ~ Fat,
            assay.type = "not_an_assay"
        ),
        "not found"
    )
    expect_error(
        prep_dipper_data(
            tse = tse_hintikka,
            formula = ~ NotAVariable,
            assay.type = "counts"
        ),
        "missing from the metadata"
    )
    expect_error(
        prep_dipper_data(formula = ~ Fat),
        "Provide either"
    )
})


test_that("prep_dipper_data rejects unsupported random effects", {

    data("VatanenT_2016_subset", package = "DiPPER")

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
