test_that("prep_dipper_data works with valid input", {

    # Load the built-in dataset
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
})
