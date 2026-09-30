library(DiPPER)

# DiPPER fit for the example in the dipper() function documentation
data("tse_hintikka")

fit_example <- dipper(
    tse = tse_hintikka,
    assay.type = "counts",
    data.type = "counts",
    formula = ~ Fat + XOS,
    read.depth = TRUE,
    niter = 400,
    chains = 2,
    cores = 2
)

usethis::use_data(fit_example, overwrite = TRUE, compress = "xz")
