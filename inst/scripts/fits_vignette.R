library(DiPPER)

# DiPPER fit for the Quick start example in the vignette
data("tse_hintikka")

fit_hintikka <- dipper(
    tse = tse_hintikka,
    assay.type = "counts",
    data.type = "counts",
    formula = ~ Fat + XOS,
    read.depth = TRUE
)

saveRDS(fit_hintikka, file = "inst/extdata/fit_hintikka.rds", compress = "xz")


# DiPPER fit for the longitudinal data example in the vignette
data("VatanenT_2016_subset")

fit_vatanen <- dipper(
    tse = VatanenT_2016_subset,
    assay.type = "relative_abundance",
    data.type = "relabundance",
    formula = ~ age_point + antibiotics + gender + (1 | subject_id),
    read.depth = FALSE
)

saveRDS(fit_vatanen, file = "inst/extdata/fit_vatanen.rds", compress = "xz")
