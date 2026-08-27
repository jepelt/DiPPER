library(curatedMetagenomicData)
library(tidyverse)
library(mia)

# 1. Fetch metadata and remove samples with missing predictor values
meta <- sampleMetadata |>
    filter(
        study_name == "VatanenT_2016",
        !is.na(infant_age),
        !is.na(antibiotics_current_use),
        !is.na(gender)
    )

# 2. Filter metadata to keep measurements from two time-points
meta_sub <- meta |>
    filter(
        (infant_age >= 35 & infant_age <= 65) |
            (infant_age >= 200 & infant_age <= 230)
    ) |>
    mutate(
        age_point = case_when(
            infant_age >= 35 & infant_age <= 65 ~ "1.5 months",
            infant_age >= 200 & infant_age <= 230 ~ "7 months"
        ),
        age_point = factor(age_point, levels = c("1.5 months", "7 months"))
    )

# 3. Just check subjects with an observation in 1st, 2nd, or both timepoints
meta_sub |>
    count(subject_id, age_point) |>
    pivot_wider(names_from = age_point, values_from = n, values_fill = 0) |>
    mutate(
        status = case_when(
            `1.5 months` > 0 & `7 months` > 0 ~ "Both timepoints",
            `1.5 months` > 0 ~ "Only 1.5 months",
            `7 months` > 0 ~ "Only 7 months"
        )
    ) |>
    count(status) |>
    print()

# 4. Fetch the relative abundance data for this subset
tse <- returnSamples(
    meta_sub,
    dataType = "relative_abundance",
    counts = FALSE,
    rownames = "short"
)

# 5. Convert percentages (0-100) to proportions (0-1)
SummarizedExperiment::assay(tse, "relative_abundance") <-
    SummarizedExperiment::assay(tse, "relative_abundance") / 100

# 6. Keep only essential columns, rename, and transform age
keep_cols <- c(
    "subject_id", "antibiotics_current_use", "age_point", "gender"
)
new_coldata <- SummarizedExperiment::colData(tse)[, keep_cols]

# Rename columns
colnames(new_coldata)[colnames(new_coldata) == "antibiotics_current_use"] <-
    "antibiotics"

SummarizedExperiment::colData(tse) <- new_coldata

# 7. Convert relevant columns to factors
tse$subject_id <- as.factor(tse$subject_id)
tse$antibiotics <- factor(tse$antibiotics, levels = c("no", "yes"))
tse$gender <- as.factor(tse$gender)

# 8. Agglomerate to Genus level
tse <- mia::agglomerateByRank(tse, rank = "genus")

VatanenT_2016_subset <- tse

# 9. Check the number of infants
stopifnot(length(unique(VatanenT_2016_subset$subject_id)) == 79)

# 10. Save to the data/ directory
usethis::use_data(VatanenT_2016_subset, overwrite = TRUE)
