# Load the original MultiAssayExperiment object
data("HintikkaXOData", package = "mia")

# Extract the microbiota TreeSummarizedExperiment
tse <- HintikkaXOData[["microbiota"]]

# Extract colData and subset only the Fat and XOS columns
new_coldata <- SummarizedExperiment::colData(HintikkaXOData)[, c("Fat", "XOS")]

# Convert XOS to character
new_coldata$XOS <- ifelse(new_coldata$XOS == 1, "Yes", "No")

# Check that rownames of new_coldata match the colnames of tse
stopifnot(identical(rownames(new_coldata), colnames(tse)))

# Replace the colData of tse with the new_coldata
SummarizedExperiment::colData(tse) <- new_coldata

# Agglomerate to Genus level
tse <- mia::agglomerateByRank(tse, rank = "Genus")

tse_hintikka <- tse

# Save the dataset to the data/ directory
usethis::use_data(tse_hintikka, overwrite = TRUE)
