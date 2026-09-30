# Create small synthetic datasets.

# Sparse count matrix: every feature is present in some but not all samples,
# so that nothing is removed by the default prevalence filtering below.
make_counts <- function(K = 12, N = 20, seed = 1, named = TRUE) {
    set.seed(seed)
    m <- matrix(0L, nrow = K, ncol = N)
    for (i in seq_len(K)) {
        present <- sample(N, sample(seq(3, N - 3), 1))
        m[i, present] <- rpois(length(present), 40) + 1L
    }
    colnames(m) <- paste0("s", seq_len(N))
    if (named) {
        rownames(m) <- paste0("taxon", seq_len(K))
    }
    m
}


# Simple metadata
make_meta <- function(N = 20, seed = 1) {
    set.seed(seed + 100)
    data.frame(
        group = factor(rep(c("Low", "High"), length.out = N),
                       levels = c("Low", "High")),
        age = rnorm(N),
        row.names = paste0("s", seq_len(N))
    )
}


# prep_dipper_data() call with the synthetic data and light filtering
prep_synthetic <- function(counts = make_counts(), meta = make_meta(),
                           formula = ~ group, ...) {
    prep_dipper_data(
        assay = counts, meta = meta, formula = formula,
        min.present = 2, min.absent = 2, verbose = FALSE, ...
    )
}
