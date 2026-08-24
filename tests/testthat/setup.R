has_cmdstan <- function() {
    has_pkg <- requireNamespace("cmdstanr", quietly = TRUE)
    if (!has_pkg) return(FALSE)

    # Check if cmdstan path is properly set
    path <- tryCatch(
        cmdstanr::cmdstan_path(),
        error = function(e) NULL
    )

    !is.null(path)
}
