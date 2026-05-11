#' Run the DiPPER model
#'
#' This is a wrapper function that sequentially runs \code{prep_dipper_data}
#' and \code{run_dipper}.
#'
#' @inheritParams prep_dipper_data
#' @inheritParams run_dipper
#' @param ... Additional arguments passed directly to \code{run_dipper}.
#'
#' @return The output of the \code{run_dipper} function.
#' @export
DiPPER <- function(tse = NULL,
                   meta = NULL,
                   abund_matrix = NULL,
                   formula,
                   var_of_interest = NULL,
                   read_depth = TRUE,
                   symmetric = F,
                   assay_type = c("counts", "relabundance", "pa"),
                   assay_name = NULL,
                   tax_rank = NULL,
                   pa_threshold = 0,
                   min_present = 5,
                   min_absent = min_present,
                   iter_sampling = 1000,
                   cores = 4,
                   seed = 1,
                   ...) {

    # Prepare the data with standard settings for covariates
    prepared_data <- prep_dipper_data(
        tse = tse,
        meta = meta,
        abund_matrix = abund_matrix,
        formula = formula,
        var_of_interest = var_of_interest,
        read_depth = read_depth,
        assay_type = assay_type,
        assay_name = assay_name,
        tax_rank = tax_rank,
        pa_threshold = pa_threshold,
        min_present = min_present,
        min_absent = min_absent,
        standardize = TRUE,
        center_dummies = TRUE
    )

    # Run the core model
    result <- run_dipper(
        prep_data = prepared_data,
        symmetric = symmetric,
        iter_sampling = iter_sampling,
        cores = cores,
        seed = 1,
        ...
    )

    return(result)
}
