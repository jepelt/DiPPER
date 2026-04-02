#' Prepare data for the DiPPER model
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param meta A data.frame containing metadata (required if tse is NULL).
#' @param abund_matrix A matrix containing counts, relative abundances, or
#'   presence/absence data (required if tse is NULL).
#' @param formula Model formula.
#' @param var_of_interest The variable of interest from the formula.
#'   If NULL, the first term of the formula is automatically used.
#' @param read_depth Logical or character. If TRUE, calculates log10 read depth
#'   from counts. If character, uses the specified metadata column.
#' @param assay_type Character. Type of the assay: "counts" (default),
#'   "relabundance", or "pa" (presence/absence).
#' @param assay_name Name of the assay in the TSE object to use. If NULL
#'   (default), it is automatically determined based on assay_type. Ignored if
#'   tse is NULL.
#' @param tax_rank A character string specifying the taxonomic rank to
#'   agglomerate by (e.g., "genus"). If NULL, no agglomeration is done.
#'   Ignored if tse is NULL.
#' @param pa_threshold Numeric or function. Threshold for presence. If strictly
#'   between 0 and 1, relative abundances are used. If >= 1 or exactly 0,
#'   counts are used. If a function (e.g., mean), it is applied per taxon.
#'   Default is 0. Ignored if assay_type is "pa".
#' @param min_present Minimum number or proportion of samples where a taxon
#'   must be present.
#' @param min_absent Minimum number or proportion of samples where a taxon
#'   must be absent. Default is min_present.
#' @param standardize Logical. Standardize numeric covariates (mean=0, sd=1).
#' @param center_dummies Logical. Mean-center all columns in design matrix.
#'
#' @export
#' @importFrom stats model.matrix update complete.cases terms
prep_dipper_data <- function(tse = NULL,
                             meta = NULL,
                             abund_matrix = NULL,
                             formula,
                             var_of_interest = NULL,
                             read_depth = TRUE,
                             assay_type = c("counts", "relabundance", "pa"),
                             assay_name = NULL,
                             tax_rank = NULL,
                             pa_threshold = 0,
                             min_present = 5,
                             min_absent = min_present,
                             standardize = TRUE,
                             center_dummies = TRUE) {

    # Set default assay_name based on matched assay_type if NULL
    if (is.null(assay_name)) {
        assay_type <- match.arg(assay_type)
        assay_name <- assay_type
    }

    assay_type <- match.arg(assay_type)
    raw_mat <- NULL
    meta_df <- NULL

    # 1. Input Extraction
    if (!is.null(tse)) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("Package 'SummarizedExperiment' is required for TSE objects.")
        }

        if (!is.null(tax_rank)) {
            if (!requireNamespace("mia", quietly = TRUE)) {
                stop(paste0(
                    "Package 'mia' is required for taxonomic agglomeration. ",
                    "Please install it or set tax_rank = NULL."
                ))
            }
            available_ranks <- mia::taxonomyRanks(tse)
            if (!tax_rank %in% available_ranks) {
                stop(sprintf(
                    "Taxonomic rank '%s' not found. Available ranks: %s",
                    tax_rank, paste(available_ranks, collapse = ", ")
                ))
            }
            tse <- mia::agglomerateByRank(tse, rank = tax_rank)
        }

        raw_mat <- SummarizedExperiment::assay(tse, assay_name)
        meta_df <- as.data.frame(SummarizedExperiment::colData(tse))

    } else if (!is.null(meta) && !is.null(abund_matrix)) {
        raw_mat <- as.matrix(abund_matrix)
        meta_df <- as.data.frame(meta)
    } else {
        stop("Provide either 'tse', OR both 'meta' and 'abund_matrix'.")
    }

    # 2. Validation and Presence/Absence Conversion
    if (anyNA(raw_mat)) {
        stop("The abundance data contains missing values (NA).")
    }

    if (assay_type == "counts") {
        if (any(raw_mat < 0) || any(raw_mat %% 1 != 0)) {
            stop("Count data must contain non-negative integers.")
        }
    } else if (assay_type == "relabundance") {
        if (any(raw_mat < 0) || any(raw_mat > 1)) {
            stop("Relative abundance values must be between 0 and 1.")
        }
        if (any(colSums(raw_mat) > 1.0001)) {
            stop("Column sums of relative abundances must be <= 1.")
        }
    } else if (assay_type == "pa") {
        if (!all(raw_mat %in% c(0, 1))) {
            stop("Presence/absence data must strictly contain 0s and 1s.")
        }
    }

    if (assay_type == "pa") {
        pa_matrix <- raw_mat
    } else {
        if (is.function(pa_threshold)) {
            row_thresh <- apply(raw_mat, 1, pa_threshold)
            pa_matrix <- ifelse(raw_mat > row_thresh, 1, 0)

        } else if (is.numeric(pa_threshold)) {
            if (pa_threshold < 0) {
                stop("pa_threshold cannot be negative.")
            } else if (pa_threshold > 0 && pa_threshold < 1) {
                if (assay_type == "counts") {
                    rel_ab <- sweep(raw_mat, 2, colSums(raw_mat), "/")
                    pa_matrix <- ifelse(rel_ab > pa_threshold, 1, 0)
                } else {
                    pa_matrix <- ifelse(raw_mat > pa_threshold, 1, 0)
                }
            } else {
                pa_matrix <- ifelse(raw_mat > pa_threshold, 1, 0)
            }
        } else {
            stop("pa_threshold must be a numeric value or a function.")
        }
    }

    if (ncol(pa_matrix) != nrow(meta_df)) {
        stop("Columns in abundance data must match rows in metadata.")
    }

    pa_cols <- colnames(pa_matrix)
    meta_rows <- rownames(meta_df)
    if (!is.null(pa_cols) && !is.null(meta_rows)) {
        if (!identical(pa_cols, meta_rows)) {
            stop("Colnames of abundance data must match rownames of metadata.")
        }
    } else {
        warning("Missing col/rownames. Assuming identically ordered samples.")
    }

    # 3. Extract var_of_interest if not explicitly provided
    if (is.null(var_of_interest)) {
        term_labels <- attr(stats::terms(formula), "term.labels")
        if (length(term_labels) == 0) {
            stop("The formula has no predictors.")
        }
        var_of_interest <- term_labels[1]
        msg <- sprintf("Using '%s' as the variable of interest.",
                       var_of_interest)
        message(msg)
    }

    # 4. Read Depth Handling
    read_depth_var <- NULL
    if (isTRUE(read_depth)) {
        if (assay_type != "counts") {
            stop(paste0(
                "read_depth = TRUE requires assay_type = 'counts'. ",
                "Provide read depth variable name manually or set FALSE."
            ))
        }
        meta_df$log10_read_depth <- log10(colSums(raw_mat))
        read_depth_var <- "log10_read_depth"

        formula_str <- deparse(formula)
        if (!grepl("log10_read_depth", formula_str)) {
            formula <- stats::update(formula, ~ . + log10_read_depth)
        }
    } else if (is.character(read_depth)) {
        if (!read_depth %in% colnames(meta_df)) {
            stop(sprintf("Read depth variable '%s' not found.", read_depth))
        }
        read_depth_var <- read_depth
    }

    # 5. Standardize and center numeric covariates
    vars_in_model <- all.vars(formula)
    for (var in vars_in_model) {
        if (is.numeric(meta_df[[var]])) {
            if (!is.null(read_depth_var) && var == read_depth_var) {
                meta_df[[var]] <- as.numeric(
                    scale(meta_df[[var]], center = TRUE, scale = FALSE)
                )
            } else if (standardize) {
                meta_df[[var]] <- as.numeric(scale(meta_df[[var]]))
            }
        }
    }

    # 6. Filtering Taxa
    N <- ncol(pa_matrix)
    th_pres <- ifelse(
        min_present > 0 && min_present < 1,
        ceiling(min_present * N),
        min_present
    )
    th_abs <- ifelse(
        min_absent > 0 && min_absent < 1,
        ceiling(min_absent * N),
        min_absent
    )

    if (th_pres > 0 || th_abs > 0) {
        prevalence <- rowSums(pa_matrix)
        keep_taxa <- prevalence >= th_pres & (N - prevalence) >= th_abs
        pa_matrix <- pa_matrix[keep_taxa, , drop = FALSE]

        if (nrow(pa_matrix) == 0) {
            stop("All taxa filtered out. Adjust min_present or min_absent.")
        }
    }

    # 7. Build Design Matrix and Check NAs
    meta_sub <- meta_df[, vars_in_model, drop = FALSE]
    complete_cases <- stats::complete.cases(meta_sub)
    if (sum(!complete_cases) > 0) {
        stop(sprintf(
            "Metadata contains NAs. %d sample(s) affected.",
            sum(!complete_cases)
        ))
    }

    X_full <- stats::model.matrix(formula, data = meta_df)
    X_design <- X_full[, -1, drop = FALSE] # Drop intercept

    # Ensure var_of_interest is in the design matrix, move to first column
    var_cols <- grep(var_of_interest, colnames(X_design))
    if (length(var_cols) == 0) {
        stop(sprintf("Variable '%s' not found.", var_of_interest))
    }

    other_cols <- setdiff(1:ncol(X_design), var_cols)
    X_design <- X_design[, c(var_cols, other_cols), drop = FALSE]

    # 8. Center all columns (including dummy variables)
    if (center_dummies) {
        col_means <- colMeans(X_design)
        X_design <- sweep(X_design, 2, col_means, FUN = "-")
    }

    # Extract levels for the variable of interest if categorical
    var_levels <- NULL
    if (is.character(meta_df[[var_of_interest]]) ||
        is.factor(meta_df[[var_of_interest]])) {
        var_levels <- levels(as.factor(meta_df[[var_of_interest]]))
    }

    # 9. Output Object
    list(
        y = as.matrix(pa_matrix),
        X = X_design,
        N = ncol(pa_matrix),
        K = nrow(pa_matrix),
        P = ncol(X_design),
        taxa_names = rownames(pa_matrix),
        sample_names = colnames(pa_matrix),
        design_matrix_cols = colnames(X_design),
        formula = formula,
        var_of_interest = var_of_interest,
        var_levels = var_levels,
        read_depth_var = read_depth_var
    )
}
