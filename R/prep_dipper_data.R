#' Prepare data for DiPPER
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param formula Model formula.
#' @param assay A matrix containing counts, relative abundances, or
#'   presence/absence data (required if tse is NULL). Rows should be
#'   taxa/features, columns should be samples.
#' @param meta A data.frame containing metadata (required if tse is NULL).
#' @param assay.type Character. Name of the assay in the TSE object to use.
#'   Required if tse is provided.
#' @param data.type Character. Type of the data: "counts" (default),
#'   "relabundance", or "pa" (presence/absence).
#' @param var.of.interest The variable of interest from the formula.
#'   If NULL, the first term of the formula is automatically used.
#' @param read.depth Logical or character. If TRUE, calculates log10 read depth
#'   from counts. If character, uses the specified metadata column.
#' @param threshold Numeric or function. Threshold for presence.
#'   If between 0 and 1, the threshold is based on relative abundance (e.g.,
#'   0.05 means that relative abundances > 0.05 are considered present).
#'   If whole number, the threshold is based on counts (e.g. 5 means that
#'   count > 5 is considered present).
#'   If a function (e.g. median), it is applied for the given assay for each
#'   feature. Default is 0. Ignored if data.type is "pa".
#' @param min.present For taxon/feature filtering. Minimum number or proportion
#'   of samples where a taxon must be present. Default is 5.
#' @param min.absent  For taxon/feature filtering. Minimum number or proportion
#'   of samples where a taxon must be absent. Default is min.present.
#'
#' @keywords internal
#' @importFrom stats model.matrix update complete.cases terms as.formula sd
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assay colData
prep_dipper_data <- function(tse = NULL,
                             formula,
                             assay = NULL,
                             meta = NULL,
                             assay.type = NULL,
                             data.type = c("counts", "relabundance", "pa"),
                             var.of.interest = NULL,
                             read.depth = TRUE,
                             threshold = 0,
                             min.present = 5,
                             min.absent = min.present) {

    data.type <- match.arg(data.type)

    raw_mat <- NULL
    meta_df <- NULL


    # 1. Input extraction ------------------------------------------------------
    if (!is.null(tse)) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop(paste0(
                "Package 'SummarizedExperiment' is required ",
                "for TSE objects."
            ))
        }

        if (is.null(assay.type)) {
            available_assays <- SummarizedExperiment::assayNames(tse)
            stop(sprintf(
                paste0(
                    "You must specify 'assay.type' when using a TSE object.\n",
                    "Available assays in your object: %s"
                ),
                paste(available_assays, collapse = ", ")
            ))
        }

        if (!assay.type %in% SummarizedExperiment::assayNames(tse)) {
            stop(sprintf(
                "Assay '%s' not found. Available assays: %s",
                assay.type,
                paste(SummarizedExperiment::assayNames(tse), collapse = ", ")
            ))
        }

        raw_mat <- SummarizedExperiment::assay(tse, assay.type)
        meta_df <- as.data.frame(SummarizedExperiment::colData(tse))

    } else if (!is.null(meta) && !is.null(assay)) {
        raw_mat <- as.matrix(assay)
        meta_df <- as.data.frame(meta)
    } else {
        stop("Provide either 'tse', OR both 'meta' and 'assay'.")
    }


    # 2. Validation and P/A conversion ----------------------------------------
    if (anyNA(raw_mat)) stop("The abundance data contains missing values (NA).")

    if (data.type == "counts") {
        if (any(raw_mat < 0) || any(raw_mat %% 1 != 0)) {
            stop("Count data must contain non-negative integers.")
        }
    } else if (data.type == "relabundance") {
        if (any(raw_mat < 0) || any(raw_mat > 1)) {
            stop("Relative abundance values must be between 0 and 1.")
        }
        if (any(colSums(raw_mat) > 1.0001)) {
            stop("Column sums of relative abundances must be <= 1.")
        }
    } else if (data.type == "pa") {
        if (!all(raw_mat %in% c(0, 1))) {
            stop("Presence/absence data must strictly contain 0s and 1s.")
        }
    }

    if (data.type == "pa") {
        pa_matrix <- raw_mat
    } else {
        if (is.function(threshold)) {
            row_thresh <- apply(raw_mat, 1, threshold)
            pa_matrix <- ifelse(raw_mat > row_thresh, 1, 0)
        } else if (is.numeric(threshold)) {
            if (threshold < 0) {
                stop("threshold cannot be negative.")
            } else if (threshold > 0 && threshold < 1) {
                if (data.type == "counts") {
                    rel_ab <- sweep(raw_mat, 2, colSums(raw_mat), "/")
                    pa_matrix <- ifelse(rel_ab > threshold, 1, 0)
                } else {
                    pa_matrix <- ifelse(raw_mat > threshold, 1, 0)
                }
            } else {
                pa_matrix <- ifelse(raw_mat > threshold, 1, 0)
            }
        } else {
            stop("threshold must be a numeric value or a function.")
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


    # 3. Parse longitudinal formula --------------------------------------------
    is_longitudinal <- FALSE
    id_var <- NULL
    subj_vec <- NULL
    S <- 0
    fixed_formula <- formula

    f_str <- paste(deparse(formula), collapse = " ")
    re_pattern <- "\\(\\s*1\\s*\\|\\s*([^\\)]+)\\s*\\)"

    if (grepl(re_pattern, f_str)) {
        is_longitudinal <- TRUE

        all_re <- gregexpr("\\([^\\)]+\\|[^\\)]+\\)", f_str)
        if (length(regmatches(f_str, all_re)[[1]]) > 1) {
            stop("Only a single random intercept (1 | id) is supported.")
        }

        id_var <- trimws(sub(paste0(".*", re_pattern, ".*"), "\\1", f_str))

        if (!id_var %in% colnames(meta_df)) {
            stop(sprintf("Subject ID variable '%s' not found.", id_var))
        }

        fixed_str <- sub(re_pattern, "", f_str)
        fixed_str <- gsub("\\+\\s*\\+", "+", fixed_str)
        fixed_str <- gsub("\\+\\s*$", "", fixed_str)
        fixed_str <- gsub("~\\s*\\+", "~", fixed_str)

        fixed_formula <- stats::as.formula(fixed_str)
    } else if (grepl("\\|", f_str)) {
        stop("Unsupported random effect structure. Use (1 | id).")
    }


    # 4. Validate formula variables and extract var.of.interest ----------------
    vars_in_model <- all.vars(fixed_formula)
    missing_vars <- setdiff(vars_in_model, colnames(meta_df))

    if (length(missing_vars) > 0) {
        stop(sprintf(
            "The following variables are missing from the metadata: %s",
            paste(missing_vars, collapse = ", ")
        ))
    }

    if (is.null(var.of.interest)) {
        term_labels <- attr(stats::terms(fixed_formula), "term.labels")
        if (length(term_labels) == 0) {
            stop("The formula has no predictors.")
        }
        var.of.interest <- term_labels[1]
        msg <- sprintf("Using '%s' as the variable of interest.",
                       var.of.interest)
        message(msg)
    }

    if (!var.of.interest %in% colnames(meta_df)) {
        stop(sprintf("Variable '%s' not found in metadata.", var.of.interest))
    }


    # 5. Read Depth Handling ---------------------------------------------------
    read.depth.var <- NULL
    if (isTRUE(read.depth)) {
        if (data.type != "counts") {
            stop(paste0(
                "read.depth = TRUE requires data.type = 'counts'. ",
                "Provide read depth variable name manually or set FALSE."
            ))
        }
        meta_df$log10_read.depth <- log10(colSums(raw_mat))
        read.depth.var <- "log10_read.depth"

        formula_str <- deparse(fixed_formula)
        if (!grepl("log10_read.depth", formula_str)) {
            formula <- stats::update(formula, ~ . + log10_read.depth)
            fixed_formula <- stats::update(
                fixed_formula, ~ . + log10_read.depth
            )
        }
    } else if (is.character(read.depth)) {
        if (!read.depth %in% colnames(meta_df)) {
            stop(sprintf("Read depth variable '%s' not found in metadata.",
                         read.depth))
        }
        read.depth.var <- read.depth
    }


    # 6. Standardize numeric covariates and track scales -----------------------
    vars_to_scale <- setdiff(vars_in_model, id_var)
    continuous_scales <- list()

    for (var in vars_to_scale) {
        if (is.numeric(meta_df[[var]])) {
            if (!is.null(read.depth.var) && var == read.depth.var) {
                meta_df[[var]] <- as.numeric(
                    scale(meta_df[[var]], center = TRUE, scale = FALSE)
                )
            } else {
                var_sd <- stats::sd(meta_df[[var]], na.rm = TRUE)
                if (var_sd > 0) {
                    continuous_scales[[var]] <- var_sd
                    meta_df[[var]] <- as.numeric(scale(meta_df[[var]]))
                }
            }
        }
    }


    # 7. Filtering taxa/features based on prevalence ---------------------------
    n_taxa_initial <- nrow(pa_matrix)
    N <- ncol(pa_matrix)
    th_pres <- ifelse(
        min.present > 0 && min.present < 1,
        ceiling(min.present * N),
        min.present
    )
    th_abs <- ifelse(
        min.absent > 0 && min.absent < 1,
        ceiling(min.absent * N),
        min.absent
    )

    if (th_pres > 0 || th_abs > 0) {
        prevalence <- rowSums(pa_matrix)
        keep_taxa <- prevalence >= th_pres & (N - prevalence) >= th_abs
        pa_matrix <- pa_matrix[keep_taxa, , drop = FALSE]

        n_taxa_final <- nrow(pa_matrix)
        if (n_taxa_final == 0) {
            stop("All taxa filtered out. Adjust min.present or min.absent.")
        }
        message(sprintf(
            "Filtering: %d taxa removed, %d taxa retained for analysis.",
            n_taxa_initial - n_taxa_final, n_taxa_final
        ))
    }


    # 8. Build design matrix and check NAs -------------------------------------
    vars_for_na_check <- vars_in_model
    if (is_longitudinal) {
        vars_for_na_check <- c(vars_for_na_check, id_var)
    }

    meta_sub <- meta_df[, vars_for_na_check, drop = FALSE]
    complete_cases <- stats::complete.cases(meta_sub)
    if (sum(!complete_cases) > 0) {
        stop(sprintf(
            "Metadata contains NAs. %d sample(s) affected.",
            sum(!complete_cases)
        ))
    }

    X_full <- stats::model.matrix(fixed_formula, data = meta_df)
    X_design <- X_full[, -1, drop = FALSE]

    var_cols <- grep(var.of.interest, colnames(X_design))
    if (length(var_cols) == 0) {
        stop(sprintf("Variable '%s' not found in design matrix.",
                     var.of.interest))
    }

    other_cols <- setdiff(1:ncol(X_design), var_cols)
    X_design <- X_design[, c(var_cols, other_cols), drop = FALSE]


    # 9. Center all columns ----------------------------------------------------
    col_means <- colMeans(X_design)
    X_design <- sweep(X_design, 2, col_means, FUN = "-")

    var_levels <- NULL
    is_categorical <- FALSE

    if (is.character(meta_df[[var.of.interest]]) ||
        is.factor(meta_df[[var.of.interest]])) {

        is_categorical <- TRUE
        var_levels <- levels(as.factor(meta_df[[var.of.interest]]))

        if (length(var_levels) > 2) {
            stop(sprintf(
                paste0(
                    "The variable of interest ('%s') has %d levels (%s).\n",
                    "DiPPER currently supports only binary or continuous ",
                    "variables.\n",
                    "Please combine levels, or subset your data to compare ",
                    "groups pairwise (e.g., %s vs %s)."
                ),
                var.of.interest,
                length(var_levels),
                paste(var_levels, collapse = ", "),
                var_levels[1],
                var_levels[2]
            ))
        } else if (length(var_levels) < 2) {
            stop(sprintf(
                "Variable of interest ('%s') must have at least two levels.",
                var.of.interest
            ))
        }
    }


    # 10. Process subject IDs for Stan -----------------------------------------
    if (is_longitudinal) {
        subj_vec <- as.integer(as.factor(meta_df[[id_var]]))
        S <- length(unique(subj_vec))
    }


    # 11. Output Object --------------------------------------------------------
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
        fixed_formula = fixed_formula,
        var.of.interest = var.of.interest,
        var_levels = var_levels,
        read.depth.var = read.depth.var,
        is_longitudinal = is_longitudinal,
        id_var = id_var,
        subj = subj_vec,
        S = S,
        continuous_scales = continuous_scales
    )
}
