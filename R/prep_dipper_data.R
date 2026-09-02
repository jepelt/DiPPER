#' Extract the abundance matrix and metadata from the user input
#'
#' Accepts either a TreeSummarizedExperiment object or a separate matrix and
#' metadata data.frame, and returns both in a common form.
#'
#' @return A list with elements `raw_mat` and `meta_df`.
#'
#' @keywords internal
#' @noRd
.dipper_extract_input <- function(tse, assay, meta, assay.type) {

    if (!is.null(tse)) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop(
                "Package 'SummarizedExperiment' is required for",
                "SummarizedExperiment and TreeSummarizedExperiment input.",
                call. = FALSE
            )
        }

        if (is.null(assay.type)) {
            available_assays <- SummarizedExperiment::assayNames(tse)
            stop(
                "You must specify 'assay.type' when using a TSE object.\n",
                "Available assays in your object: ",
                paste(available_assays, collapse = ", "),
                call. = FALSE
            )
        }

        if (!assay.type %in% SummarizedExperiment::assayNames(tse)) {
            stop(
                "Assay '", assay.type, "' not found. Available assays: ",
                paste(SummarizedExperiment::assayNames(tse), collapse = ", "),
                call. = FALSE
            )
        }

        return(list(
            raw_mat = SummarizedExperiment::assay(tse, assay.type),
            meta_df = as.data.frame(SummarizedExperiment::colData(tse))
        ))
    }

    if (!is.null(meta) && !is.null(assay)) {
        return(list(
            raw_mat = as.matrix(assay),
            meta_df = as.data.frame(meta)
        ))
    }

    stop("Provide either 'tse', OR both 'meta' and 'assay'.", call. = FALSE)
}


#' Validate the abundance matrix against the declared data type
#'
#' @return Invisibly `NULL`. Called for the side effect of stopping on
#'   invalid input.
#'
#' @keywords internal
#' @noRd
.dipper_validate_abundance <- function(raw_mat, data.type) {

    if (anyNA(raw_mat)) {
        stop(
            "The abundance data contains missing values (NA).",
            call. = FALSE
        )
    }

    if (data.type == "counts") {
        if (any(raw_mat < 0) || any(raw_mat %% 1 != 0)) {
            stop(
                "Count data must contain non-negative integers.",
                call. = FALSE
            )
        }
    } else if (data.type == "relabundance") {
        if (any(raw_mat < 0) || any(raw_mat > 1)) {
            stop(
                "Relative abundance values must be between 0 and 1.",
                call. = FALSE
            )
        }
        if (any(colSums(raw_mat) > 1.0001)) {
            stop(
                "Column sums of relative abundances must be <= 1.",
                call. = FALSE
            )
        }
    } else if (data.type == "pa") {
        if (!all(raw_mat %in% c(0, 1))) {
            stop(
                "Presence/absence data must strictly contain 0s and 1s.",
                call. = FALSE
            )
        }
    }

    invisible(NULL)
}


#' Convert an abundance matrix into a presence/absence matrix
#'
#' @return A numeric matrix of zeros and ones.
#'
#' @keywords internal
#' @noRd
.dipper_to_presence_absence <- function(raw_mat, data.type, threshold) {

    if (data.type == "pa") {
        return(raw_mat)
    }

    if (is.function(threshold)) {
        row_thresh <- apply(raw_mat, 1, threshold)
        return(ifelse(raw_mat > row_thresh, 1, 0))
    }

    if (!is.numeric(threshold)) {
        stop(
            "threshold must be a numeric value or a function.",
            call. = FALSE
        )
    }

    if (threshold < 0) {
        stop("threshold cannot be negative.", call. = FALSE)
    }

    if (threshold > 0 && threshold < 1 && data.type == "counts") {
        rel_ab <- sweep(raw_mat, 2, colSums(raw_mat), "/")
        return(ifelse(rel_ab > threshold, 1, 0))
    }

    ifelse(raw_mat > threshold, 1, 0)
}


#' Check that abundance columns and metadata rows refer to the same samples
#'
#' @return Invisibly `NULL`.
#'
#' @keywords internal
#' @noRd
.dipper_check_alignment <- function(pa_matrix, meta_df) {

    if (ncol(pa_matrix) != nrow(meta_df)) {
        stop(
            "Columns in abundance data must match rows in metadata.",
            call. = FALSE
        )
    }

    pa_cols <- colnames(pa_matrix)
    meta_rows <- rownames(meta_df)

    if (!is.null(pa_cols) && !is.null(meta_rows)) {
        if (!identical(pa_cols, meta_rows)) {
            stop(
                "Colnames of abundance data must match rownames of metadata.",
                call. = FALSE
            )
        }
    } else {
        warning(
            "Missing col/rownames. Assuming identically ordered samples.",
            call. = FALSE
        )
    }

    invisible(NULL)
}


#' Split a model formula into its fixed part and a random intercept
#'
#' Only a single random intercept of the form `(1 | id)` is supported.
#'
#' @return A list with elements `is_longitudinal`, `id_var` and
#'   `fixed_formula`.
#'
#' @keywords internal
#' @noRd
.dipper_parse_random_effect <- function(formula, meta_df) {

    f_str <- paste(deparse(formula), collapse = " ")
    re_pattern <- "\\(\\s*1\\s*\\|\\s*([^\\)]+)\\s*\\)"

    if (!grepl(re_pattern, f_str)) {
        if (grepl("\\|", f_str)) {
            stop(
                "Unsupported random effect structure. Use (1 | id).",
                call. = FALSE
            )
        }
        return(list(
            is_longitudinal = FALSE,
            id_var = NULL,
            fixed_formula = formula
        ))
    }

    all_re <- gregexpr("\\([^\\)]+\\|[^\\)]+\\)", f_str)
    if (length(regmatches(f_str, all_re)[[1]]) > 1) {
        stop(
            "Only a single random intercept (1 | id) is supported.",
            call. = FALSE
        )
    }

    id_var <- trimws(sub(paste0(".*", re_pattern, ".*"), "\\1", f_str))

    if (!id_var %in% colnames(meta_df)) {
        stop("Subject ID variable '", id_var, "' not found.", call. = FALSE)
    }

    fixed_str <- sub(re_pattern, "", f_str)
    fixed_str <- gsub("\\+\\s*\\+", "+", fixed_str)
    fixed_str <- gsub("\\+\\s*$", "", fixed_str)
    fixed_str <- gsub("~\\s*\\+", "~", fixed_str)

    list(
        is_longitudinal = TRUE,
        id_var = id_var,
        fixed_formula = stats::as.formula(fixed_str)
    )
}


#' Check the formula variables and resolve the variable of interest
#'
#' @return The name of the variable of interest.
#'
#' @keywords internal
#' @noRd
.dipper_resolve_var_of_interest <- function(fixed_formula,
                                            meta_df,
                                            var.of.interest) {

    vars_in_model <- all.vars(fixed_formula)
    missing_vars <- setdiff(vars_in_model, colnames(meta_df))

    if (length(missing_vars) > 0) {
        stop(
            "The following variables are missing from the metadata: ",
            paste(missing_vars, collapse = ", "),
            call. = FALSE
        )
    }

    if (is.null(var.of.interest)) {
        term_labels <- attr(stats::terms(fixed_formula), "term.labels")
        if (length(term_labels) == 0) {
            stop("The formula has no predictors.", call. = FALSE)
        }
        var.of.interest <- term_labels[1]
        message("Using '", var.of.interest, "' as the variable of interest.")
    }

    if (!var.of.interest %in% colnames(meta_df)) {
        stop(
            "Variable '", var.of.interest, "' not found in metadata.",
            call. = FALSE
        )
    }

    var.of.interest
}


#' Add the read depth covariate to the metadata and the model formulas
#'
#' @return A list with elements `meta_df`, `formula`, `fixed_formula` and
#'   `read.depth.var`.
#'
#' @keywords internal
#' @noRd
.dipper_add_read_depth <- function(read.depth,
                                   data.type,
                                   raw_mat,
                                   meta_df,
                                   formula,
                                   fixed_formula) {

    read.depth.var <- NULL

    if (isTRUE(read.depth)) {
        if (data.type != "counts") {
            stop(
                "read.depth = TRUE requires data.type = 'counts'. ",
                "Provide read depth variable name manually or set FALSE.",
                call. = FALSE
            )
        }
        meta_df$log10_read.depth <- log10(colSums(raw_mat))
        read.depth.var <- "log10_read.depth"

        formula_str <- paste(deparse(fixed_formula), collapse = " ")
        if (!grepl("log10_read.depth", formula_str)) {
            formula <- stats::update(formula, ~ . + log10_read.depth)
            fixed_formula <- stats::update(
                fixed_formula, ~ . + log10_read.depth
            )
        }
    } else if (is.character(read.depth)) {
        if (!read.depth %in% colnames(meta_df)) {
            stop(
                "Read depth variable '", read.depth,
                "' not found in metadata.",
                call. = FALSE
            )
        }
        read.depth.var <- read.depth
    }

    list(
        meta_df = meta_df,
        formula = formula,
        fixed_formula = fixed_formula,
        read.depth.var = read.depth.var
    )
}


#' Standardize numeric covariates and record the scales used
#'
#' The read depth covariate is only centered, not scaled, so that its
#' coefficient stays on the log10 scale.
#'
#' @return A list with elements `meta_df` and `continuous_scales`.
#'
#' @keywords internal
#' @noRd
.dipper_scale_covariates <- function(meta_df, vars_to_scale, read.depth.var) {

    continuous_scales <- list()

    for (var in vars_to_scale) {
        if (!is.numeric(meta_df[[var]])) {
            next
        }

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

    list(meta_df = meta_df, continuous_scales = continuous_scales)
}


#' Filter taxa by prevalence
#'
#' Thresholds given between 0 and 1 are interpreted as proportions of samples,
#' larger values as absolute sample counts.
#'
#' @return The filtered presence/absence matrix.
#'
#' @keywords internal
#' @noRd
.dipper_filter_taxa <- function(pa_matrix, min.present, min.absent) {

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

    if (th_pres <= 0 && th_abs <= 0) {
        return(pa_matrix)
    }

    prevalence <- rowSums(pa_matrix)
    keep_taxa <- prevalence >= th_pres & (N - prevalence) >= th_abs
    pa_matrix <- pa_matrix[keep_taxa, , drop = FALSE]

    n_taxa_final <- nrow(pa_matrix)
    if (n_taxa_final == 0) {
        stop(
            "All taxa filtered out. Adjust min.present or min.absent.",
            call. = FALSE
        )
    }

    message(
        "Filtering: ", n_taxa_initial - n_taxa_final, " taxa removed, ",
        n_taxa_final, " taxa retained for analysis."
    )

    pa_matrix
}


#' Build the centered design matrix
#'
#' The intercept is dropped, the columns of the variable of interest are moved
#' to the front, and all columns are mean-centered.
#'
#' @return The design matrix.
#'
#' @keywords internal
#' @noRd
.dipper_build_design <- function(fixed_formula, meta_df, var.of.interest) {

    X_full <- stats::model.matrix(fixed_formula, data = meta_df)
    X_design <- X_full[, -1, drop = FALSE]

    var_cols <- grep(var.of.interest, colnames(X_design))
    if (length(var_cols) == 0) {
        stop(
            "Variable '", var.of.interest, "' not found in design matrix.",
            call. = FALSE
        )
    }

    other_cols <- setdiff(seq_len(ncol(X_design)), var_cols)
    X_design <- X_design[, c(var_cols, other_cols), drop = FALSE]

    col_means <- colMeans(X_design)
    sweep(X_design, 2, col_means, FUN = "-")
}


#' Determine and validate the levels of a categorical variable of interest
#'
#' @return The levels of the variable of interest, or `NULL` if it is not
#'   categorical.
#'
#' @keywords internal
#' @noRd
.dipper_var_levels <- function(meta_df, var.of.interest) {

    var_value <- meta_df[[var.of.interest]]

    if (!is.character(var_value) && !is.factor(var_value)) {
        return(NULL)
    }

    var_levels <- levels(as.factor(var_value))

    if (length(var_levels) > 2) {
        stop(
            "The variable of interest ('", var.of.interest, "') has ",
            length(var_levels), " levels (",
            paste(var_levels, collapse = ", "), ").\n",
            "DiPPER currently supports only binary or continuous ",
            "variables.\n",
            "Please combine levels, or subset your data to compare groups ",
            "pairwise (e.g., ", var_levels[1], " vs ", var_levels[2], ").",
            call. = FALSE
        )
    }

    if (length(var_levels) < 2) {
        stop(
            "Variable of interest ('", var.of.interest,
            "') must have at least two levels.",
            call. = FALSE
        )
    }

    var_levels
}


#' Prepare data for DiPPER
#'
#' @param tse A (Tree)SummarizedExperiment object.
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


    # 1. Input extraction ------------------------------------------------------
    input <- .dipper_extract_input(tse, assay, meta, assay.type)
    raw_mat <- input$raw_mat
    meta_df <- input$meta_df


    # 2. Validation and P/A conversion ----------------------------------------
    .dipper_validate_abundance(raw_mat, data.type)
    pa_matrix <- .dipper_to_presence_absence(raw_mat, data.type, threshold)
    .dipper_check_alignment(pa_matrix, meta_df)


    # 3. Parse longitudinal formula --------------------------------------------
    re <- .dipper_parse_random_effect(formula, meta_df)
    is_longitudinal <- re$is_longitudinal
    id_var <- re$id_var
    fixed_formula <- re$fixed_formula
    subj_vec <- NULL
    S <- 0


    # 4. Validate formula variables and extract var.of.interest ----------------
    vars_in_model <- all.vars(fixed_formula)
    var.of.interest <- .dipper_resolve_var_of_interest(
        fixed_formula, meta_df, var.of.interest
    )


    # 5. Read Depth Handling ---------------------------------------------------
    rd <- .dipper_add_read_depth(
        read.depth, data.type, raw_mat, meta_df, formula, fixed_formula
    )
    meta_df <- rd$meta_df
    formula <- rd$formula
    fixed_formula <- rd$fixed_formula
    read.depth.var <- rd$read.depth.var


    # 6. Standardize numeric covariates and track scales -----------------------
    scaled <- .dipper_scale_covariates(
        meta_df, setdiff(vars_in_model, id_var), read.depth.var
    )
    meta_df <- scaled$meta_df
    continuous_scales <- scaled$continuous_scales


    # 7. Filtering taxa/features based on prevalence ---------------------------
    pa_matrix <- .dipper_filter_taxa(pa_matrix, min.present, min.absent)


    # 8. Build design matrix and check NAs -------------------------------------
    vars_for_na_check <- vars_in_model
    if (is_longitudinal) {
        vars_for_na_check <- c(vars_for_na_check, id_var)
    }

    meta_sub <- meta_df[, vars_for_na_check, drop = FALSE]
    complete_cases <- stats::complete.cases(meta_sub)
    if (sum(!complete_cases) > 0) {
        stop(
            "Metadata contains NAs. ", sum(!complete_cases),
            " sample(s) affected.",
            call. = FALSE
        )
    }

    X_design <- .dipper_build_design(fixed_formula, meta_df, var.of.interest)


    # 9. Levels of the variable of interest ------------------------------------
    var_levels <- .dipper_var_levels(meta_df, var.of.interest)


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
