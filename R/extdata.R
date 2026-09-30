#' Pre-computed DiPPER fits
#'
#' Two \code{dipper_fit} objects stored as \code{.rds} files under
#' \code{inst/extdata}. They are loaded by the package vignette when CmdStan
#' is unavailable, so that the vignette can be built and display real results
#' also on machines without a Stan toolchain.
#'
#' @return Each file contains a list of class \code{dipper_fit}. See the
#'   documentation of \code{dipper()} for a description of the elements.
#'   The files are read with \code{readRDS()}.
#'
#' @format \describe{
#'   \item{fit_hintikka.rds}{Fitted to the \code{\link{tse_hintikka}} dataset
#'   with the formula \code{~ Fat + XOS}, \code{assay.type = "counts"},
#'   \code{data.type = "counts"} and \code{read.depth = TRUE}.}
#'   \item{fit_vatanen.rds}{Fitted to the \code{\link{VatanenT_2016_subset}}
#'   dataset with the formula \code{~ age_point + antibiotics + gender +
#'   (1 | subject_id)}, \code{assay.type = "relative_abundance"},
#'   \code{data.type = "relabundance"} and \code{read.depth = FALSE}.}
#' }
#'
#' @source Generated from the bundled datasets. See
#'   \code{inst/scripts/fits_vignette.R} for the code that produced them.
#'
#' @seealso \code{\link{tse_hintikka}}, \code{\link{VatanenT_2016_subset}}
#'   and \code{\link{fit_example}} for the datasets and the example fit
#'   stored in \code{data/}.
#'
#' @examples
#' fit <- readRDS(
#'     system.file("extdata", "fit_hintikka.rds", package = "DiPPER")
#' )
#'
#' print(fit)
#' head(summary(fit))
#'
#' @name DiPPER-extdata
#' @aliases fit_hintikka fit_vatanen
#' @keywords datasets
NULL
