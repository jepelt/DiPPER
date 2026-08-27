#' Gut microbiome dataset from Hintikka et al. (2021)
#'
#' A \code{TreeSummarizedExperiment} object containing a subset of gut
#' microbiome data from an experiment studying the effects of a high-fat
#' diet and xylo-oligosaccharide (XOS) supplementation on rats (N = 20 + 20).
#' The abundance data are sequencing counts and they have been agglomerated
#' to the genus level.z
#'
#' @format A \code{TreeSummarizedExperiment} object.
#' \describe{
#'   \item{assay}{Contains the \code{counts} matrix.}
#'   \item{colData}{Contains the sample metadata, including \code{Fat}
#'   (diet group, factor with levels Low and High) and \code{XOS}
#'   (xylo-oligosaccharide supplementation).}
#' }
#' @source Derived from the \code{HintikkaXOData} dataset in the
#'   \code{mia} package. See \code{inst/scripts/tse_hintikka.R} for how the
#'   subset was derived. Original study:
#'   \url{https://www.mdpi.com/1660-4601/18/8/4049}
"tse_hintikka"

#' Example DiPPER model fit object
#'
#' A pre-calculated \code{dipper_fit} object based on tse_hintikka data. To be
#' used in examples to avoid long MCMC sampling times during package checks.
#' This object was created to enable running the examples of \code{summary},
#' \code{plot} and \code{print} without using CmdStan. The number of iterations
#' \code{niter = 400} and chains \code{chains = 2} used to create this object is
#' too low for reliable inference in practice.
#' @format A list of class \code{dipper_fit} containing:
#' \describe{
#'   \item{draws}{A matrix of posterior draws, with one row per MCMC draw and
#'   one column per taxon (the beta parameters).}
#'   \item{dipper_data}{A list containing the prepared data passed to Stan.}
#'   \item{symmetric}{Logical indicating if a symmetric Laplace prior for
#'   the parameters of interest was used.}
#'   \item{diagnostics}{A list of MCMC diagnostics (divergent transitions,
#'   maximum treedepth hits, sampling time, convergence statistics,
#'   and diagnostics level).}
#' }
#' @source Generated using the \code{tse_hintikka} dataset with
#'   \code{niter = 400} and \code{chains = 2}.
"fit_example"

#' Longitudinal gut microbiome dataset from Vatanen et al. (2016)
#'
#' A \code{TreeSummarizedExperiment} object containing a longitudinal subset of
#' infant gut microbiome data from Vatanen et al. (2016). The data consists of
#' 79 infants with measurements at 1.5 and/or 7 months of age. The abundance
#' data are relative abundances (proportions) and they have been agglomerated
#' to the genus level.
#'
#' @format A \code{TreeSummarizedExperiment} object.
#' \describe{
#'   \item{assay}{Contains the \code{relative_abundance} matrix (proportions).}
#'   \item{colData}{Contains the sample metadata: \code{subject_id} (factor),
#'   \code{antibiotics} (factor: no, yes), \code{age_point} (factor:
#'   1.5 months, 7 months), and \code{gender} (factor: male, female).}
#'   \item{rowTree}{Contains the phylogenetic tree.}
#' }
#' @source Derived from the \code{VatanenT_2016} study, obtained via the
#'   \code{curatedMetagenomicData} package. See
#'   \code{inst/scripts/VatanenT_2016_subset.R} for how the subset was
#'   derived. Original study:
#'   \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC4950857/}
"VatanenT_2016_subset"
