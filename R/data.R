#' Gut microbiome dataset from Hintikka et al. (2021)
#'
#' A \code{TreeSummarizedExperiment} object containing a subset of gut
#' microbiome data from an experiment studying the effects of a high-fat
#' diet and xylo-oligosaccharide (XOS) supplementation on rats (N = 20 + 20).
#' The abundance data are sequencing counts and they have been agglomerated
#' to the genus level.
#'
#' @format A \code{TreeSummarizedExperiment} object.
#' \describe{
#'   \item{assay}{Contains the \code{counts} matrix.}
#'   \item{colData}{Contains the sample metadata, including \code{Fat}
#'   (diet group, factor with levels Low and High) and \code{XOS}
#'   (xylo-oligosaccharide supplementation).}
#' }
#' @source \url{https://www.mdpi.com/1660-4601/18/8/4049}
"tse_hintikka"

#' Example model fit object for DiPPER
#'
#' A pre-calculated \code{dipper_fit} object to be used in examples and
#' vignettes to avoid long MCMC sampling times during package checks.
#'
#' @format A list of class \code{dipper_fit} containing:
#' \describe{
#'   \item{stanfit}{The \code{CmdStanMCMC} object returned by CmdStanR.}
#'   \item{dipper_data}{A list containing the prepared data passed to Stan.}
#'   \item{symmetric}{Logical indicating if a symmetric Laplace prior for
#'   differential prevalence parameters was used.}
#' }
#'
#' @source Generated using the \code{tse_hintikka} dataset with
#'   \code{niter = 400} and \code{chains = 1}.
#'
#' @usage data("fit_example")
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
#' @source Data obtained via the \code{curatedMetagenomicData} R package.
#' Original study: \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC4950857/}
"VatanenT_2016_subset"
