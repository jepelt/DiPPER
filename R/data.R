#' Rat microbiome data: low/high-fat diet and XOS supplementation
#'
#' A \code{TreeSummarizedExperiment} object containing microbiome data from a
#' rat experiment studying the effects of a high-fat diet and prebiotic
#' xylo-oligosaccharide (XOS) supplementation.
#'
#' The dataset has been agglomerated to the Genus level using
#' \code{mia::agglomerateByRank()}.
#'
#' @format A \code{TreeSummarizedExperiment} object with 40 samples (rats).
#' \describe{
#'   \item{colData}{Contains experimental metadata. The included variables are
#'   the diet group (\code{Fat}: High/Low) and supplementation (\code{XOS}:
#'   No / Yes).}
#'   \item{assays}{Contains the \code{counts} assay with microbial abundances.}
#' }
#' @source Hintikka L. et al. (2021). Xylo-Oligosaccharides in Prevention of
#' Hepatic Steatosis and Adipose Tissue Inflammation: Associating Taxonomic and
#' Metabolomic Patterns in Fecal Microbiomes with Biclustering.
#' \emph{International Journal of Environmental Research and Public Health},
#' 18(8):4049; https://doi.org/10.3390/ijerph18084049.
"tse_hintikka"

#' Longitudinal infant/child gut microbiome data
#'
#' A \code{TreeSummarizedExperiment} object containing a lightweight subset
#' (40 subjects, 1-3 measurements per subject) of the Vatanen et al. (2016)
#' longitudinal infant/child microbiome dataset. This dataset tracks the
#' maturation of the infant intestinal microbiome over the first
#' year of life.
#'
#' The dataset has been agglomerated to
#' the Genus level using \code{mia::agglomerateByRank()}. The assay contains
#' relative abundances expressed as proportions (0-1).
#'
#' @format A \code{TreeSummarizedExperiment} object with 76 samples.
#' \describe{
#'   \item{colData}{Contains experimental metadata. Key variables include
#'   \code{subject_id} (unique individual identifier), \code{age}
#'   (in months), \code{antibiotics} (no/yes), and \code{gender}
#'   (male/female).}
#'   \item{assays}{Contains the \code{relative_abundance} assay.}
#' }
#' @source Vatanen T. et al. (2016). Variation in Microbiome LPS Immunogenicity
#' Contributes to Autoimmunity in Humans. \emph{Cell},
#' 165(4):842-853; https://doi.org/10.1016/j.cell.2016.04.007.
"VatanenT_2016_subset"
