#' C. elegans SNV set
#'
#' A data set that contains bi-allelic markers in the C. elegans population.
#' Values of -1 corresponds to the REF allele, 1 correspond to the ALT allele
#'
#' @format A data frame with 37851 rows, 4 variables, and 249 strains (as variables):
#' \describe{
#'   \item{CHROM}{Chromosome of SNV}
#'   \item{POS}{Physical position of SNV}
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternate allele}
#'   \item{NIC511}{The first of 249 strains with genotypes}
#'   ...
#' }
#' @source \url{https://elegansvariation.org/}
"snps"

#' C. elegans relatedness matrix
#'
#' An N x N matrix containing the additive relatedness matrix in the C. elegans population
#'
#' @format A data frame with 249 rows, 249 variables:
#' \describe{
#'   \item{Strain1}{The relatedness of Strain1 to all other strains}
#'   \item{Strain2}{The relatedness of Strain2 to all other strains}
#'   ...
#' }
#' @source \url{https://elegansvariation.org/}
"kinship"
