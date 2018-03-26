#' Perform GWA analysis
#'
#' \code{perform_mapping} runs association mapping on a phenotype
#'
#' @param phenotype a data frame with columns: strain and phenotype [\strong{Default:} \code{"wormbase_gene"}]
#' @param genotype a genotype matrix in the format [\strong{Default:} \code{cegwas2::snps}]
#' @param kinship a N x N relatedness matrix in the format [\strong{Default:} \code{cegwas2::kinship}]
#' @param P3D TRUE/FALSE - TRUE refers to the EMMAx algorithm, FALSE refers to EMMA algorithm
#' [\strong{Default:} \code{FALSE}]
#' @param MAF a value ranging from 0 - 1 that determines the minimum minor allele frequencey
#' a marker must have to be used in association mapping [\strong{Default:} \code{0.05}]
#' @examples get_db()
#' @return a dataframe with the following columns
#' \itemize{
#'      \item \strong{CHROM} - Chromosome name
#'      \item \strong{POS} - Physical position of marker
#'      \item \strong{marker} - Marker name
#'      \item \strong{trait} - Trait name of input phenotype
#'      \item \strong{BF} - Bonferroni-corrected p-value threshold
#'      \item \strin{log10p} - [\code{-log10}] transformation of the p-value for the indicated marker
#' }
#' @export

perform_mapping <- function(phenotype = NULL,
                            genotype = cegwas2::snps,
                            kinship = cegwas2::kinship,
                            P3D = FALSE,
                            min.MAF = 0.05) {

    # Clean phenotypes
    Y = na.omit(phenotype)
    colnames(Y) <- c("strain", "trait")

    # Sort kinship matrix to to be in the same order as phenotyped strains
    kinship_sorted = kinship[sort(row.names(kinship)),sort(colnames(kinship))]
    # Prune kinship matrix to only include phenotyped individuals
    K = kinship_sorted[row.names(kinship_sorted) %in% Y$strain, colnames(kinship_sorted) %in% Y$strain]
    # Sort kinship matrix to to be in the same order as phenotyped strains
    markers = genotype %>%
        tidyr::unite(marker, CHROM, POS, remove = FALSE ) %>%
        dplyr::select(marker, CHROM, POS, everything(), -REF, -ALT)

    markers_sorted <- data.frame(markers[,1:3],
                                 markers[, sort(names(markers[names(markers)%in%Y$strain]))])

    # ID markers that have a MAF outside the user-defined range
    keepMarkers <- data.frame(
        MAF = apply(markers_sorted[,4:ncol(markers_sorted)],
                    MARGIN = 1,
                    FUN = function(x){
                        x[x==-1] <- 0
                        return(sum(x, na.rm = T)/length(x))}
        )) %>%
        dplyr::mutate(marker = markers_sorted$marker) %>%
        dplyr::filter(MAF >= min.MAF) %>%
        dplyr::filter(MAF <= 1 - min.MAF)

    # Remove markers identified to be out of MAF range
    M <- markers_sorted %>%
        dplyr::filter(marker %in% keepMarkers$marker)

    # Perform mapping
    gwa_results = rrBLUP::GWAS(pheno = data.frame(Y),
                               geno = data.frame(M),
                               K = K,
                               n.PC = 0,
                               min.MAF = min.MAF,
                               n.core = parallel::detectCores(),
                               P3D = P3D)

    # Process mapping results
    gwa_results_pr <- gwa_results %>%
        dplyr::rename(log10p = trait) %>%
        dplyr::filter(log10p != 0) %>%
        dplyr::mutate(BF = log10(.05/n()),
                      trait = colnames(Y)[2]) %>%
        dplyr::select(CHROM, POS, marker, trait, BF, log10p)

    return(gwa_results_pr)
}



