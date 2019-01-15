#' Perform GWA analysis
#'
#' \code{perform_mapping} runs association mapping on a phenotype
#'
#' @param phenotype a data frame with columns: strain and phenotype [\strong{Default:} \code{"wormbase_gene"}]
#' @param genotype a genotype matrix in the format CHROM, POS, REF, ALT, strain1, ... , strainN
#' [\strong{Default:} \code{cegwas2::snps}]
#' @param kinship a N x N relatedness matrix in the format of [\strong{Default:} \code{cegwas2::kinship}]
#' @param P3D TRUE/FALSE - TRUE refers to the EMMAx algorithm  FALSE refers to EMMA algorithm
#' @param n.PC Integer describing the number of principal components
#' of the genotype matrix to use as fixed effects, [\strong{Default:0}]
#' @param min.MAF a value ranging from 0 - 1 that determines
#' the minimum minor allele frequencey a marker must have to be
#' used in association mapping [\strong{Default:} \code{0.05}]
#' @param map_by_chrom TRUE/FALSE - BLUP residual mappings from [\strong{Bloom, J. S. et al. 2015}]
#' [\strong{Default:} \code{FALSE}]
#' @param mapping_cores Integer, number of cores to assign to mapping
#' @param FDR_threshold threshold for calculated FDR and Bonferroni, [\strong{Default:0.05}]
#' @return a dataframe with the following columns
#' \itemize{
#'      \item \strong{CHROM} - Chromosome name
#'      \item \strong{POS} - Physical position of marker
#'      \item \strong{marker} - Marker name
#'      \item \strong{trait} - Trait name of input phenotype
#'      \item \strong{BF} - Bonferroni-corrected p-value threshold
#'      \item \strong{log10p} - [\code{-log10}] transformation of the p-value for the indicated marker
#'      \item \strong{pval} - p-value for indicated marker
#'      \item \strong{Zscore} - standard score for all p-values
#'      \item \strong{qvalue} - p-values adjusted by FDR
#' }
#' @export
perform_mapping <- function(phenotype = NULL,
                            genotype = cegwas2::snps,
                            kinship = cegwas2::kinship,
                            P3D = FALSE,
                            n.PC = 0,
                            min.MAF = 0.05,
                            map_by_chrom = FALSE,
                            mapping_cores = parallel::detectCores(),
                            FDR_threshold = 0.05) {

    # Vairable descriptions:
    # Y = phenotype
    # K = kinship matrix
    # M = genotype matrix

    if (any(colnames(genotype)[1:4] != c("CHROM", "POS", "REF", "ALT"))){
        stop(message(glue::glue("The genotype matrix is not formatted correctly. Please refer to documentation")))
    }

    if (any(colnames(kinship) != rownames(kinship))){
        stop(message(glue::glue("The kinship matrix is not formatted correctly. Please refer to documentation")))
    }

    #=============================#
    # Clean Phenotypes            #
    #=============================#
    Y = na.omit(phenotype)
    colnames(Y) <- c("strain", "trait")

    #=============================#
    # Process Kinship             #
    #=============================#
    # Sort kinship matrix to to be in the same order as phenotyped strains
    kinship_sorted = kinship[sort(row.names(kinship)), sort(colnames(kinship))]
    # Prune kinship matrix to only include phenotyped individuals
    K = kinship_sorted[row.names(kinship_sorted) %in% Y$strain,
                       colnames(kinship_sorted) %in% Y$strain]

    #=============================#
    # Process Markers             #
    #=============================#
    # Sort kinship matrix to to be in the same order as phenotyped strains
    markers = genotype %>%
        tidyr::unite(marker, CHROM, POS, remove = FALSE ) %>%
        dplyr::select(marker, CHROM, POS, dplyr::everything(), -REF, -ALT)

    markers_sorted <- data.frame(markers[,1:3],
                                 markers[, sort(names(markers[names(markers) %in% Y$strain]))])

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

    #=============================#
    # Perform Mapping             #
    #=============================#
    if (map_by_chrom) {
        # need 6 different kinship matrices
        # chrom2:x, map on one
        # chrom1:5, map on x, etc
        # then subtract blups from phenotype and map on
        # chromosome that you are not accounting for relationship

        # I think these are two ways to correct phenotype for relatedness
        yfit <- rrBLUP::mixed.solve(y = Y$trait, K = K)
        kfit <- rrBLUP::kin.blup(data = data.frame(Y),
                                 geno = "strain",
                                 pheno = "trait",
                                 K = K)

        # then get blup residuals
        y_blup <- Y$trait - yfit$u
        # then map outside the EMMA framework, similar to linkage mapping
        # (−n(ln(1−r2)/2ln(10)))
        # repeat for all chromosomes
    } else {
        # Perform mapping
        gwa_results = rrBLUP::GWAS(pheno = data.frame(Y),
                                   geno = data.frame(M),
                                   K = K,
                                   n.PC = 0,
                                   min.MAF = min.MAF,
                                   n.core = mapping_cores,
                                   P3D = P3D,
                                   plot = FALSE)
    }

    #=============================#
    # Process Mapping             #
    #=============================#
    gwa_results_pr <- gwa_results %>%
        dplyr::rename(log10p = trait) %>%
        dplyr::filter(log10p != 0) %>%
        dplyr::mutate(BF = -log10(FDR_threshold/n()), # set bonferroni threshold
                      trait = colnames(Y)[2]) %>%
        dplyr::select(CHROM, POS, marker, trait, BF, log10p) %>%
        dplyr::mutate(pval = 10^-log10p) %>% # convert log10p to p
        dplyr::mutate(Zscore = (pval - mean(pval)) / sd(pval), # calculate z-score
                      qvalue = qvalue(pval)) %>% # calculate q-value
        dplyr::arrange(CHROM,POS)

    #=============================#
    # Calculate FDR               #
    #=============================#
    q.ans <- qvalue(gwa_results_pr$pval)
    qv.sc <- cbind(q.ans, gwa_results_pr$log10p)
    qv.sc <- qv.sc[order(qv.sc[,1]),]

    if (qv.sc[1,1] < FDR_threshold) {
        avg.score <- tapply(qv.sc[,2], qv.sc[,1], mean)  #take mean of log10p for every unique qvalue
        qvals <- as.numeric(rownames(avg.score)) # unique q values
        x <- which.min(abs(qvals - FDR_threshold)) # find smalles q-value
        first <- max(1, x-2)
        last <- min(x+2, length(qvals))
        if ((last - first) < 4) {last <- first + 3}
        splin <- smooth.spline(x=qvals[first:last], y=avg.score[first:last], df=3)
        gwa_results_pr$fdr <- predict(splin, x=FDR_threshold)$y
    } else {
        gwa_results_pr$fdr <- NA
    }

    return(gwa_results_pr)
}

# Taken from rrBLUP GWAS function ~
# https://github.com/cran/rrBLUP/blob/master/R/GWAS.R
qvalue <- function(p) {
    smooth.df = 3

    if (min(p) < 0 || max(p) > 1) {
        print("ERROR: p-values not in valid range.")
        return(0)
    }

    lambda <- seq(0, 0.90, 0.05)
    m <- length(p)

    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
        pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
    }

    spi0 <- smooth.spline(lambda, pi0, df=smooth.df)
    pi0 <- predict(spi0, x=max(lambda))$y
    pi0 <- min(pi0, 1)

    if (pi0 <= 0) {
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
        return(0)
    }

    # The estimated q-values calculated here
    u <- order(p)

    # Ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
        idx <- sort.list(x)

        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)

        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl

        return(tbl)
    }

    v <- qvalue.rank(p)

    qvalue <- pi0*m*p/v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)

    for(i in (m-1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i+1]], 1)
    }

    return(qvalue)
}
