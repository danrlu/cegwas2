
perform_mapping <- function(phenotype = pr_phenotypes[,1:2],
                            genotype = cegwas2::snps,
                            kinship = cegwas2::kinship,
                            P3D = FALSE,
                            model = "additive",
                            MAF = 0.05,
                            FDR = 0.05) {

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

    # ID markers
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

    M <- markers_sorted %>%
        dplyr::filter(marker %in% keepMarkers$marker)

    gwa_results = rrBLUP::GWAS(pheno = data.frame(Y),
                               geno = data.frame(M),
                               K = K,
                               n.PC = 0,
                               min.MAF = min.MAF,
                               n.core = parallel::detectCores(),
                               P3D = P3D)

    gwa_results_pr <- gwa_results %>%
        dplyr::rename(log10p = trait) %>%
        dplyr::filter(log10p != 0) %>%
        dplyr::mutate(BF = log10(.05/n()),
                      trait = colnames(Y)[2]) %>%
        dplyr::select(CHROM, POS, marker, trait, BF, log10p)

    gwa_results_pr_emma <- gwa_results %>%
        dplyr::rename(log10p = trait) %>%
        dplyr::filter(log10p != 0) %>%
        dplyr::mutate(BF = log10(.05/n()),
                      trait = colnames(Y)[2]) %>%
        dplyr::select(CHROM, POS, marker, trait, BF, log10p)

    return(gwa_results_pr)

    # GWAS WITH SOMMER PACKAGE - SEEMS MUCH SLOWER THAN rrBLUP
    # # Remove NAs from phenotype set
    # Y = na.omit(phenotype)
    # colnames(Y) <- c("strain", "trait")
    #
    # # Sort kinship matrix to to be in the same order as phenotyped strains
    # kinship_sorted = kinship[sort(row.names(kinship)),sort(colnames(kinship))]
    # # Prune kinship matrix to only include phenotyped individuals
    # K = kinship_sorted[row.names(kinship_sorted) %in% Y$strain, colnames(kinship_sorted) %in% Y$strain]
    # # Sort kinship matrix to to be in the same order as phenotyped strains
    # markers = genotype %>%
    #     tidyr::unite(marker, CHROM, POS) %>%
    #     dplyr::select(marker, everything(), -REF, -ALT)
    #
    # markers_sorted <- data.frame(marker = markers$marker,
    #                              markers[,sort(names(markers[,2:ncol(markers)]))])
    #
    # row.names(markers_sorted) = markers_sorted$marker
    #
    # M <- markers_sorted %>%
    #     dplyr::select(-marker)
    #
    # M <- M[,colnames(M) %in% Y$strain]
    #
    # keepMarkers <- data.frame(
    #     MAF = apply(M,
    #                 MARGIN = 1,
    #                 FUN = function(x){
    #                     x[x==-1] <- 0
    #                     return(sum(x, na.rm = T)/length(x))}
    #     )) %>%
    #     dplyr::mutate(marker = row.names(.)) %>%
    #     dplyr::filter(MAF >= min.MAF)
    #
    # M <- M[row.names(M) %in% keepMarkers$marker, ] %>%
    #     t()
    #
    # Za <- model.matrix(~strain-1, Y)
    # ETA <- list(add=list(Z = Za,
    #                      K = K))

    # gwa_result <- sommer::GWAS(Y = Y$trait,
    #                            Z = ETA,
    #                            M = M[,1:100],
    #                            method = "EMMA",
    #                            draw = TRUE,
    #                            P3D = FALSE,
    #                            models="additive",
    #                            ploidy=2,
    #                            min.MAF=0.05,
    #                            gwas.plots = FALSE,
    #                            fdr.level=0.05)

    # p_values <- data.frame(t(gwa_result$M.scores$score)) %>%
    #     dplyr::mutate(marker = row.names(t(gwa_result$M.scores$score))) %>%
    #     dplyr::select(marker, log10p = T1) %>%
    #     dplyr::mutate(BF = -log10(0.05/n()),
    #                   FDR = gwa_result[["fdr"]])

}



