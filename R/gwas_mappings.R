
perform_mapping <- function(phenotype = pr_phenotypes[,1:2],
                            genotype = cegwas2::snps,
                            kinship = cegwas2::kinship,
                            method = "EMMA",
                            model = "additive",
                            MAF,
                            FDR) {

    # Remove NAs from phenotype set
    Y = na.omit(phenotype)
    colnames(Y) <- c("strain", "trait")

    # Sort kinship matrix to to be in the same order as phenotyped strains
    kinship_sorted = kinship[sort(row.names(kinship)),sort(colnames(kinship))]
    # Prune kinship matrix to only include phenotyped individuals
    K = kinship_sorted[row.names(kinship_sorted) %in% Y$strain, colnames(kinship_sorted) %in% Y$strain]
    # Sort kinship matrix to to be in the same order as phenotyped strains
    markers = genotype %>%
        tidyr::unite(marker, CHROM, POS) %>%
        dplyr::select(marker, everything(), -REF, -ALT)

    markers_sorted <- data.frame(marker = markers$marker,
                                 markers[,sort(names(markers[,2:ncol(markers)]))])

    row.names(markers_sorted) = markers_sorted$marker

    M <- markers_sorted %>%
        dplyr::select(-marker) %>%
        t()

    M <- M[row.names(M) %in% Y$strain,]

    Za <- model.matrix(~strain-1, Y)
    ETA <- list(add=list(Z = Za,
                         K = K))

    gwa_result <- sommer::GWAS(Y = Y$trait,
                               Z = ETA,
                               M = M,
                               method = "EMMA",
                               draw = TRUE,
                               P3D = FALSE,
                               models="additive",
                               ploidy=2,
                               min.MAF=0.05,
                               gwas.plots = FALSE,
                               fdr.level=0.05)

}



