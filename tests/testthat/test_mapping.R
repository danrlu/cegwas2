df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE))


test_that("Test EMMAx mapping", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                                 summarize_replicates = "mean",
                                                 prune_method = "BAMF",
                                                 remove_outliers = TRUE)

    gmap <- cegwas2::perform_mapping(phenotype = pr_phenotypes[,1:2],
                            genotype = cegwas2::snps,
                            kinship = cegwas2::kinship,
                            P3D = TRUE,
                            min.MAF = 0.05)
})




