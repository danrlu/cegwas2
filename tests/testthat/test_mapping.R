testthat::context("tests/testthat/test_mapping.R")

df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE))


test_that("Test EMMAx mapping", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                                 summarize_replicates = "mean",
                                                 prune_method = "BAMF",
                                                 remove_outliers = TRUE)
    gmap <- cegwas2::perform_mapping(phenotype = pr_phenotypes[20:240,c(1,2)],
                            P3D = TRUE,
                            min.MAF = 0.1)

    expect_true(min(gmap$qvalue) < 0.05)
})


test_that("Test EMMA mapping with subset of strains for speed", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                                 summarize_replicates = "mean",
                                                 prune_method = "BAMF",
                                                 remove_outliers = TRUE)
    gmap <- cegwas2::perform_mapping(phenotype = pr_phenotypes[20:110,c(1,2)],
                                     min.MAF = 0.1)

    expect_false(min(gmap$qvalue) < 0.05)
})


test_that("Test SNV matrix format check ", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                                 summarize_replicates = "mean",
                                                 prune_method = "BAMF",
                                                 remove_outliers = TRUE)
    expect_error(cegwas2::perform_mapping(phenotype = pr_phenotypes[20:110,c(1,2)],
                                          genotype = cegwas2::snps[,4:ncol(cegwas2::snps)]))
})


test_that("Test kinship matrix format check", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                                 summarize_replicates = "mean",
                                                 prune_method = "BAMF",
                                                 remove_outliers = TRUE)
    expect_error(cegwas2::perform_mapping(phenotype = pr_phenotypes[20:110,c(1,2)],
                                          kinship = cegwas2::kinship[,1:239]))
})
