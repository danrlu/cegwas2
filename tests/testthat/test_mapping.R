testthat::context("tests/testthat/test_mapping.R")

df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE))

pr_phenotypes <- process_phenotypes(df = df,
                                    summarize_replicates = "mean",
                                    prune_method = "BAMF",
                                    remove_outliers = TRUE)

to_map = pr_phenotypes[20:240,c(1,2)]

test_that("Test EMMAx mapping", {
    gmap <- perform_mapping(phenotype = to_map, P3D = TRUE, min.MAF = 0.1, mapping_cores = 1)

    expect_true(min(gmap$qvalue) < 0.05)
})

to_map1 <- pr_phenotypes[20:110,c(1,2)]

test_that("Test EMMA mapping with subset of strains for speed", {

    gmap <- perform_mapping(phenotype = to_map1, min.MAF = 0.1, mapping_cores = 1)

    expect_false(min(gmap$qvalue) < 0.05)
})

bad_genotype <- cegwas2::snps[,4:ncol(cegwas2::snps)]

test_that("Test SNV matrix format check ", {
    expect_error(cegwas2::perform_mapping(phenotype = to_map1, genotype = bad_genotype))
})

bad_kinship <- cegwas2::kinship[,1:239]

test_that("Test kinship matrix format check", {
    expect_error(cegwas2::perform_mapping(phenotype = to_map1, kinship = bad_kinship))
})
