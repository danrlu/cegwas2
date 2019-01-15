testthat::context("tests/testthat/test_object.R")

options(stringsAsFactors = FALSE)

df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE)) %>%
    dplyr::select(strain, trait1)


test_trait <- ceGWAS$new(phenotype = df)

test_that("Test ceGWAS doesn't change input data", {

    remakedf <- data.frame(strain = as.character(test_trait$strains),
                           trait1 = test_trait$phenotype)
    dftest <- data.frame(na.omit(df))

    expect_true(identical(dftest,remakedf))

})

test_trait$set_markers(genotype_matrix = cegwas2::snps,
                       kinship_matrix = cegwas2::kinship)
test_trait$run_mapping(P3D = TRUE)

blup_message <- capture.output(gmap <- perform_mapping(phenotype = test_trait$processed_phenotype,
                                                       P3D = TRUE, min.MAF = 0.1, mapping_cores = 1))

test_that("Test ceGWAS mappings", {

    expect_equal(colnames(test_trait$mapping), colnames(gmap))
})
