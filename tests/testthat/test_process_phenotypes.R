# A test VCF file with no annotation field
df <- system.file("extdata",
                  "test_phenotype.tsv",
                  package = "cegwas2",
                  mustWork = TRUE)

test_that("Check number of rows in test phenotype set", {
    expect_true(nrow(df) == 7469)
})
