# A test VCF file with no annotation field
df <- system.file("extdata",
                  "test_phenotype.tsv",
                  package = "cegwas2",
                  mustWork = TRUE)

test_that("Check number of rows in test phenotype set", {
    expect_true(nrow(df) == 7469)
})



test_that("Test process_phenotypes output with BAMF remove outliers", {
    pr_phenotypes <- process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "BAMF",
                                        remove_outliers = TRUE)
    expect_equal(nrow(pr_phenotypes), 240)
    expect_equal(ncol(pr_phenotypes), 5)
    expect_equal(sum(is.na(pr_phenotypes)), 4)
})


test_that("Test process_phenotypes output with BAMF keep outliers", {
    pr_phenotypes <- process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "BAMF",
                                        remove_outliers = FALSE)
    expect_equal(nrow(pr_phenotypes), 243)
    expect_equal(ncol(pr_phenotypes), 6)
    expect_equal(sum(is.na(pr_phenotypes)), 12)
    expect_equal(pr_phenotypes$strain[duplicated(pr_phenotypes$strain)], c("JU2862", "LSJ1", "QG2075"))
    expect_equal(nrow(dplyr::filter(pr_phenotypes, outlier)), 3)
})

test_that("Test process_phenotypes output with Z keep outliers", {
    pr_phenotypes <- process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "Z",
                                        remove_outliers = TRUE) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., any_vars(is.na(.)) )


    expect_equal(nrow(pr_phenotypes), 9)
    expect_equal(ncol(pr_phenotypes), 5)
    expect_equal(sum(is.na(pr_phenotypes)), 14)
    expect_equal(pr_phenotypes$strain, c("CB4853", "JU1400", "JU2581", "JU2862", "JU751", "LSJ1", "MY518", "NIC514", "QG2075"))
})


test_that("Test process_phenotypes output with Z keep outliers", {
    pr_phenotypes <- process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "TUKEY",
                                        remove_outliers = FALSE) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., any_vars(is.na(.)) )


    expect_equal(nrow(pr_phenotypes), 9)
    expect_equal(ncol(pr_phenotypes), 5)
    expect_equal(sum(is.na(pr_phenotypes)), 14)
    expect_equal(pr_phenotypes$strain, c("CB4853", "JU1400", "JU2581", "JU2862", "JU751", "LSJ1", "MY518", "NIC514", "QG2075"))
})


