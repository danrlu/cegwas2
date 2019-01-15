testthat::context("tests/testthat/test_process_phenotypes.R")

df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE))

testthat::test_that("Check number of rows in test phenotype set", {
    testthat::expect_true(nrow(df) == 7469)
    testthat::expect_true(sum(is.na(df)) == 22398)
})


testthat::test_that("Test process_phenotypes output with BAMF remove outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "BAMF",
                                        remove_outliers = TRUE)
    testthat::expect_equal(nrow(pr_phenotypes), 240)
    testthat::expect_equal(ncol(pr_phenotypes), 5)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 4)
})


testthat::test_that("Test process_phenotypes output with BAMF keep outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "BAMF",
                                        remove_outliers = FALSE)
    testthat::expect_equal(nrow(pr_phenotypes), 243)
    testthat::expect_equal(ncol(pr_phenotypes), 6)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 12)
    testthat::expect_equal(pr_phenotypes$strain[duplicated(pr_phenotypes$strain)],
                 c("JU2862", "LSJ1", "QG2075"))
    testthat::expect_equal(nrow(dplyr::filter(pr_phenotypes, outlier)), 3)
})

testthat::test_that("Test process_phenotypes output with Z keep outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "Z",
                                        remove_outliers = TRUE,
                                        threshold = 3) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., dplyr::any_vars(is.na(.)) )


    testthat::expect_equal(nrow(pr_phenotypes), 9)
    testthat::expect_equal(ncol(pr_phenotypes), 5)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 14)
    testthat::expect_equal(pr_phenotypes$strain,
                 c("CB4853", "JU1400", "JU2581", "JU2862", "JU751",
                   "LSJ1", "MY518", "NIC514", "QG2075"))
})


testthat::test_that("Test process_phenotypes output with TUKEY keep outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "TUKEY",
                                        remove_outliers = FALSE,
                                        threshold = 2) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., dplyr::any_vars(is.na(.)) )


    testthat::expect_equal(nrow(pr_phenotypes), 30)
    testthat::expect_equal(ncol(pr_phenotypes), 6)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 60)
})

testthat::test_that("Test process_phenotypes output with TUKEY keep outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "TUKEY",
                                        remove_outliers = FALSE,
                                        threshold = 2)


    testthat::expect_equal(nrow(pr_phenotypes), 255)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 60)
})


testthat::test_that("Test process_phenotypes output with MAD keep outliers", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "mean",
                                        prune_method = "MAD",
                                        remove_outliers = FALSE,
                                        threshold = 2) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., dplyr::any_vars(is.na(.)) )


    testthat::expect_equal(nrow(pr_phenotypes), 114)
    testthat::expect_equal(ncol(pr_phenotypes), 6)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 228)
})

testthat::test_that("Test process_phenotypes output with MAD keep outliers - median", {
    pr_phenotypes <- cegwas2::process_phenotypes(df = df,
                                        summarize_replicates = "median",
                                        prune_method = "MAD",
                                        remove_outliers = FALSE,
                                        threshold = 2) %>%
        dplyr::rowwise() %>%
        dplyr::filter_all( ., dplyr::any_vars(is.na(.)) )


    testthat::expect_equal(nrow(pr_phenotypes), 104)
    testthat::expect_equal(ncol(pr_phenotypes), 6)
    testthat::expect_equal(sum(is.na(pr_phenotypes)), 208)
})



