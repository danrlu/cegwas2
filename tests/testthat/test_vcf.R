
# A test VCF file with no annotation field
test_vcf_file <- system.file("extdata",
                             "test.vcf.gz",
                             package = "cegwas2",
                             mustWork = TRUE)

test_that("Show info is null", {
    expect_true(is.null(query_vcf(vcf = test_vcf_file)))
})


test_that("Test that a genotype matches what we expect", {
    expect_equal(
        nrow(query_vcf("I:1-10000","II:1-100000") %>%
        dplyr::filter(SAMPLE == "JU311",
                      CHROM == "I",
                      POS == 811,
                      a1 == "G",
                      a2 == "G")),
        1)
})

test_that("Impact filtering for low and high", {
    expect_true(
        all((query_vcf("MtDNA:1-100000",
                   impact = c("LOW", "HIGH")) %>%
             dplyr::select(impact) %>%
             dplyr::distinct() %>%
             dplyr::arrange(impact))$impact == c("HIGH", "LOW"))
    )
})
