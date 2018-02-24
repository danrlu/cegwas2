
# A test VCF file with no annotation field
test_vcf_file <- system.file("extdata",
                             "test.vcf.gz",
                             package = "cegwas2",
                             mustWork = TRUE)

test_that("Show info is null", {
    expect_null(query_vcf())
})


test_that("Test that a genotype matches what we expect", {
    expect_equal(
        nrow(query_vcf("I:1-10000",
                       "II:1-100000",
                       format = c("TGT", "GT"),
                       samples = "CB4856") %>%
                 dplyr::filter(CHROM == "II",
                               genotype == 2,
                               POS == 96264))
        , 1)
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


test_that("Query gene ID", {
    expect_true(
        all((query_vcf("WBGene00014450"))$gene_id == "WBGene00014450")
    )
})


test_that("Query locus ID", {
    expect_true(
        all((query_vcf("pot-2"))$gene_name == "pot-2")
    )
})


test_that("Query chromosome", {
    expect_true(
        nrow(query_vcf("MtDNA")) > 10000
    )
})


test_that("Query specific sample", {
    expect_true(
        all(unique((query_vcf("V:1-10000", samples = "CB4856"))$SAMPLE) == "CB4856")
    )
})


test_that("Fetch INFO and FORMAT columns", {
    expect_true(
        all(
            c("DP", "AD") %in% names(query_vcf("I:1-10000", info = c("DP"), format = c("AD")))
        )
    )
})


test_that("Use impact = ALL", {
    expect_true(
        setequal(
        unique(
            query_vcf("I:1-100000",
                      impact = "ALL")$impact),
        c(NA, "MODIFIER", "LOW", "MODERATE", "HIGH")
        )
    )
})

#===============#
# Error queries #
#===============#

test_that("Invalid INFO column", {
    expect_error(
        query_vcf("pot-2",
                  info = c("Not a column")),
        regexp = "^VCF does not have specified info fields:"
        )
})


test_that("Invalid FORMAT column", {
    expect_error(
        query_vcf("pot-2",
                  format = c("Not a column")),
                  regexp = "^VCF does not have specified format fields:"
        )
})

test_that("Invalid IMPACT filter", {
    expect_error(
        query_vcf("pot-2",
                  impact = c("Not an impact")),
        regexp = "^Invalid impact")
})


test_that("Invalid sample specified", {
    expect_error(
        query_vcf("pot-2",
                  samples = c("CB4856", "Not a sample")),
        regexp = "^Specified.*",
    )
})
