
test_vcf_file <- system.file("extdata",
                             "test.vcf.gz",
                             package = "cegwas2",
                             mustWork = TRUE)


test_that("str_length is number of characters", {
    expect_equal(str_length("a"), 1)
    expect_equal(str_length("ab"), 2)
    expect_equal(str_length("abc"), 3)
})
