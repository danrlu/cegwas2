

test_that("Test wormbase_gene table", {
    expect_true(
            any(
                get_db("wormbase_gene") %>%
                dplyr::filter(chrom == 'MtDNA') %>%
                dplyr::collect() %>%
                dplyr::select(locus) %>%
                unique() %>%
                dplyr::pull(locus) == "nduo-6"
            )
    )
})


test_that("Test strain table", {
    expect_true(
            any(
                get_db("strain") %>%
                    dplyr::filter(isotype == 'N2') %>%
                    dplyr::collect() %>%
                    dplyr::pull(previous_names) == "GA1"
            )
        )
})


test_that("Test homolog table", {
    expect_true(
        any(
            get_db("homologs") %>%
                dplyr::filter(gene_name == 'acdh-7') %>%
                dplyr::collect() %>%
                dplyr::pull(homolog_gene) == "ACADM"
        )
    )
})
