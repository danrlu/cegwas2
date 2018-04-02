# cegwas-2

[![Travis-CI Build Status](https://travis-ci.org/AndersenLab/cegwas2.svg?branch=master)](https://travis-ci.org/AndersenLab/cegwas2)
[![codecov](https://codecov.io/gh/AndersenLab/cegwas2/branch/master/graph/badge.svg)](https://codecov.io/gh/AndersenLab/cegwas2)

### Initial Setup

This documentation can be removed later.

```
devtools::create("cegwas2")
```

### Coverage

![coverage](https://codecov.io/gh/AndersenLab/cegwas2/graphs/tree.svg)


### Testing

```
testthat::auto_test_package()
```


### Linting

Code linting is performed using [lintr](https://github.com/jimhester/lintr). 

Linting locally

```
lintr::lint_package()
```

In some cases, you may need to exclude a code chunk from linting.
This should be rare, but you can do it with the following:

```
# Begin Exclude Linting"

code_I_want_to_not_lint <- function(x) {
    x + 1
}

# End Exclude Linting"
```

```
# Example usage
# Load example phenotype set
df <- data.table::fread(system.file("extdata",
                                    "test_phenotype.tsv",
                                    package = "cegwas2",
                                    mustWork = TRUE)) %>%
    dplyr::select(strain, trait1)

# Initialize ceGWAS R6 object
test_trait <- ceGWAS$new(phenotype = df)
# Set genotype and kinship matices
test_trait$set_markers(genotype_matrix = cegwas2::snps,
                       kinship_matrix = cegwas2::kinship)                  
# EMMAx mapping
test_trait$run_mapping(P3D = TRUE)
# EMMA mapping
test_trait$run_mapping(P3D = FALSE)
```

