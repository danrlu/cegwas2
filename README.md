# cegwas-2

[![Travis-CI Build Status](https://travis-ci.org/AndersenLab/cegwas2.svg?branch=master)](https://travis-ci.org/AndersenLab/cegwas2)
[![codecov](https://codecov.io/gh/AndersenLab/cegwas2/branch/master/graph/badge.svg)](https://codecov.io/gh/AndersenLab/cegwas2)

### Initial Setup

This documentation can be removed later.

```
devtools::create("cegwas2")
```

### Coverage

![coverage](https://codecov.io/gh/AndersenLab/cegwas2/commit/111aa38cd0f4a5010ef4334ea29761e83aeac3b6/graphs/tree.svg)


### Testing

```r
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

