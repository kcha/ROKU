<!-- README.md is generated from README.Rmd. Please edit that file -->



ROKU
====

An R implementation of ROKU, a method for identifying tissue-specific genes.

Details of the method are described in:

Kadota, K., Ye, J., Nakai, Y., Terada, T., Shimizu, K., 2006. [ROKU: a novel method for identification of tissue-specific genes](http://www.biomedcentral.com/1471-2105/7/294). BMC Bioinformatics 7, 294.

Dependencies
------------

ROKU uses the `affy` package, which can be installed from Bioconductor:

``` {.r}
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
```

Usage
-----

ROKU requires as input a data frame of expression values where each row is a gene and each column is a tissue/sample. For example:

``` {.r}
head(psidata)
```

To find tissue-specific genes:

``` {.r}
rk <- ROKU(psidata)
df <- ROKU_as_df(rk)
head(df)
```

Genes with a non-zero U statistic are considered outliers.
