<!-- README.md is generated from README.Rmd. Please edit that file -->
ROKU
====

An R implementation of ROKU, a method for identifying tissue-specific genes.

Details of the method are described in:

Kadota, K., Ye, J., Nakai, Y., Terada, T., Shimizu, K., 2006. [ROKU: a novel method for identification of tissue-specific genes](http://www.biomedcentral.com/1471-2105/7/294). BMC Bioinformatics 7, 294.

**Disclaimer**

This package was written as an personal exercise to practice implementing a published algorithm in R, and to learn now to create an R package. While some effort has been made to ensure accuracy of the algorithm, the package has not been thoroughly tested and documented. Use at your own risk (see [LICENSE](LICENSE)).

Dependencies
------------

ROKU uses the `affy` package, which can be installed from Bioconductor:

``` r
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
```

Installation
------------

``` r
install.packages("devtools")
devtools::install_github("kcha/ROKU")
```

Usage
-----

ROKU requires as input a data frame of expression values where each row is a gene and each column is a tissue/sample. For example:

``` r
head(psidata)
#>      Neural1    Neural2   Neural3      ESC1       ESC2     Kidney
#> A  27.125019  27.106195  25.32516  94.19426  98.596466  91.348439
#> B   2.032126   3.627567   3.74161  68.60577  69.543868   1.924062
#> C  96.864982  97.846862  97.00417  97.23877  97.383610   0.000000
#> D 100.000000 100.000000 100.00000 100.00000 100.000000 100.000000
#> E  56.707902  43.974580  49.61398  39.20734  64.567588  58.050368
#> F   9.867039   5.026210   5.85521   5.16909   6.879749  83.756046
#>         Liver    Testis     Muscle     Heart
#> A  99.7869522  96.74810  98.433381  92.74806
#> B   0.7742677   6.59369   3.072821   6.44550
#> C  93.5602998  25.00000  91.780465  50.00000
#> D 100.0000000 100.00000 100.000000 100.00000
#> E  51.1538912  56.41729  56.192312  44.57487
#> F  90.2871259  91.50642  87.251575  94.61886
```

To find tissue-specific genes:

``` r
rk <- ROKU(psidata)
df <- ROKU_as_df(rk)
head(df)
#>    Entropy Entropy.Normalized Outlier.Detection.Method   Neural1
#> A 2.098470          0.6317024                      AIC -68.73251
#> B 1.564075          0.4708333                      AIC   0.00000
#> C 1.902679          0.5727633                      AIC   0.00000
#> D       NA                 NA                     <NA>   0.00000
#> E 3.059714          0.9210658                      AIC   0.00000
#> F 3.318326          0.9989157                      AIC   0.00000
#>      Neural2    Neural3      ESC1     ESC2     Kidney     Liver    Testis
#> A -68.751338 -70.532370   0.00000  0.00000   0.000000  0.000000   0.00000
#> B   0.000000   0.000000  65.11057 66.04866   0.000000  0.000000   0.00000
#> C   0.000000   0.000000   0.00000  0.00000 -96.003539 -2.443239 -71.00354
#> D   0.000000   0.000000   0.00000  0.00000   0.000000  0.000000   0.00000
#> E  -9.128835  -3.489432 -13.89607 11.46417   4.946952 -1.949524   0.00000
#> F   0.000000   0.000000   0.00000  0.00000  35.969799 42.500879  43.72017
#>      Muscle      Heart
#> A  0.000000   0.000000
#> B  0.000000   0.000000
#> C -4.223074 -46.003539
#> D  0.000000   0.000000
#> E  0.000000  -8.528546
#> F 39.465328  46.832609
```

Genes with a non-zero values are considered outliers.
