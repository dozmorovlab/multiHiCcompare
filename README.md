# [multiHiCcompare](https://dozmorovlab.github.io/multiHiCcompare/)

## Overview

`multiHiCcompare` is an extension of the original [HiCcompare R package](http://bioconductor.org/packages/HiCcompare/). `multiHiCcompare` provides functions for the joint normalization and comparison of complex Hi-C experiments. `multiHiCcompare` operates on processed Hi-C data in the form of sparse upper triangular matrices.

`multiHiCcompare` accepts four-column text files storing chromatin interaction matrices in a sparse matrix format. There are many sources of public Hi-C data such as the [Aiden Lab](http://aidenlab.org/data.html) (`.hic` files) and the [Mirnylab FTP site](http://cooler.readthedocs.io/en/latest/index.html) (`.cool` files). `multiHiCcompare` performs differential chromatin interaction analysis between two biological conditions, one or multiple Hi-C matrices per condition.  

`multiHiCcompare` implements a cyclic loess joint normalization algorithm to remove bias between multiple Hi-C datasets and prepare them for comparison. `multiHiCcompare` also provides a general linear model-based difference detection method, implemented in the `edgeR` R package. 

The main functions of `multiHiCcompare` are:

- `cyclic_loess()` and `fastlo()` for performing joint normalization.
- `hic_exactTest()` and `hic_glm()` for performing comparisons between experimental conditions.

Some example Hi-C data are included in the package. Refer to the `multiHiCcompare` vignette for full usage instructions, `vignette("multiHiCcompare")`

## Installation

First, make sure you have all dependencies installed in R.

``` r
install.packages(c('dplyr', 'data.table', 'devtools', 'qqman', 'metap', 'pheatmap', 'pbapply'))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BiocParallel", "HiCcompare", "edgeR", "GenomicRanges", "GenomeInfoDbData"))
```

To install `multiHiCcompare` from Bioconductor, use the following commands.

``` r
# It is recommended to use the GitHub release until the next Bioconductor update
# BiocManager::install("multiHiCcompare")
# library(multiHiCcompare)
```

Or to install the latest version of `multiHiCcompare` directly from the GitHub.

``` r
library(devtools)
# Stable release
install_github('dozmorovlab/multiHiCcompare', build_vignettes = TRUE)
# Developmental version
# install_github('jstansfield0/multiHiCcompare', build_vignettes = TRUE)
library(multiHiCcompare)
```

## Usage

To use `multiHiCcompare`, you will first need to obtain some Hi-C data. Data is available from the sources listed in the overview, along with many others. You will need to extract the data and read it into R as either a 3 column sparse upper triangular matrix and then combine it with an additional column for the chromosome. Hi-C data ready to be used in `multiHiCcompare` should look like the following:

``` r
  chr  region1  region2 IF
1  22 16000000 16000000 11
2  22 16100000 16100000  1
3  22 16200000 16200000  3
4  22 16300000 16300000 15
5  22 16400000 16400000  3
6  22 16400000 16500000  1
...
```

The four columns correspond to the chromosome, the start location of the first interacting region (in base pairs), the start location of the second interacting region, and the interaction frequency (IF) of the interaction. 

Below is an example analysis of a Hi-C experiment.

``` r
# load data
library(multiHiCcompare)
data("r1", "r2", "r3", "r4")

# make hicexp object
hicexp <- make_hicexp(r1, r2, r3, r4, groups = c(1, 1, 2, 2))

# jointly normalize data
hicexp <- cyclic_loess(hicexp)

# compare groups
hicexp <- hic_exactTest(hicexp)

# view manhattan plot of results
manhattan_hicexp(hicexp)
```

## Citation

Stansfield, John C, Kellen G Cresswell, and Mikhail G Dozmorov. [MultiHiCcompare: Joint Normalization and Comparative Analysis of Complex Hi-C Experiments](https://doi.org/10.1093/bioinformatics/btz048). _Bioinformatics_, January 22, 2019. 

## Learn more

Stansfield, John C., Duc Tran, Tin Nguyen, and Mikhail G. Dozmorov. [R Tutorial: Detection of Differentially Interacting Chromatin Regions From Multiple Hi-C Datasets](https://doi.org/10.1002/cpbi.76). Current Protocols in Bioinformatics, May 2019

[HiCcompareWorkshop](https://github.com/mdozmorov/HiCcompareWorkshop) - "Detection of Differentially Interacting Chromatin Regions From Multiple Hi-C Datasets" workshop presented on Bioconductor 2020 conference


## Contributions & Support

[issue](https://github.com/dozmorovlab/HiCcompare/issues) for any of these or contact the author directly: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu)

## Contributors

Authors: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu) & [@mdozmorov](https://github.com/mdozmorov) (mikhail.dozmorov@vcuhealth.org)