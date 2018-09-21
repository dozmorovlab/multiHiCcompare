# multiHiCcompare

## Overview

`multiHiCcompare` is an extension of the original `HiCcompare` R package (http://bioconductor.org/packages/HiCcompare/ or https://github.com/dozmorovlab/HiCcompare). `multiHiCcompare` provides functions for the joint normalization and comparison of complex Hi-C experiments. `multiHiCcompare` operates on processed Hi-C data in the form of sparse upper triangular matrices. 

`multiHiCcompare` accepts four-column text files storing chromatin interaction matrices in a sparse matrix format. There are many sources of public Hi-C data such as http://aidenlab.org/data.html and http://cooler.readthedocs.io/en/latest/index.html. `multiHiCcompare` is designed to the give the user the ability to perform a comparative analysis on the 3D structure of the genomes of cells in different biological states or under different experimental conditions. `multiHiCcompare` implements a cyclic loess joint normalization algorithm to remove bias between multiple Hi-C datasets and prepare them for comparison. `multiHiCcompare` also provides a general linear model based difference detection method that makes use of the `edgeR` R package's framework. 

The main functions of `multiHiCcompare` are:

- `cyclic_loess()` and `fastlo()` for performing joint normalization.
- `hic_exactTest()` and `hic_glm()` for performing comparisons between experimental conditions.

Some example Hi-C data along with an extensive vignette documenting the usage of `multiHiCcompare` are included in the package. 

Read the full paper describing the methods used in `multiHiCcompare` here: <insert link>.

## Installation

First make sure you have all dependencies installed in R.

```
install.packages(c('dplyr', 'data.table', 'devtools', 'qqman', 'metap', 'pheatmap', 'pbapply'))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BiocParallel", "HiCcompare", "edgeR", "GenomicRanges", "GenomeInfoDbData"))
```

To install HiCcompare from bioconductor open R and enter the following commands. Currently it is recommended to use the GitHub release or the development version of the bioconductor release.

```
# It is recommended to use the github release until the next Bioconductor update
# BiocManager::install("multiHiCcompare")
# library(multiHiCcompare)
```

Or to install the latest version of HiCcompare directly from the github release open R and enter the following commands.

```
library(devtools)
install_github('jstansfield0/multiHiCcompare', build_vignettes = TRUE)
library(multiHiCcompare)
```

## Usage

To use `multiHiCcompare` you will first need to obtain some Hi-C data. Data is available from the sources listed in the overview along with many others. You will need to extract the data and read it into R as either a 3 column sparse upper triangular matrix and then combine it with an additional column for the chromosome. Hi-C data ready to be used in `multiHiCcompare` should look like the following:

```
  chr  region1  region2 IF
1  22 16000000 16000000 11
2  22 16100000 16100000  1
3  22 16200000 16200000  3
4  22 16300000 16300000 15
5  22 16400000 16400000  3
6  22 16400000 16500000  1
...
```

The four columns correspond to the chromosome, the start location of the first interacting region (in basepairs), the start location of the second interacting region, and the interaction frequency (IF) of the interaction. 

Below is an example analysis of a Hi-C experiment.

```
# load data
library(multiHiCcompare)
data("r1", "r2", "r3", "r4")

# make hicexp object
hicexp <- make_hicexp(r1, r2, r3, r4, groups = c(1,1,2,2))

# jointly normalize data
hicexp <- cyclic_loess(hicexp)

# compare groups
hicexp <- hic_exactTest(hicexp)

# view manhattan plot of results
manhattan_hicexp(hicexp)
```

Refer to the `multiHiCcompare` vignette for full usage instructions. For a full explanation of the methods used in `multiHiCcompare` see the manuscript here <link>.

To view the usage vignette:

browseVignettes("multiHiCcompare")

## Citation
Please cite `multiHiCcompare` if you use it in your analysis.

<link>

## Contributions & Support
Suggestions for new features and bug reports are welcome. Please create a new issue for any of these or contact the author directly: @jstansfield0 (stansfieldjc@vcu.edu)

## Contributors
Authors: @jstansfield0 (stansfieldjc@vcu.edu) & @mdozmorov (mikhail.dozmorov@vcuhealth.org)