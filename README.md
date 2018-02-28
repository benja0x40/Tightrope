Tightrope
================================================================================
**Normalization of ChIP-seq experiments involving global variations of chromatin marks**

### Overview

`Tightrope` is an R package for the normalization and representation
of ChIP-seq enrichments relative to the experimental background, which
provides an empirical solution to the challenging task of estimating
normalization factors in conditions involving global variations
of chromatin marks.

`Tightrope` is an `R` package proposing a ChIP-seq normalization method, named
Background Read Density (BRD)[<sup>1,2</sup>](#1), capable of accurate
estimation of normalization factors in conditions involving global variations
of chromatin marks.

Documentation and source code examples showing how to use `Tightrope` 
to normalize a ChIP-seq dataset with the BRD method can be found
in the package vignettes and the reference manual.

#### Rationale of the BRD normalization

As a prerequisite, the ChIP-seq dataset must include background profiles
(e.g. input DNA, ChIP-seq performed with IgG or in KO conditions), preferably
generated using the same sequencing setup as for the experimental profiles.

The BRD normalization is build upon the two following assumptions.

1. For all considered experimental conditions, some portions of the genome are
   invariably devoid of the immuno-precipitated chromatin mark.
   These invariable genomic portions are named **background candidates**.
   
2. Background candidates can be distinguished as a density mode of
   transformed ChIP-seq read count distributions

#### Main features

* Generation of read count matrixes from mapped reads (bam files)
* Estimation of ChIP-seq normalization factors with the BRD method[<sup>5, 6</sup>](#5)
* Density plots for read count distributions

### <a name="install"></a>Installation

Run the `R` code below to install `Tightrope`.

```R
library("devtools")
# Graphic and color mapping functions
install_github("benja0x40/Barbouille", dependencies = T)
install_github("benja0x40/Tightrope", dependencies = T)
```

If the installation fails, try to install dependencies as indicated
in the following section.

#### Dependencies

  - [R environment](https://www.r-project.org/) version 3.4 or newer
  - To develop and execute `R` scripts we recommend using [RStudio](https://www.rstudio.com/products/rstudio/download)
  - CRAN packages `devtools`, `stringr`, `triangle`, `caTools`, `ica`, `FNN`, `igraph`, `mixtools`
  - [Bioconductor](http://www.bioconductor.org/) packages `Rsamtools`, `GenomicAlignments`, `GenomicRanges`, `GenomicFeatures`, `rtracklayer`, `biomaRt`
  - GitHub package [Barbouille](https://github.com/benja0x40/Barbouille)

After installing the R environment (and RStudio), run the `R` code below
to install all dependencies as well as `Tightrope`.

```R
# Setting value below to TRUE will reinstall all required packages (optional)
reinstall <- FALSE

# Detect already installed packages
pkg <- ifelse(reinstall, c(), installed.packages()[, "Package"])

# CRAN packages
lst <- c("devtools", "stringr", "triangle", "caTools", "ica", "FNN", "igraph", "mixtools")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) install.packages(lst, dependencies = T)

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
lst <- c("Rsamtools", "GenomicAlignments", "GenomicRanges", "GenomicFeatures", "rtracklayer", "biomaRt")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) biocLite(lst)

# GitHub packages
library("devtools")
install_github("benja0x40/Barbouille", dependencies = T)
install_github("benja0x40/Tightrope")
```

### Acknowledgements

Thanks to Itys Comet, Sachin Pundhir and Albin Sandelin for their comments and
suggestions during the development of the BRD normalization, and to Sudeep
Sahadevan and Jens Vilstrup Johansen from the bioinformatics core facility
at the [Biotech Research and Innovation Centre](http://www.bric.ku.dk). 
