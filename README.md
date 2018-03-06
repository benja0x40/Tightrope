Tightrope
================================================================================
**Normalization of ChIP-seq experiments involving global variations of chromatin marks**

### Overview

`Tightrope` is an `R` package proposing a ChIP-seq normalization method, named
Background Read Density (BRD)[<sup>1,2</sup>](#1), capable of accurate
estimation of normalization factors in conditions involving global variations
of chromatin marks.

Documentation and source code examples can be found in the package
[vignette]("./vignettes/BRD.pdf") and reference manual.

#### Main features

* Generation of read count matrixes from mapped reads (bam files)
* Estimation of ChIP-seq normalization factors with the BRD method
* Plot functions to control the consistency of normalization factors

### <a name="install"></a>Installation

Run the `R` code below to install `Tightrope`.

```R
library("devtools")
# Plotting and color mapping functions
install_github("benja0x40/Barbouille", dependencies = T)
# BRD implementation and related functions
install_github("benja0x40/Tightrope", dependencies = T)
```

If the installation fails, try to install dependencies as indicated
in the following section.

#### Dependencies

  - [R environment](https://www.r-project.org/) version 3.4 or newer
  - To develop and execute `R` scripts we recommend using [RStudio](https://www.rstudio.com/products/rstudio/download)
  - CRAN packages `devtools`, `stringr`, `triangle`, `caTools`, `ica`, `FNN`, `igraph`, `mixtools`, `data.table`
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
lst <- c("devtools", "stringr", "triangle", "caTools", "ica", "FNN", "igraph", "mixtools", "data.table")
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

# Acknowledgements

Thanks to Itys Comet, Daria Shlyueva, Aliaksandra Radzisheuskaya,
Sachin Pundhir and Albin Sandelin for their comments and suggestions
during the development of the BRD method,
and to Jens Vilstrup Johansen and Sudeep Sahadevan
from the bioinformatics core facility at the
[Biotech Research and Innovation Centre](http://www.bric.ku.dk)
for their support.

### References

<a name="1"></a>1. Mohammad et al. 2017 - *EZH2 is a potential therapeutic target for H3K27M-mutant pediatric gliomas.*  
[publisher](https://dx.doi.org/10.1038/nm.4293) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28263309)

<a name="2"></a>2. Leblanc B., Mohammad F., Hojfeldt J., Helin K. - *Normalization of ChIP-seq experiments involving global variations of chromatin marks.*  
(in preparation)
