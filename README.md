Tightrope
================================================================================

## Overview

`Tightrope` is an R package for the normalization and analysis
of ChIP-seq data on histone marks and histone variants.

## Examples

Here are some examples...

Further documentation can be found in the package vignette and reference manual.

## Installation

Run the `R` code below to install `Tightrope`.

```R
library("devtools")
install_github("benja0x40/Tightrope")
```

If the installation fails, try to install dependencies as indicated below.

### Dependencies

  - [R environment](https://www.r-project.org/) version 3.4 or newer
  - CRAN packages `devtools`, `igraph`, `FNN`, `triangle`, `mixtools`
  - [Bioconductor](http://www.bioconductor.org/) packages `GenomicAlignments`, `GenomicRanges`
  - GitHub `R` package
    [Barbouille](https://github.com/benja0x40/Barbouille)

Run the `R` code below to install all dependencies.

```R
# Setting value below to TRUE will reinstall all required packages (optional)
reinstall <- FALSE

# Detect already installed packages
pkg <- ifelse(reinstall, c(), installed.packages()[, "Package"])

# CRAN packages
lst <- c("devtools", "igraph", "FNN", "triangle", "mixtools")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) install.packages(lst, dependencies = T)

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
lst <- c("GenomicAlignments", "GenomicRanges")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) biocLite(lst)

# GitHub package
library("devtools")
lst <- paste0("benja0x40/", c("Barbouille", "Tightrope"))
lst <- setdiff(lst, pkg)
if(length(lst) > 0) lapply(lst, install_github)
```

## Contact

Benjamin Leblanc (benjaminolivierleblanc@gmail.com)

## References

Perec 1980 - *Experimental demonstration of the tomatotopic organization in the Soprano (Cantatrix sopranica L.).*  
[publisher](http://dx.doi.org/10.2307/3684039) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/)

