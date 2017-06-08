Tightrope
================================================================================

Normalization and analysis of ChIP-seq on variant histones and histone marks.

### Package installation ###

#### Prerequisites ####

  - [R environment](https://www.r-project.org/) version 3.x
  - CRAN packages `devtools`, `igraph`, `FNN`, `triangle`, `mixtools`
  - [Bioconductor](http://www.bioconductor.org/) packages
    `GenomicAlignments`, `GenomicRanges`
  - [GitHub](https://github.com/benja0x40/) package `Barbouille`
  
Run the R code below to install CRAN and Bioconductor packages dependencies
for Tightrope.

```R
# Already installed
pkg <- installed.packages()[, "Package"]

# CRAN packages
lst <- c(
"devtools", "igraph", "FNN", "triangle", "mixtools"
)
lst <- setdiff(lst, pkg)
if(length(lst) > 0) {
  install.packages(lst, repos = "https://cloud.r-project.org/")
}

# Bioconductor packages
lst <- c("GenomicAlignments", "GenomicRanges")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(lst)
}
```

#### Installations from github ####

Run the bash code below to build package Tightrope from github.

```bash
# Clone github repository
cd ~/DataImportTools
git clone git@github.com:benja0x40/Tightrope.git

# Update cloned repository
cd ~/DataImportTools/Tightrope
git pull

# Build package
cd ..
R CMD build Tightrope
```
Run the R code below to install Tightrope.

```r
# When package will be public
# library("devtools")
# install_github("benja0x40/Tightrope")

# Using manually built package archive
install.packages("Tightrope_0.1.0.tar.gz")
```

