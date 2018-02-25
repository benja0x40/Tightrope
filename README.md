Tightrope
================================================================================
**Normalization of ChIP-seq experiments involving global variations of chromatin marks**

### Overview

`Tightrope` is an R package for the normalization and representation
of ChIP-seq enrichments relative to the experimental background, which
provides an empirical solution to the challenging task of estimating
normalization factors in conditions involving global variations
of chromatin marks.

Accurate normalization of such experiments is increasingly needed in the field
of chromatin genomics. Knowing that basic sequencing depth correction
can be dramatically erroneous, this situation led to the recent
introduction of spike-in protocols[<sup>1-4</sup>](#1) proposing
an experimental solution to ChIP-seq normalization issues.
Seeking an alternative strategy due to reliability concerns with spike-in 
results while we were studying H3K27me3 profiles in a cellular model
of DIPG[<sup>5</sup>](#5), we designed a naive computational approach named
Background Read Density (BRD).
After our initial validation and successful application of the BRD normalization
for H3K27me3, we introduced several key improvements to achieve accurate
normalizations with various experimental setups and chromatin marks.
We are currently working on the manuscript describing the BRD method
and presenting comparative analyses between spike-in
and BRD normalizations for different datasets[<sup>6</sup>](#6).

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
install_github("benja0x40/Barbouille") # Graphic and color mapping functions
install_github("benja0x40/Tightrope")
```

If the installation fails, try to install dependencies as indicated
in the following section.

#### Dependencies

  - [R environment](https://www.r-project.org/) version 3.4 or newer
  - To develop and execute `R` scripts we recommend using [RStudio](https://www.rstudio.com/products/rstudio/download)
  - CRAN packages `devtools`, `triangle`, `caTools`, `ica`, `FNN`, `igraph`, `mixtools`
  - [Bioconductor](http://www.bioconductor.org/) packages `Rsamtools`, `GenomicAlignments`, `GenomicRanges`
  - GitHub package [Barbouille](https://github.com/benja0x40/Barbouille)

After installing the R environment (and RStudio), run the `R` code below
to install all dependencies as well as `Tightrope`.

```R
# Setting value below to TRUE will reinstall all required packages (optional)
reinstall <- FALSE

# Detect already installed packages
pkg <- ifelse(reinstall, c(), installed.packages()[, "Package"])

# CRAN packages
lst <- c("devtools", "triangle", "caTools", "ica", "FNN", "igraph", "mixtools")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) install.packages(lst, dependencies = T)

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
lst <- c("Rsamtools", "GenomicAlignments", "GenomicRanges")
lst <- setdiff(lst, pkg)
if(length(lst) > 0) biocLite(lst)

# GitHub packages
library("devtools")
install_github("benja0x40/Barbouille")
install_github("benja0x40/Tightrope")
```

### Acknowledgements

Thanks to Itys Comet, Sachin Pundhir and Albin Sandelin for their comments and
suggestions during the development of the BRD normalization, and to Sudeep
Sahadevan and Jens Vilstrup Johansen from the bioinformatics core facility
at the [Biotech Research and Innovation Centre](http://www.bric.ku.dk). 

### References

<a name="1"></a>1. Bonhoure et al. 2014 - *Quantifying ChIP-seq data: a spiking method providing an internal reference for sample-to-sample normalization*  
[publisher](https://dx.doi.org/10.1101/gr.168260.113) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/24709819)

<a name="2"></a>2. Orlando et al. 2014 - *Quantitative ChIP-Seq normalization reveals global modulation of the epigenome.*  
[publisher](https://dx.doi.org/10.1016/j.celrep.2014.10.018) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/25437568)

<a name="3"></a>3. Hu et al. 2015 - *Biological chromodynamics: a general method for measuring protein occupancy across the genome by calibrating ChIP-seq.*  
[publisher](https://dx.doi.org/10.1093/nar/gkv670) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/26130708)

<a name="4"></a>4. Egan et al. 2016 - *An Alternative Approach to ChIP-Seq Normalization Enables Detection of Genome-Wide Changes in Histone H3 Lysine 27 Trimethylation upon EZH2 Inhibition.*  
[publisher](https://dx.doi.org/10.1371/journal.pone.0166438) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27875550)

<a name="5"></a>5. Mohammad et al. 2017 - *EZH2 is a potential therapeutic target for H3K27M-mutant pediatric gliomas.*  
[publisher](https://dx.doi.org/10.1038/nm.4293) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28263309)

Leblanc B., Mohammad F., Hojfeldt J., Helin K. - *Normalization of ChIP-seq experiments involving global variations of chromatin marks.*  
(in preparation)
