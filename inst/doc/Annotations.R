## ----hidden_dependencies, include=FALSE----------------------------------
library(knitr)
suppressPackageStartupMessages(library(Tightrope))

## ----hidden_configuration, include=FALSE---------------------------------
CFG <- list(
  UpdatePackageData = F
)

## ----eval=FALSE----------------------------------------------------------
#  library(Tightrope)  # Load package
#  data("EGA91_human") # Load R data with gene features in Human
#  EGA91_human$summary # Summary of loaded gene features

## ----eval=FALSE----------------------------------------------------------
#  EGA91_human$GNU # Gene units

## ----echo=FALSE----------------------------------------------------------
# Simplify EGA91_human$GNU for cleaner/concise output in the pdf vignette
grg <- EGA91_human$GNU
mcols(grg) <- mcols(grg)[, 3, drop = F]
grg
suppressWarnings(rm(grg))

## ----echo=FALSE----------------------------------------------------------
kable(EGA91_human$summary, caption = "Summary of gene features in EGA91_human")

## ------------------------------------------------------------------------
# Human genome annotations (Ensembl GRCh38.p10 release 91)
data("EGA91_human")    # Gene features from Ensembl      (environment)
data("CGI_hg38")       # CpG Islands from UCSC           (GRanges)
data("hg38_blacklist") # Blacklisted regions from ENCODE (GRanges)

## ----echo=FALSE----------------------------------------------------------
kable(EGA91_mouse$summary, caption = "Summary of gene features in EGA91_mouse")

## ------------------------------------------------------------------------
# Mouse genome annotations (Ensembl GRCm38.p5 release 91)
data("EGA91_mouse")    # Gene features from Ensembl      (environment)
data("CGI_mm10")       # CpG Islands from UCSC           (GRanges)
data("mm10_blacklist") # Blacklisted regions from ENCODE (GRanges)

## ----eval=FALSE----------------------------------------------------------
#  library("devtools") # (devtools can be installed from CRAN repositories)
#  install_github("benja0x40/LittleThumb")

## ----include=FALSE-------------------------------------------------------
suppressPackageStartupMessages(library(LittleThumb))

## ----eval=FALSE----------------------------------------------------------
#  library(LittleThumb)

## ----eval=FALSE----------------------------------------------------------
#  # =============================================================================.
#  # Human genome annotations (Ensembl GRCh38.p10 release 91)
#  # -----------------------------------------------------------------------------.
#  organism <- list(
#    name    = "Homo sapiens",
#    taxid   = 9606,
#    UCSC    = list(genome = "hg38"),
#    Ensembl = list(release = 91),
#    BSg_pkg = "BSgenome.Hsapiens.UCSC.hg38",
#    TxDb = list(
#      source = "Ensembl GRCh38.p10 release 91 (December 2017)",
#      portal = "http://Dec2017.archive.ensembl.org/Homo_sapiens",
#      gtf = paste0(
#        "ftp://ftp.ensembl.org/pub/",
#        "release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz"
#      )
#    ),
#    BioMart = list(
#      host     = "www.ensembl.org", # Dec2017.archive.ensembl.org
#      database = "ensembl",
#      dataset  = "hsapiens_gene_ensembl"
#    )
#  )
#  
#  # Gene features from Ensembl
#  MakeObj(EGA91_human, {
#    EGA91_human <- BuildGeneFeatures(organism)
#  })
#  
#  # CpG Islands from UCSC
#  MakeObj(CGI_hg38, {
#    CGI_hg38 <- ImportCpGIslandExt(organism$UCSC$genome)
#  })

## ----include=FALSE-------------------------------------------------------
if(CFG$UpdatePackageData) {
  devtools::use_data(EGA91_human, overwrite = T)
  devtools::use_data(CGI_hg38, overwrite = T)
}

## ----eval=FALSE----------------------------------------------------------
#  # =============================================================================.
#  # Human genome annotations (Ensembl GRCh37.p13 release 75)
#  # -----------------------------------------------------------------------------.
#  organism <- list(
#    name  = "Homo sapiens",
#    taxid = 9606,
#    UCSC    = list(genome = "hg19"),
#    Ensembl = list(release = 75),
#    BSg_pkg = "BSgenome.Hsapiens.UCSC.hg19",
#    TxDb = list(
#      source = "Ensembl GRCh37.p13 release 75 (February 2014)",
#      portal = "http://Feb2014.archive.ensembl.org/Homo_sapiens",
#      gtf = paste0(
#        "ftp://ftp.ensembl.org/pub/",
#        "release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
#      )
#    ),
#    BioMart = list(
#      host     = "Feb2014.archive.ensembl.org", # "grch37.ensembl.org"
#      database = "ensembl",
#      dataset  = "hsapiens_gene_ensembl"
#    )
#  )
#  
#  # Gene features from Ensembl
#  MakeObj(EGA75_human, {
#    EGA75_human <- BuildGeneFeatures(organism)
#  })
#  
#  # CpG Islands from UCSC
#  MakeObj(CGI_hg19, {
#    CGI_hg19 <- ImportCpGIslandExt(organism$UCSC$genome)
#  })

## ----include=FALSE-------------------------------------------------------
if(CFG$UpdatePackageData) {
  devtools::use_data(EGA75_human, overwrite = T)
  devtools::use_data(CGI_hg19, overwrite = T)
}

## ----eval=FALSE, include=FALSE, echo=FALSE-------------------------------
#  kable(EGA75_human$summary, caption = "Summary of gene features in EGA75_human")

## ----eval=FALSE----------------------------------------------------------
#  # =============================================================================.
#  # Mouse genome annotations (Ensembl GRCm38.p5 release 91)
#  # -----------------------------------------------------------------------------.
#  organism <- list(
#    name    = "Mus musculus",
#    taxid   = 10090,
#    UCSC    = list(genome = "mm10"),
#    Ensembl = list(release = 91),
#    BSg_pkg = "BSgenome.Mmusculus.UCSC.mm10",
#    TxDb = list(
#      source = "Ensembl GRCm38.p5 release 91 (December 2017)",
#      portal = "http://Dec2017.archive.ensembl.org/Mus_musculus",
#      gtf = paste0(
#        "ftp://ftp.ensembl.org/pub/",
#        "release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz"
#      )
#    ),
#    BioMart = list(
#      host     = "www.ensembl.org", # Dec2017.archive.ensembl.org
#      database = "ensembl",
#      dataset  = "mmusculus_gene_ensembl"
#    )
#  )
#  
#  # Gene features from Ensembl
#  MakeObj(EGA91_mouse, {
#    EGA91_mouse <- BuildGeneFeatures(organism)
#  })
#  
#  # CpG Islands from UCSC
#  MakeObj(CGI_mm10, {
#    CGI_mm10 <- ImportCpGIslandExt(organism$UCSC$genome)
#  })

## ----include=FALSE-------------------------------------------------------
if(CFG$UpdatePackageData) {
  devtools::use_data(EGA91_mouse, overwrite = T)
  devtools::use_data(CGI_mm10, overwrite = T)
}

## ----eval=FALSE----------------------------------------------------------
#  # =============================================================================.
#  # Mouse genome annotations (Ensembl NCBIM37 release 67)
#  # -----------------------------------------------------------------------------.
#  organism <- list(
#    name    = "Mus musculus",
#    taxid   = 10090,
#    UCSC    = list(genome = "mm9"),
#    Ensembl = list(release = 67),
#    BSg_pkg = "BSgenome.Mmusculus.UCSC.mm9",
#    TxDb = list(
#      source = "Ensembl NCBIM37 release 67 (May 2012)",
#      portal = "http://May2012.archive.ensembl.org/Mus_musculus",
#      gtf = paste0(
#        "ftp://ftp.ensembl.org/pub/",
#        "release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz"
#      )
#    ),
#    BioMart = list(
#      host     = "May2012.archive.ensembl.org",
#      database = "ensembl",
#      dataset  = "mmusculus_gene_ensembl"
#    )
#  )
#  
#  # Gene features from Ensembl
#  MakeObj(EGA67_mouse, {
#    EGA67_mouse <- BuildGeneFeatures(organism)
#  })
#  
#  # CpG Islands from UCSC
#  MakeObj(CGI_mm9, {
#    CGI_mm9 <- ImportCpGIslandExt(organism$UCSC$genome)
#  })

## ----include=FALSE-------------------------------------------------------
if(CFG$UpdatePackageData) {
  devtools::use_data(EGA67_mouse, overwrite = T)
  devtools::use_data(CGI_mm9, overwrite = T)
}

## ----eval=FALSE, include=FALSE, echo=FALSE-------------------------------
#  kable(EGA67_mouse$summary, caption = "Summary of gene features in EGA67_mouse")

## ----eval=FALSE----------------------------------------------------------
#  ANN <- EGA91_mouse
#  CGI <- CGI_mm10
#  BLK <- mm10_blacklist
#  
#  # Make 4kb windows every 1kb over the whole genome
#  WGT <- GenomicTiling(ANN$seqinfo, s = 1000, w = 4000)
#  
#  # Functional annotation of genomic bins
#  grp <- factor(
#    rep(1, length(WGT)), levels = 1:5, ordered = T,
#    labels = c("intergenic", "genic", "TSS", "CGI", "blacklist")
#  )
#  grp[which(countOverlaps(WGT, ANN$GNU) > 0)] <- "genic"
#  grp[which(countOverlaps(WGT, ANN$TSS) > 0)] <- "TSS"
#  grp[which(countOverlaps(WGT, CGI) > 0)] <- "CGI"
#  grp[which(countOverlaps(WGT, BLK) > 0)] <- "blacklist"
#  
#  WGT$category <- grp

## ----eval=FALSE----------------------------------------------------------
#  ANN <- EGA91_mouse
#  CGI <- CGI_mm10
#  BLK <- mm10_blacklist
#  
#  # Define promoter regions
#  TSS <- resize(ANN$TSS, width = 4000, fix = "center", use.names = F)
#  
#  # Define CGI regions as 4kb windows centered at CpG-islands annotated by UCSC
#  CGI <- resize(CGI, width = 4000, fix = "center", use.names = F)
#  CGI <- with(ANN, CleanupGRanges(CGI, seqinfo, organism = organism$name))
#  

## ----eval=FALSE, include=FALSE-------------------------------------------
#  kable(
#    data.frame(
#      `R object` = c("GNU", "TSS", "CGI", "WGT"),
#      description = c("Gene units", "Promoter regions", "CpG-islands", "Whole genome tiling"),
#      `number of regions` = sapply(list(ann, TSS, CGI, WGT), length),
#      `genomic coverage (bp)` = sapply(lapply(list(ANN$GNU, TSS, CGI, WGT), coverage), CoveredLength),
#      stringsAsFactors = F, check.names = F
#    ), caption = "Regions of interest for H3K27me3"
#  )

