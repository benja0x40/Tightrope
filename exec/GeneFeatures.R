library(Tightrope)

# HUMAN ########################################################################

# =============================================================================.
# EGA91_human
# -----------------------------------------------------------------------------.
organism <- list(
  name    = "Homo sapiens",
  taxid   = 9606,
  UCSC    = list(
    genome = "hg38",
    cpgIslandExt = paste0(
      "http://hgdownload.soe.ucsc.edu/goldenPath/", "hg38",
      "/database/cpgIslandExt.txt.gz"
    )
  ),
  Ensembl = list(release = 91),
  BSg_pkg = "BSgenome.Hsapiens.UCSC.hg38",
  TxDb = list(
    source = "Ensembl GRCh38.p10 release 91 (December 2017)",
    portal = "http://Dec2017.archive.ensembl.org/Homo_sapiens",
    gtf = paste0(
      "ftp://ftp.ensembl.org/pub/",
      "release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz"
    )
  ),
  BioMart = list(
    host     = "www.ensembl.org", # Dec2017.archive.ensembl.org
    database = "ensembl",
    dataset  = "hsapiens_gene_ensembl"
  )
)
# -----------------------------------------------------------------------------.
EGA91_human <- BuildGeneFeatures(organism)
SaveObj(EGA91_human)
# -----------------------------------------------------------------------------.
download.file(organism$UCSC$cpgIslandExt, destfile = "CGI_hg38.txt.gz")
CGI_hg38 <- ImportCpGIslandExt("CGI_hg38.txt.gz", seqinfo = EGA91_human$seqinfo)
SaveObj(CGI_hg38)
# -----------------------------------------------------------------------------.
devtools::use_data(EGA91_human, overwrite = T)
devtools::use_data(CGI_hg38, overwrite = T)

# =============================================================================.
# EGA75_human
# -----------------------------------------------------------------------------.
organism <- list(
  name  = "Homo sapiens",
  taxid = 9606,
  UCSC = list(
    genome = "hg19",
    cpgIslandExt = paste0(
      "http://hgdownload.soe.ucsc.edu/goldenPath/", "hg19",
      "/database/cpgIslandExt.txt.gz"
    )
  ),
  Ensembl = list(release = 75),
  BSg_pkg = "BSgenome.Hsapiens.UCSC.hg19",
  TxDb = list(
    source = "Ensembl GRCh37.p13 release 75 (February 2014)",
    portal = "http://Feb2014.archive.ensembl.org/Homo_sapiens",
    gtf = paste0(
      "ftp://ftp.ensembl.org/pub/",
      "release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
    )
  ),
  BioMart = list(
    host     = "Feb2014.archive.ensembl.org", # "grch37.ensembl.org"
    database = "ensembl",
    dataset  = "hsapiens_gene_ensembl"
  )
)
# -----------------------------------------------------------------------------.
EGA75_human <- BuildGeneFeatures(organism)
SaveObj(EGA75_human)
# -----------------------------------------------------------------------------.
download.file(organism$UCSC$cpgIslandExt, destfile = "CGI_hg19.txt.gz")
CGI_hg19 <- ImportCpGIslandExt("CGI_hg19.txt.gz", seqinfo = EGA75_human$seqinfo)
SaveObj(CGI_hg19)
# -----------------------------------------------------------------------------.
devtools::use_data(EGA75_human, overwrite = T)
devtools::use_data(CGI_hg19, overwrite = T)

# MOUSE ########################################################################

# =============================================================================.
# EGA91_mouse
# -----------------------------------------------------------------------------.
organism <- list(
  name    = "Mus musculus",
  taxid   = 10090,
  UCSC    = list(
    genome = "mm10",
    cpgIslandExt = paste0(
      "http://hgdownload.soe.ucsc.edu/goldenPath/", "mm10",
      "/database/cpgIslandExt.txt.gz"
    )
  ),
  Ensembl = list(release = 91),
  BSg_pkg = "BSgenome.Mmusculus.UCSC.mm10",
  TxDb = list(
    source = "Ensembl GRCm38.p5 release 91 (December 2017)",
    portal = "http://Dec2017.archive.ensembl.org/Mus_musculus",
    gtf = paste0(
      "ftp://ftp.ensembl.org/pub/",
      "release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz"
    )
  ),
  BioMart = list(
    host     = "www.ensembl.org", # Dec2017.archive.ensembl.org
    database = "ensembl",
    dataset  = "mmusculus_gene_ensembl"
  )
)
# -----------------------------------------------------------------------------.
EGA91_mouse <- BuildGeneFeatures(organism)
SaveObj(EGA91_mouse)
# -----------------------------------------------------------------------------.
download.file(organism$UCSC$cpgIslandExt, destfile = "CGI_mm10.txt.gz")
CGI_mm10 <- ImportCpGIslandExt("CGI_mm10.txt.gz", seqinfo = EGA91_mouse$seqinfo)
SaveObj(CGI_mm10)
# -----------------------------------------------------------------------------.
devtools::use_data(EGA91_mouse, overwrite = T)
devtools::use_data(CGI_mm10, overwrite = T)

# =============================================================================.
# EGA67_mouse
# -----------------------------------------------------------------------------.
organism <- list(
  name    = "Mus musculus",
  taxid   = 10090,
  UCSC    = list(
    genome = "mm9",
    cpgIslandExt = paste0(
      "http://hgdownload.soe.ucsc.edu/goldenPath/", "mm9",
      "/database/cpgIslandExt.txt.gz"
    )
  ),
  Ensembl = list(release = 67),
  BSg_pkg = "BSgenome.Mmusculus.UCSC.mm9",
  TxDb = list(
    source = "Ensembl NCBIM37 release 67 (May 2012)",
    portal = "http://May2012.archive.ensembl.org/Mus_musculus",
    gtf = paste0(
      "ftp://ftp.ensembl.org/pub/",
      "release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz"
    )
  ),
  BioMart = list(
    host     = "May2012.archive.ensembl.org",
    database = "ensembl",
    dataset  = "mmusculus_gene_ensembl"
  )
)
# -----------------------------------------------------------------------------.
EGA67_mouse <- BuildGeneFeatures(organism, gene_name = "mgi_symbol")
SaveObj(EGA67_mouse)
# -----------------------------------------------------------------------------.
download.file(organism$UCSC$cpgIslandExt, destfile = "CGI_mm9.txt.gz")
CGI_mm9 <- ImportCpGIslandExt("CGI_mm9.txt.gz", seqinfo = EGA67_mouse$seqinfo)
SaveObj(CGI_mm9)
# -----------------------------------------------------------------------------.
devtools::use_data(EGA67_mouse, overwrite = T)
devtools::use_data(CGI_mm9, overwrite = T)

