library(Tightrope)

# =============================================================================.
# EGA67_mouse
# -----------------------------------------------------------------------------.
organism <- list(
  name    = "Mus musculus",
  taxid   = 10090,
  UCSC    = list(
    genome = "mm9",
    cpgIslandExt = "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cpgIslandExt.txt.gz"
  ),
  Ensembl = list(release = 91),
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
devtools::use_data(EGA67_mouse)
devtools::use_data(CGI_mm9)

# =============================================================================.
# EGA91_mouse
# -----------------------------------------------------------------------------.
organism <- list(
  name    = "Mus musculus",
  taxid   = 10090,
  UCSC    = list(
    genome = "mm10",
    cpgIslandExt = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz"
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
devtools::use_data(EGA91_mouse)
devtools::use_data(CGI_mm10)
