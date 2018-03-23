# =============================================================================.
#' Build an R environment containing genomic features from Ensembl
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportCpGIslandExt},
#'   \link{ImportGenomicRanges},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param organism
#' list.
#'
#' @param gene_id
#' character (default = "ensembl_gene_id").
#'
#' @param gene_name
#' character (default = "hgnc_symbol").
#'
#' @return
#' \code{BuildGeneFeatures} returns an \link{environment}.
# -----------------------------------------------------------------------------.
#' @export
BuildGeneFeatures <- function(
  organism, gene_id = "ensembl_gene_id", gene_name = "hgnc_symbol"
) {

  # ---------------------------------------------------------------------------.
  # Prepare environment
  EGA <- new.env(parent = globalenv(), size = NA)
  EGA$timestamp <- format(Sys.time(), '%A %d.%m.%Y %H:%M:%S')
  EGA$organism  <- organism
  EGA$summary    <- data.frame(
    matrix(c(
      "INTERGENIC", "intergenic regions",
      "GNU",        "gene units",
      "TXU",        "transcribed units",
      "TSS",        "transcription start sites",
      "TTS",        "transcription termination sites",
      "CCS",        "constitutive coding sequences (i.e. 'purely' exon)",
      "CNC",        "constitutive intragenic non-coding sequences",
      "ACS",        "alternative coding sequence (not CCS or CNC)"
    ), nrow = 8, byrow = T, dimnames = list(NULL, c("object", "description"))),
    stringsAsFactors = F
  )
  # ---------------------------------------------------------------------------.
  message("Retrieving chromosome information from UCSC...")
  EGA$organism$UCSC$chrominfo <- fetchExtendedChromInfoFromUCSC(
    genome = EGA$organism$UCSC$genome
  )
  EGA$seqinfo <- with(
    EGA$organism$UCSC$chrominfo, Seqinfo(
      seqnames = UCSC_seqlevel, seqlengths = UCSC_seqlength,
      isCircular = circular, genome = EGA$organism$UCSC$genome
    )
  )
  # ---------------------------------------------------------------------------.
  # Make TXDB directly from Ensembl annotations (gtf file download)
  message("Downloading genes from Ensembl and making TxDb...")
  TXDB <- with(
    EGA$organism, makeTxDbFromGFF(
      TxDb$gtf, organism = name, dataSource = TxDb$source
    )
  )
  seqlevelsStyle(TXDB) <- "UCSC"
  TXDB <- keepStandardChromosomes(TXDB, species = gsub(" ", "_", organism$name))
  # ---------------------------------------------------------------------------.
  message("Building gene units...")
  EGA$GNU <- genes(TXDB, columns = c("gene_id", "cds_id"))
  EGA$GNU$gene_name    <- ""
  EGA$GNU$gene_biotype <- ""
  # ---------------------------------------------------------------------------.
  message("Retrieving gene annotations from Ensembl/BioMart...")
  mrt <- with(
    EGA$organism$BioMart, useMart(
      host = host, biomart = database, dataset = dataset
    )
  )
  lbl <- c(gene_id, gene_name, "gene_biotype")
  res <- getBM(
    attributes = lbl, filters = gene_id, values = EGA$GNU$gene_id, mart = mrt
  )
  # Update gene annotation
  idx <- match(res[[gene_id]], EGA$GNU$gene_id)
  EGA$GNU$gene_name[idx]    <- res[[gene_name]]
  EGA$GNU$gene_biotype[idx] <- res$gene_biotype
  seqinfo(EGA$GNU) <- with(EGA, seqinfo[seqlevels(GNU)])
  # ---------------------------------------------------------------------------.
  message("Building transcribed units...")
  EGA$TXU <- transcripts(TXDB, columns = c("gene_id", "tx_name"))
  EGA$TXU$gene_id <- as.character(EGA$TXU$gene_id)
  EGA$TXU$transcript_id <- EGA$TXU$tx_name
  idx <- match(EGA$TXU$gene_id, EGA$GNU$gene_id)
  EGA$TXU$gene_biotype <- EGA$GNU$gene_biotype[idx]
  EGA$TXU$gene_name    <- EGA$GNU$gene_name[idx]
  EGA$TXU$tx_name <- NULL
  seqinfo(EGA$TXU) <- with(EGA, seqinfo[seqlevels(TXU)])
  # ---------------------------------------------------------------------------.
  message("Extracting Transcription Start/Termination Sites (TSS and TTS)...")
  EGA$TSS <- promoters(EGA$TXU, upstream = 0, downstream = 0)
  EGA$TTS <- promoters(invertStrand(EGA$TXU), upstream = 0, downstream = 0)
  # ---------------------------------------------------------------------------.
  message("Extracting coding and non-coding genic regions...")
  rex <- FlatGRanges(cds(TXDB), EGA$seqinfo)
  rin <- FlatGRanges(intronicParts(TXDB), EGA$seqinfo)
  r5u <- FlatGRanges(fiveUTRsByTranscript(TXDB), EGA$seqinfo)
  r3u <- FlatGRanges(threeUTRsByTranscript(TXDB), EGA$seqinfo)
  EGA$CCS <- GenomicRanges::setdiff(
    rex, reduce(GenomicRanges::union(rin, GenomicRanges::union(r5u, r3u)))
  )
  EGA$CNC <- GenomicRanges::setdiff(FlatGRanges(EGA$GNU, EGA$seqinfo), rex)
  EGA$ACS <- GenomicRanges::setdiff(
    FlatGRanges(EGA$GNU, EGA$seqinfo), with(EGA, reduce(union(CCS, CNC)))
  )
  # ---------------------------------------------------------------------------.
  message("Building intergenic regions...")
  EGA$INTERGENIC  <- EGA$GNU
  strand(EGA$INTERGENIC) <- "*"
  EGA$INTERGENIC <- gaps(reduce(EGA$INTERGENIC))
  EGA$INTERGENIC <- EGA$INTERGENIC[strand(EGA$INTERGENIC) == "*"]
  # ---------------------------------------------------------------------------.
  EGA$summary$size <- sapply(
    EGA$summary$object, function(x) length(EGA[[x]])
  )
  EGA$summary$coverage <- sapply(
    EGA$summary$object, function(x) GenomicCoverage(EGA[[x]])
  )
  # ---------------------------------------------------------------------------.
  EGA
}
