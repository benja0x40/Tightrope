# =============================================================================.
#' Whole genome tiling with regular bins (i.e. running windows)
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportGenomicRanges},
#'   \link{ImportCpGIslandExt},
#'   \link{BuildGeneFeatures},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param g
#' \link{Seqinfo} object for the considered organism.
#'
#' @param s
#' integer, spacing between genomic bins in bp.
#'
#' @param w
#' integer, width of genomic bins in bp.
#'
#' @return
#' \code{GenomicTiling} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @export
GenomicTiling <- function(g, s, w = NULL) {
  if(is.null(w)) w <- s
  grg <- tileGenome(g, tilewidth = s, cut.last.tile.in.chrom = T)
  if(s != w) {
    grg <- suppressWarnings(resize(grg, width = w, use.names = F))
    grg <- trim(grg, use.names = F)
    grg <- grg[width(grg) == w]
  }
  grg
}

# =============================================================================.
#' Import tab delimited text files representing genomic interval data
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportCpGIslandExt},
#'   \link{BuildGeneFeatures},
#'   \link{GenomicTiling},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @inheritParams df2grg
#'
#' @param fpath
#' path to a text file defining genomic ranges.
#'
#' @param seqinfo
#' \link{Seqinfo} object for the considered organism.
#'
#' @param header
#' logical indicating if column names are present (default = F).
#'
#' @param ...
#' optional parameters forwarded to the \link{read.delim} function.
#'
#' @return
#' \code{ImportGenomicRanges} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @export
ImportGenomicRanges <- function(
  fpath, chr = 1, start = 2, end = 3, strand = NULL, seqinfo = NULL,
  xidx = NULL, xlbl = NULL, header = F, ...
) {
  con <- file(fpath)
  grg <- read.delim(con, header = header, stringsAsFactors=F, ...)
  grg <- df2grg(grg, chr, start, end, strand, xidx, xlbl)
  if(! is.null(seqinfo)) seqinfo(grg) <- seqinfo[seqlevels(grg)]
  grg
}

# =============================================================================.
#' Cleanup genomic ranges
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportGenomicRanges},
#'   \link{ImportCpGIslandExt}
# -----------------------------------------------------------------------------.
#' @description
#' Cleanup genomic ranges such as imported ChIP-seq peaks from bed files by
#' removing ranges corresponding to non canonical chromosomes.
#'
#' @param grg
#' a \link{GRanges} object.
#'
#' @param seqinfo
#' \link{Seqinfo} object for the considered organism.
#'
#' @param organism
#' name of the organism (e.g. "Homo sapiens", "Mus musculus").
#'
#' @param blacklist
#' an optional \link{GRanges} object used to filter out any genomic ranges
#' overlaping with the blacklisted one (default = NULL, none).
#'
#' @param keep.strand
#' logical, when TRUE strand information is retained, when FALSE (default)
#' strand information is discarded (i.e. strand = "*").
#'
#' @return
#' \code{CleanupGRanges} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @export
CleanupGRanges <- function(
  grg, seqinfo, organism, blacklist = NULL, keep.strand = F
) {

  if(! keep.strand) strand(grg)  <- "*"
  seqinfo(grg) <- seqinfo[seqlevels(grg)]
  species <- gsub(" ", "_", organism)

  grg <- keepStandardChromosomes(
    grg, species = species, pruning.mode = "coarse"
  )

  if(! is.null(blacklist)) {
    grg <- grg[countOverlaps(grg, blacklist) == 0]
  }

  grg
}

# =============================================================================.
#' Total coverage of genomic intervals
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{CoveredLength}
# -----------------------------------------------------------------------------.
#' @param x
#' a \link{GRanges} object.
#'
#' @param u
#' numeric, base unit value (e.g 1 => bp, 1E3 => kb, 1E6 => Mb).
#' By default u = 1 (bp).
#'
#' \code{GenomicCoverage} returns the total coverage of x.
# -----------------------------------------------------------------------------.
#' @export
GenomicCoverage <- function(x, u = 1) {
  x <- width(FlatGRanges(x, seqinfo(x)))
  x <- sum(as.numeric(x))
  if(u != 1) x <- round(x / u)
  x
}

# HIDDEN #######################################################################

# =============================================================================.
#' Compute the genomic length covered by sequence reads
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{GenomicCoverage}
# -----------------------------------------------------------------------------.
#' @param cvg
#' a coverage object, see \link{GenomicAlignments::coverage}.
#'
#' @return
#' \code{CoveredLength} returns the genomic length covered by sequencing reads.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
CoveredLength <- function(cvg) {
  s <- 0L
  for(chr in names(cvg)) {
    i <- which(runValue(cvg[[chr]]) > 0)
    s <- s + sum(as.numeric(runLength(cvg[[chr]])[i]))
  }
  s
}

# =============================================================================.
#' Reduce genomic intervals
# -----------------------------------------------------------------------------.
#' @param x
#' a \link{GRanges} object.
#'
#' @param seqinfo
#' a \link{Seqinfo} object.
#'
#' @return
#' \code{FlatGRanges} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
FlatGRanges <- function(x, seqinfo) {
  x <- unlist(x)
  strand(x) <- "*"
  x <- reduce(x)
  seqinfo(x) <- seqinfo[seqlevels(x)]
  x
}

# =============================================================================.
#' Consensus genomic ranges
# -----------------------------------------------------------------------------.
#' @inheritParams GenomicRanges::countOverlaps
#'
#' @param grl
#' a list of \link{GRanges} objects refering to the same genome build.
#'
#' @param seqinfo
#' a \link{Seqinfo} object.
#'
#' @return
#' \code{ConsensusGRanges} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ConsensusGRanges <- function(grl, minoverlap, ignore.strand = T, seqinfo = NULL) {
  # Merge peak lists
  r <- GRanges(seqinfo = seqinfo)
  for(lbl in names(grl)) {
    r <- GenomicRanges::union(r, grl[[lbl]], ignore.strand = ignore.strand)
  }
  # Split peaks such that none can overlap more than one peak from another set
  r <- disjoin(r)
  # Track original coordintates of disjoin peaks
  r$initial_coordinates <- paste(
    seqnames(r), start(r), end(r), width(r), sep = "_"
  )
  # ---------------------------------------------------------------------------.
  for(lbl in names(grl)) {
    # Tag overlapping peaks
    ovl <- GenomicRanges::countOverlaps(
      r, grl[[lbl]], ignore.strand = ignore.strand, minoverlap = minoverlap
    )
    mcols(r)[[lbl]] <- ovl > 0
    # Refine coordinates of overlapping peaks (safe intersect thanks to disjoin)
    idx <- GenomicRanges::findOverlaps(
      r, grl[[lbl]], ignore.strand = ignore.strand, minoverlap = minoverlap
    )
    qhx <- queryHits(idx)
    shx <- subjectHits(idx)
    cmn <- GenomicRanges::pintersect(
      r[qhx], grl[[lbl]][shx], ignore.strand = ignore.strand,
      drop.nohit.ranges = F
    )
    cmn$hit <- NULL
    r[qhx] <- cmn
  }
  r
}

# NOT EXPORTED #################################################################

# =============================================================================.
#' Convert \code{data.frame} to \link{GRanges}
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportGenomicRanges}
# -----------------------------------------------------------------------------.
#' @param df
#' a data.frame defining genomic ranges.
#'
#' @param chr
#' column for chromsome identifiers (default = 1).
#'
#' @param start
#' column for start position (default = 2).
#'
#' @param end
#' column for end position (default = 3).
#'
#' @param strand
#' column for strand (default = NULL, e.g. none).
#'
#' @param xidx
#' list of extra column indexes to be imported.
#' By default all available columns are imported.
#'
#' @param xlbl
#' character vector defining the name of imported extra columns
#' (default = NULL).
#'
#' @return
#' \code{df2grg} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
df2grg <- function(
  df, chr = 1, start = 2, end = 3, strand = NULL, xidx = NULL, xlbl = NULL
) {

  if(is.character(chr))    chr    <- match(chr, colnames(df))
  if(is.character(start))  start  <- match(start, colnames(df))
  if(is.character(end))    end    <- match(end, colnames(df))
  if(is.character(strand)) strand <- match(strand, colnames(df))

  # Copy all extra data by default
  if(is.null(xidx)) {
    xidx <- which(! (1:ncol(df)) %in% c(chr, start, end, strand))
  }

  df <- df[,c(chr, start, end, strand, xidx)]
  grg <- GRanges(
    seqnames = df[,1],
    ranges = IRanges(start=df[,2], end=df[,3], names=NULL)
  )

  if(! is.null(strand)) {
    strand(grg) <- df[,4]
  }

  mcols(grg) <- df[,(5-is.null(strand)):ncol(df)]
  if(! is.null(xlbl)) {
    colnames(mcols(grg)) <- xlbl
  }

  grg
}
