# =============================================================================.
#' Whole genome tiling with regular bins (i.e. running windows)
# -----------------------------------------------------------------------------.
#' @param g
#' a \link{Seqinfo} object.
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

# HIDDEN #######################################################################

# =============================================================================.
#' Compute the genomic length covered by sequencing reads
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
#' Cleanup genomic ranges
# -----------------------------------------------------------------------------.
#' @description
#' Cleanup genomic ranges (e.g. for imported ChIP-seq peaks from bed files).
#'
#' @param grg
#' a \link{GRanges} object.
#'
#' @param seqinfo
#' a \link{Seqinfo} object.
#'
#' @param organism
#' character.
#'
#' @param blacklist
#' a \link{GRanges} object (default = NULL, none).
#'
#' @param keep.strand
#' logical (default = F, no).
#'
#' @return
#' \code{CleanupGRanges} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
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
#' Total coverage of genomic intervals
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
#' @keywords internal
#' @export
GenomicCoverage <- function(x, u = 1) {
  x <- width(FlatGRanges(x, seqinfo(x)))
  x <- sum(as.numeric(x))
  if(u != 1) x <- round(x / u)
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
