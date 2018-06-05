## ----hidden_dependencies, include=FALSE----------------------------------
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(shades))

suppressPackageStartupMessages(library(LittleThumb))
suppressPackageStartupMessages(library(Barbouille))
suppressPackageStartupMessages(library(Tightrope))

## ----hidden_functions, include=FALSE-------------------------------------
# =============================================================================.
# Markdown links
# -----------------------------------------------------------------------------.
MarkdownLink <- function(x, txt = NULL, lnk = "") {
  if(is.null(txt)) txt <- x
  lnk <- paste0("[", txt, "](", lnk, x, ")")
  chk <- x == ""
  lnk[chk] <- txt[chk]
  lnk
}

## ----hidden_configuration, include=FALSE---------------------------------
# =============================================================================.
#
# -----------------------------------------------------------------------------.
CFG <- new.env()

CFG$path <- "~/OneDrive/_Work_/@BRIC - Copenhagen/ChIP-seq spike-in & Tightrope/_R_ANALYSES_BRD/analyses/Tightrope/cache/"

CFG$images <- list(
  
  build = F,
  
  path    = "images/BRD/",
  name    = "img",
  counter = 0,
  
  resolution = 300,
  width      = 4,
  height     = 4.25,
  units      = "in",
  
  par = list(pch = 20, mar = c(4.5, 4.5, 4.0, 2.0))
)

## ----overview, echo=FALSE, out.width="80%", fig.align='center'-----------
include_graphics("images/Overview.png")

## ------------------------------------------------------------------------
 # Load Tightrope package (this will produce numerous messages)
library(Tightrope)

## ----eval=FALSE----------------------------------------------------------
#  help(package = "Tightrope")

## ----eval=FALSE----------------------------------------------------------
#  vignette("BRD", package = "Tightrope")

## ----include=FALSE-------------------------------------------------------
library(LittleThumb)
LittleThumb(rootpath = "AutoSaved")

## ----eval=FALSE----------------------------------------------------------
#  vignette("Annotations", package = "Tightrope")

## ------------------------------------------------------------------------
# Load mouse genome annotations (mm10)
data("EGA91_mouse")    # Gene features from Ensembl      (environment)
data("CGI_mm10")       # CpG Islands from UCSC           (GRanges)
data("mm10_blacklist") # Blacklisted regions from ENCODE (GRanges)

## ----message=FALSE-------------------------------------------------------
# Load whole genome tiling object WGT or create it and save it if necessary
MakeObj(WGT, {
  # Make 4kb windows every 1kb over the whole genome (GRanges object)
  WGT <- GenomicTiling(EGA91_mouse$seqinfo, s = 1000, w = 4000)
  
  WGT <- CleanupGRanges(
    WGT, EGA91_mouse$seqinfo, "Mus musculus", mm10_blacklist
  )
  
  # Functional annotation of genomic bins
  grp <- factor(
    rep(1, length(WGT)), levels = 1:5, ordered = T,
    labels = c("intergenic", "genic", "TSS", "CGI", "blacklist")
  )
  grp[which(countOverlaps(WGT, EGA91_mouse$GNU) > 0)] <- "genic"
  grp[which(countOverlaps(WGT, EGA91_mouse$TSS) > 0)] <- "TSS"
  grp[which(countOverlaps(WGT, CGI_mm10) > 0)] <- "CGI"
  grp[which(countOverlaps(WGT, mm10_blacklist) > 0)] <- "blacklist"
  
  WGT$category <- grp
})

## ----eval=FALSE----------------------------------------------------------
#  # Define bam file paths and the corresponding names of experimental conditions
#  bam.files <- c(
#    H3WT_DMSO_H3K27me3_R1 = "~/sequencing/sample1.bam",
#    K27M_DMSO_H3K27me3_R1 = "~/sequencing/sample2.bam",
#    H3WT_DMSO_Input_R1    = "~/sequencing/sample3.bam",
#    K27M_DMSO_Input_R1    = "~/sequencing/sample4.bam"
#    # etc.
#  )
#  conditions <- names(bam.files)

## ----eval=FALSE, message=FALSE-------------------------------------------
#  # Define accepted reads (single-end sequencing)
#  bam.flag <- scanBamFlag(isDuplicate = F, isUnmappedQuery = F)
#  
#  # Read counts over genomic tiling bins WGT (single-end sequencing)
#  COUNTS <- ReadCountMatrix(
#      bam.files, WGT, paired = F, names = conditions,
#      param = ScanBamParam(flag = bam.flag)
#  )

## ----eval=FALSE----------------------------------------------------------
#  # Define accepted reads (paired-end sequencing)
#  bam.flag <- scanBamFlag(
#    isDuplicate = F, isUnmappedQuery = F, isProperPair = T
#  )
#  
#  # Read counts over genomic tiling bins WGT (paired-end sequencing)
#  COUNTS <- ReadCountMatrix(
#      bam.files, WGT, paired = T, names = conditions,
#      param = ScanBamParam(flag = bam.flag)
#  )

## ----eval=FALSE, include=FALSE-------------------------------------------
#  # Define ChIP and Input conditions
#  chip <- c("H3K27me3_WT", "H3K27me3_K27M")
#  ctrl <- c("Input_WT", "Input_K27M")
#  
#  # Precomputed number of mapped reads from chromatin spike-in
#  spikein.reads <- c(
#    Input_WT      = 308729,
#    H3K27me3_WT   = 1640186,
#    Input_K27M    = 209380,
#    H3K27me3_K27M = 9088680
#  )
#  
#  # Scaling factors based on chromatin spike-in
#  spikein.factors <- mean(spikein.reads) / spikein.reads

## ----reshape_dataset, eval=FALSE, message=FALSE, include=FALSE-----------
#  MakeObj(NSC_K27M, {
#    # Rename genomic annotations
#    ANN <- EGA91_mouse
#    CGI <- CGI_mm10
#    BLK <- mm10_blacklist
#  
#    # Cleanup (optional)
#    # rm(EGA91_mouse, CGI_mm10, mm10_blacklist)
#  
#    data("NSC_K27M_METADATA") # Annotation of experimental conditions
#    data("NSC_K27M_COUNTS")   # Precomputed read count matrixes
#    data("NSC_K27M_SPIKEIN")  # Number of Drosophila reads (spike-in)
#  
#    # Load precomputed whole genome tiling object WGT
#    WGT <- readRDS("~/Google Drive/Work/@BRIC - Copenhagen/ChIP-seq spike-in & Tightrope/_R_ANALYSES_BRD/analyses/NextSeq_K27M_spike_in/data/WGT.rdata")
#  
#    # Functional annotation of genomic bins
#    grp <- factor(
#      rep(1, length(WGT)), levels = 1:5, ordered = T,
#      labels = c("intergenic", "genic", "TSS", "CGI", "blacklist")
#    )
#    grp[which(countOverlaps(WGT, ANN$GNU) > 0)] <- "genic"
#    grp[which(countOverlaps(WGT, ANN$TSS) > 0)] <- "TSS"
#    grp[which(countOverlaps(WGT, CGI) > 0)] <- "CGI"
#    grp[which(countOverlaps(WGT, BLK) > 0)] <- "blacklist"
#  
#    WGT$category <- grp
#  
#    NSC_K27M <- new.env()
#    NSC_K27M$METADATA <- NSC_K27M_METADATA
#    NSC_K27M$SPIKEIN  <- NSC_K27M_SPIKEIN
#    NSC_K27M$COUNTS   <- NSC_K27M_COUNTS$WGT
#    NSC_K27M$WGT      <- WGT
#  
#    # Cleanup
#    rm(grp, NSC_K27M_METADATA, NSC_K27M_COUNTS, NSC_K27M_SPIKEIN)
#    gc()
#  })

## ------------------------------------------------------------------------
# Load the NSC dataset (replicate of experiments from Mohammad et al. 2017)
data("NSC_K27M")
library(data.table)

## ----eval=FALSE----------------------------------------------------------
#  NSC_K27M$METADATA # Annotation of experimental conditions (data.table object)
#  NSC_K27M$SPIKEIN  # Number of Drosophila reads (chromatin spike-in)
#  NSC_K27M$COUNTS   # Precomputed read count matrixes over genomic tiling bins
#  NSC_K27M$WGT      # Whole genome tiling bins (GRanges object)

## ----echo=FALSE----------------------------------------------------------
kable(NSC_K27M$METADATA[, -(1:3)])

## ----message=FALSE-------------------------------------------------------
COUNTS <- NSC_K27M$COUNTS
WGT    <- NSC_K27M$WGT

## ----echo=FALSE----------------------------------------------------------
tbl <- colSums(COUNTS)
kable(
  data.frame(
    sample_id = names(tbl),
    number    = as.vector(tbl)
  ),
  caption = "Total read counts"
)
suppressWarnings(rm(tbl))

## ------------------------------------------------------------------------
chip <- NSC_K27M$METADATA[antibody == "H3K27me3"]$sample_id
ctrl <- NSC_K27M$METADATA[antibody == "Input"]$sample_id

## ------------------------------------------------------------------------
# log2 transformed counts with dithering
MakeObj(l2c, {
  l2c <- log2(DitherCounts(COUNTS))
  l2c <- RemoveInputBias(l2c, controls = ctrl, as.log2 = T)
})

# RPM
MakeObj(rpm, {
  rpm <- log2(colSums(COUNTS))
  rpm <- t(t(l2c) - rpm + mean(rpm))
})

## ----eval=FALSE----------------------------------------------------------
#  PlotCountDistributions(
#    rpm[, chip], pops = WGT$category, sampling = 5E5, main = "Genomic bins"
#  )
#  chk <- FiniteValues(rpm[, ctrl])
#  abline(h = median(rpm[chk, ctrl]), lwd = 1.5, lty = 3)

## ----echo=FALSE, message=FALSE, out.width="60%", fig.align="center"------
img <- list("GenomicBins", h = 6, w = 8)
MkImg({
  PlotCountDistributions(
    rpm[, chip], pops = WGT$category, sampling = 5E5, main = "Genomic bins"
  )
  chk <- FiniteValues(rpm[, ctrl])
  abline(h = median(rpm[chk, ctrl]), lwd = 1.5, lty = 3)
})

## ----echo=FALSE, message=FALSE, out.width="25%"--------------------------
img <- list("ScatterPlots")

grd <- "hcl.duo.light"
grp <- "intergenic"
clr <- c("grey", shades::gradient("magma", steps = nlevels(WGT$category)))
names(clr) <- c(levels(WGT$category), NA)

rng <- range(rpm[FiniteValues(rpm), ])
x <- chip[c(1, 5)]
y <- chip[c(2, 6)]
xlab <- "WT"
ylab <- "K27M"

idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    rpm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "genic"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    rpm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "TSS"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    rpm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "CGI"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    rpm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})

## ----eval=FALSE----------------------------------------------------------
#  # Prepare figure layout
#  layout(matrix(1:6, 2, 3, byrow = F))
#  
#  # Run BRD with default parameters and plot the density graph
#  PlotBRD(
#    BRD(cnt = COUNTS[WGT$category == "genic", ], controls = ctrl, sampling = 1E5),
#    plots = c("density", "thresholds"), title = "using genes"
#  )
#  PlotBRD(
#    BRD(cnt = COUNTS[WGT$category == "TSS", ], controls = ctrl, sampling = 1E5),
#    plots = c("density", "thresholds"), title = "using promoters"
#  )
#  PlotBRD(
#    BRD(cnt = COUNTS[WGT$category == "CGI", ], controls = ctrl, sampling = 1E5),
#    plots = c("density", "thresholds"), title = "using CpG-islands"
#  )

## ----echo=FALSE, message=FALSE, fig.align="center", out.width="75%"------
img <- list("PreviewBRD", h = 6.375, w = 9)
MkImg({
  layout(matrix(1:6, 2, 3, byrow = F))
  PlotBRD(
    BRD(cnt = COUNTS[WGT$category == "genic", ], controls = ctrl, sampling = 1E5),
    plots = c("density", "thresholds"), title = "using genes"
  )
  PlotBRD(
    BRD(cnt = COUNTS[WGT$category == "TSS", ], controls = ctrl, sampling = 1E5),
    plots = c("density", "thresholds"), title = "using promoters"
  )
  PlotBRD(
    BRD(cnt = COUNTS[WGT$category == "CGI", ], controls = ctrl, sampling = 1E5),
    plots = c("density", "thresholds"), title = "using CpG-islands"
  )
})

## ----brd_estimation, eval = FALSE----------------------------------------
#  # Run BRD using read counts over gene units (unimodal density)
#  brd.gnu <- BRD(
#    cnt = COUNTS[WGT$category == "genic", ], controls = ctrl, ncl = 1,
#    sampling = 3E5, dither = 10
#  )
#  # Run BRD using read counts over promoters (bimodal density)
#  brd.tss <- BRD(
#    cnt = COUNTS[WGT$category == "TSS", ], controls = ctrl, ncl = 1,
#    sampling = 3E5, dither = 10
#  )
#  # Run BRD using read counts over CpG-islands (unimodal density)
#  brd.cgi <- BRD(
#    cnt = COUNTS[WGT$category == "CGI", ], controls = ctrl, ncl = 1,
#    sampling = 3E5, dither = 10
#  )

## ----include=FALSE-------------------------------------------------------
# Run BRD using read counts over gene units (unimodal density)
MakeObj(brd.gnu, {
  brd.gnu <- BRD(
    cnt = COUNTS[WGT$category == "genic", ], controls = ctrl, ncl = 1,
    sampling = 3E5, dither = 10
  )
})
# Run BRD using read counts over promoters (bimodal density)
MakeObj(brd.tss, {
  brd.tss <- BRD(
    cnt = COUNTS[WGT$category == "TSS", ], controls = ctrl, ncl = 1,
    sampling = 3E5, dither = 10
  )
})
# Run BRD using read counts over CpG-islands (unimodal density)
MakeObj(brd.cgi, {
  brd.cgi <- BRD(
    cnt = COUNTS[WGT$category == "CGI", ], controls = ctrl, ncl = 1,
    sampling = 3E5, dither = 10
  )
})

## ----eval=FALSE----------------------------------------------------------
#  # Value of BRD scaling factors
#  ScalingFactors(brd.gnu) # Estimated using gene units
#  ScalingFactors(brd.tss) # Estimated using promoters
#  ScalingFactors(brd.cgi) # Estimated using CpG-islands

## ----echo=FALSE----------------------------------------------------------
tbl <- data.frame(
  genes = ScalingFactors(brd.gnu),
  promoters = ScalingFactors(brd.tss),
  CGIs = ScalingFactors(brd.cgi)
)

tbl$`delta max. (%)` <- apply(
  tbl[, 1:3], 1, function(x) round(100 * diff(range(x)) / mean(x), 1)
)
tbl$`TSS vs CGI (%)` <- apply(
  tbl[, 2:3], 1, function(x) round(100 * diff(range(x)) / mean(x), 1)
)

kable(tbl, caption = "BRD scaling factors")

## ----echo=FALSE, fig.height=6, out.width="45%"---------------------------

# Scaling factors based on chromatin spike-in
spikein.factors <- mean(NSC_K27M$SPIKEIN) / NSC_K27M$SPIKEIN
spikein.factors <- spikein.factors[colnames(COUNTS)]

x <- factor(rownames(tbl), levels = rownames(tbl))
chk <- x %in% chip

# spikein.factors <- spikein.factors / median(spikein.factors[chk]) * median(as.matrix(tbl[chk, 1:3]))

par(mar = par()$mar + c(8, 0, 0, 6))
clr <- c(grey(0), rgb(0.9, 0, 0), rgb(1, 0.5, 0), rgb(0, 0.5, 1))
names(clr) <- c(colnames(tbl)[1:3], "spike-in")

plot(x, rep(0, length(x)), ylim = c(0, 1.1 * max(tbl[, 1:3], spikein.factors)), las = 2, type = "n", col = NA, border = NA, ylab = "scaling factor", main = "unmatched")

for(idx in list(1:4, 5:8, 9:12, 13:16)) {
  points(x[idx], tbl$genes[idx], type = "b", col = clr[1], lwd = 2, pch = 20)
  points(x[idx], tbl$promoters[idx], type = "b", col = clr[2], lwd = 2, pch = 20)
  points(x[idx], tbl$CGIs[idx], type = "b", col =  clr[3], lwd = 2, pch = 20)
  points(x[idx], spikein.factors[idx], type = "b", col = clr[4], lwd = 2, pch = 20)
}

legend(
  "topright", inset = c(-0.2, 0), legend = names(clr), fill = clr, bty = "n", xpd = T
)

# Match WT DMSO factors
for(idx in list(1:4, 5:8, 9:12, 13:16)) {
  spikein.factors[idx] <- spikein.factors[idx] / spikein.factors[idx[1]] * tbl$CGIs[idx[1]]
}

plot(x, rep(0, length(x)), ylim = c(0, 1.1 * max(tbl[, 1:3], spikein.factors)), las = 2, type = "n", col = NA, border = NA, ylab = "scaling factor", main = "WT DMSO matched")

for(idx in list(1:4, 5:8, 9:12, 13:16)) {
  points(x[idx], tbl$genes[idx], type = "b", col = clr[1], lwd = 2, pch = 20)
  points(x[idx], tbl$promoters[idx], type = "b", col = clr[2], lwd = 2, pch = 20)
  points(x[idx], tbl$CGIs[idx], type = "b", col =  clr[3], lwd = 2, pch = 20)
  points(x[idx], spikein.factors[idx], type = "b", col = clr[4], lwd = 2, pch = 20)
}

legend(
  "topright", inset = c(-0.2, 0), legend = names(clr), fill = clr, bty = "n", xpd = T
)


## ----eval=FALSE----------------------------------------------------------
#  # Prepare figure layout
#  layout(matrix(1:12, 4, 3, byrow = F))
#  
#  # Plot BRD control graphs using genes, promoters and CpG-islands
#  PlotBRD(brd.gnu, title = "using genes")       # Left panels
#  PlotBRD(brd.tss, title = "using promoters")   # Center panels
#  PlotBRD(brd.cgi, title = "using CpG-islands") # Right panels

## ----echo=FALSE, message=FALSE, out.width="90%", fig.align="center"------
img <- list("TestBRD", h = 12.75, w = 9)
MkImg({
  layout(matrix(1:12, 4, 3, byrow = F))
  PlotBRD(brd.gnu, title = "using genes")       # Left panels
  PlotBRD(brd.tss, title = "using promoters")   # Center panels
  PlotBRD(brd.cgi, title = "using CpG-islands") # Right panels
})

## ------------------------------------------------------------------------
# l2c <- log2(DitherCounts(COUNTS))
# 
# chk <- NormalizeCountMatrix(l2c, brd.cgi, as.log2 = T)
# nrm <- RemoveInputBias(chk, controls = ctrl, as.log2 = T)
# 
# rmb <- RemoveInputBias(l2c, controls = ctrl, as.log2 = T)
# rmb <- NormalizeCountMatrix(rmb, brd.cgi, as.log2 = T)
# 
# layout(matrix(1:6, 2, 3, byrow = T))
# trash <- SideBySideDensity(
#   chk[WGT$category == "CGI", ], mapper = cmf, las = 2, x.labels = F, main = "chk"
# )
# trash <- SideBySideDensity(
#   nrm[WGT$category == "CGI", ], mapper = cmf, las = 2, x.labels = F, main = "nrm"
# )
# trash <- SideBySideDensity(
#   rmb[WGT$category == "CGI", ], mapper = cmf, las = 2, x.labels = F, main = "rmb"
# )


## ------------------------------------------------------------------------
# Normalize genomic read counts using BRD scaling factors from CGI bins
MakeObj(nrm, {
  nrm <- NormalizeCountMatrix(l2c, brd.cgi, as.log2 = T)
})

## ----eval=FALSE----------------------------------------------------------
#  # Show genomic read count distributions after BRD normalization
#  PlotCountDistributions(
#    nrm[, chip], pops = WGT$category, sampling = 5E5,
#    main = "Genomic bins"
#  )
#  

## ----echo=FALSE, message=FALSE, out.width="60%", fig.align="center"------
img <- list("GenomicBinsNorm", h = 6, w = 8)
MkImg({
  PlotCountDistributions(
    nrm[, chip], pops = WGT$category, sampling = 5E5,
    main = "Genomic bins"
  )
  chk <- FiniteValues(nrm[, ctrl])
  abline(h = median(nrm[chk, ctrl]), lwd = 1.5, lty = 3)
})

## ----echo=FALSE, message=FALSE, out.width="25%"--------------------------
img <- list("ScatterPlotsNorm")

grd <- "hcl.duo.light"
grp <- "intergenic"
clr <- c("grey", shades::gradient("magma", steps = nlevels(WGT$category)))
names(clr) <- c(levels(WGT$category), NA)

rng <- range(nrm[FiniteValues(nrm), ])
x <- chip[c(1, 5)]
y <- chip[c(2, 6)]
xlab <- "WT"
ylab <- "K27M"

idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    nrm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "genic"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    nrm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "TSS"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    nrm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})
grp <- "CGI"
idx <- which(WGT$category == grp)
MkImg({
  ScatterMaps(
    nrm[idx, chip], rng, x = x, y = y, bins = 705, smoothing = 25,
    gradient = grd, colors = list(m = clr[grp], p = clr[grp]),
    xlab = xlab, ylab = ylab, main = paste0(grp, "")
  )
  abline(-1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 1, 1, lwd = 1.5, col = grey(0.5, 0.4))
  abline( 0, 1, lwd = 1.5, lty = 3)
})

## ----eval=FALSE, include=FALSE-------------------------------------------
#  grp <- as.integer(1 + (WGT$category == "CGI"))
#  layout(matrix(1:4, 2, 2, byrow = T))
#  ScatterMaps(
#    l2c[, chip], x = c(1,5), y = c(2,6), pops = grp,
#    gradient = "hcl.duo.light", scales = "balanced",
#    render = "maximum",
#    bins = 705, smoothing = 25, extend = c(1, 1.1),
#    colors = list(p = clr[c(2, 4)])
#  )
#  ScatterMaps(
#    l2c[, chip], x = c(1,5), y = c(2,6), pops = grp,
#    gradient = "hcl.duo.light", scales = "balanced",
#    render = "prevalence", scoring = "glf",
#    bins = 705, smoothing = 25, extend = c(1, 1.1),
#    colors = list(p = clr[c(2, 4)])
#  )

## ----eval=FALSE----------------------------------------------------------
#  # Define column identifiers corresponding to Input profiles
#  ctrl <- c("Input_WT", "Input_K27M")
#  
#  # Preliminary application of the BRD method based on the read count matrix cnt
#  brd <- BRD(cnt, controls = ctrl)

## ----eval=FALSE----------------------------------------------------------
#  PlotBRD(brd) # Generate BRD control graphs (e.g. bivariate density distribution)

## ----eval=FALSE----------------------------------------------------------
#  # Adjusted application of the BRD method for a bimodal bivariate density
#  brd <- BRD(cnt, controls = ctrl, ncl = 2)

## ----eval=FALSE----------------------------------------------------------
#  # Apply BRD normalization factors to raw read counts
#  nrm <- t(t(cnt) * 2^brd$normfactors)

