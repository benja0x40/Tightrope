# LIBRARIES ####################################################################

library(Tightrope)

# FUNCTIONS ####################################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
ExperimentAnnotations <- function(wks_name, xps = T) {
  tbl <- get(wks_name)$ann[xps, ]
  if(any(grepl("^geo\\.", colnames(tbl)))) {
    skip <- c(
      "experiment", "submission", "study", "sample", "run", "molecule"
    )
  } else {
    skip <- c(
      "project", "SampleNumber", "IsPairedEnd", "Read1", "Read2", "NumReadsPF"
    )
  }
  skip <- c(
    skip, "mapping_command", "basemount.host", "mapping_options", "mapping_index"
  )
  rex <- "^(geo|library_|platform_|instrument_)|(file|path)"
  skip <- colnames(tbl) %in% skip | grepl(rex, colnames(tbl))
  tbl <- tbl[, ! skip]
  tbl
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
InsertAntibody <- function(viz, abd) {
  for(i in 1:length(viz)) {
    viz[[i]] <- gsub("X---X", abd, viz[[i]])
  }
  viz
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
MakeSettingsForAntibody <- function(cfg, antibody, ann) {
  for(opt in cfg$abd.prm) {
    if(is.list(cfg[[opt]])) {
      cfg[[opt]] <- cfg[[opt]][[antibody]]
    } else {
      cfg[[opt]] <- cfg[[opt]][antibody]
    }
  }
  cfg$viz <- InsertAntibody(cfg$viz, antibody)
  cfg$xps <- ann$antibody %in% c(antibody, cfg$abd.chk)
  cfg$nxp <- sum(cfg$xps)
  cfg
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
ShowDataset <- function(wks_name) {
  if(exists(wks_name)) {
    knitr::kable(ExperimentAnnotations(wks_name, xps = CFG_DTS[[wks_name]]$xps))
  }
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
ShowSettings <- function(cfg) {
  cfg <- cfg[-(1:3)]
  cfg$bdt <- paste(cfg$bdt, collapse = " >< ", sep = "")
  tbl <- matrix(unlist(cfg), 1, dimnames = list(NULL, names(cfg)))
  knitr::kable(tbl)
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
FindOutliers <- function(
  cnt, widths = NULL, minwidth = 100, maxoccur = 100, p = 0.99, plot = F
) {

  if(length(p) == 1) p <- rep(p, 2)

  if(is.null(widths)) {
    r <- rowMeans(log10(cnt + 0.5))
    outliers <- r > mean(r[rankstat(r) > p[2]])
  } else {
    w <- sort(table(widths), decreasing = T)
    w <- names(w)[1:min(which(w < maxoccur))]
    w <- widths %in% as.numeric(w) | widths < minwidth

    rx <- rankstat(log10(widths))
    ry <- rankstat(rowMeans(log10(cnt + 0.5)))
    xy <- cbind((rx + ry) / 2, (ry - rx))

    l <- sqrt(rx * (1 - ry))
    h <- sqrt(ry * (1 - rx))
    ql <- quantile(l, probs = p[1])
    qh <- quantile(h, probs = p[2])
    kl <- mean(l[l > ql])
    kh <- mean(h[h > qh])

    outliers <- w | l > kl | h > kh

    if(plot) {
      layout(matrix(1:4, 2, 2, byrow = T))
      par(pch = 20, cex = 0.5)
      plot(xy, xlim = 0:1, ylim = c(-1, 1), col = rgb(h, l, 0))
      plot(xy, xlim = 0:1, ylim = c(-1, 1), col = grey(0, 0.1))
      points(xy[l > kl, ], col = rgb(0, 1, 0), pch = 20)
      points(xy[h > kh, ], col = rgb(1, 0, 0), pch = 20)
      plot(xy, xlim = 0:1, ylim = c(-1, 1), col = grey(0, 0.1))
      points(xy[w, ], col = rgb(1, 0, 1), pch = 20)
      plot(xy[! outliers, ], xlim = 0:1, ylim = c(-1, 1), col = grey(0, 0.1))
    }
  }

  outliers
}
# INITIALIZATIONS ##############################################################

# =============================================================================.
# Dataset
# -----------------------------------------------------------------------------.
# ESC_BRD, NSC_K27M, E14_EPZ, Lu_et_al, Orlando_et_al
load_dataset   <- "Orlando_et_al" # "ALL"
# c("ESC_BRD", "NSC_K27M", "E14_EPZ", "Lu_et_al", "Orlando_et_al")
target_dataset <- "Orlando_et_al"
# -----------------------------------------------------------------------------.
# H3K27me3, H3K27me2, H3K27me1, Suz12, H3K36me3
# antibody <- "H3K27me2"
# =============================================================================.
# Glopal settings
base_path <- "/media/benjamin/USB16GB/LT_WORKS/" # mistral
base_path <- "/Volumes/USB16GB/LT_WORKS/"        # iMac
# -----------------------------------------------------------------------------.
save.results <- T
with.legend  <- T
with.axes    <- T
# -----------------------------------------------------------------------------.
smobs   <- T
dither  <- 5
zscore  <- T
knn     <- 300
mincs   <- 10
# -----------------------------------------------------------------------------.

CFG_DTS <- list()

# IMPORTATIONS #################################################################

# =============================================================================.
# ROI.LST
# -----------------------------------------------------------------------------.
ROI.LST <- readRDS(paste0(base_path, "NextSeq_E14EPZ/_RDATA_/ROI.LST.rdata"))

# =============================================================================.
# ESC_BRD
# -----------------------------------------------------------------------------.
if(any(c("ESC_BRD", "E14_EPZ", "ALL") %in% load_dataset)) {
  ESC_BRD <- list()

  data_path <- paste0(base_path, "NextSeq_ESC_BRD_R1_TEST/")
  ESC_BRD$R1 <- list(
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  data_path <- paste0(base_path, "NextSeq_ESC_BRD_R2_TEST/")
  ESC_BRD$R2 <- list(
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  # Reshape to join replicates
  ESC_BRD <- list(
    name = "ESC_BRD",
    data_path =   data_path <- paste0(base_path, "NextSeq_ESC_BRD_RDATA/"),
    ann = rbind(
      ESC_BRD$R1$ann,
      ESC_BRD$R2$ann
    ),
    ROI = readRDS(paste0(data_path, "ROI.LST.rdata")),
    CNT = with(ESC_BRD, JoinColumns(R1$CNT, R2$CNT))
  )

  if(! with(ESC_BRD, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in ESC_BRD")
  }
  # ---------------------------------------------------------------------------.
  ESC_BRD <- list2env(ESC_BRD)
}

# =============================================================================.
# NSC_K27M
# -----------------------------------------------------------------------------.
if(any(c("NSC_K27M", "ALL") %in% load_dataset)) {
  data_path <- paste0(base_path, "NextSeq_K27M_spike_in/")
  NSC_K27M <- list(
    name = "NSC_K27M",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    ROI = readRDS(paste0(data_path, "_RDATA_/ROI.LST.rdata")),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  if(! with(NSC_K27M, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in NSC_K27M")
  }
  # ---------------------------------------------------------------------------.
  NSC_K27M <- list2env(NSC_K27M)
}

# =============================================================================.
# E14_EPZ
# -----------------------------------------------------------------------------.
if(any(c("E14_EPZ", "ALL") %in% load_dataset)) {
  data_path <- paste0(base_path, "NextSeq_E14EPZ/")
  E14_EPZ <- list(
    name = "E14_EPZ",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    ROI = readRDS(paste0(data_path, "_RDATA_/ROI.LST.rdata")),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  k <- c("ESC_Input_0_R1") # "ESC_Input_0_R1", "ESC_Input_0_R2"
  E14_EPZ$CNT <- JoinColumns(E14_EPZ$CNT, ExtractColumns(ESC_BRD$CNT, k))

  k <- match(k, ESC_BRD$ann$lt_id)
  lst <- intersect(colnames(E14_EPZ$ann), colnames(ESC_BRD$ann))

  E14_EPZ$ann_orignal <- E14_EPZ$ann
  E14_EPZ$ann <- rbind(E14_EPZ$ann_orignal[, lst], ESC_BRD$ann[k, lst])

  rm(k, lst)

  if(! with(E14_EPZ, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in E14_EPZ")
  }
  # ---------------------------------------------------------------------------.
  E14_EPZ <- list2env(E14_EPZ)
}

# =============================================================================.
# Lu_et_al
# -----------------------------------------------------------------------------.
if(any(c("Lu_et_al", "ALL") %in% load_dataset)) {
  data_path <- paste0(base_path, "GSE63195_Lu_et_al_2016/")
  Lu_et_al <- list(
    name = "Lu_et_al",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    ROI = readRDS(paste0(data_path, "_RDATA_/ROI.LST.rdata")),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  if(! with(Lu_et_al, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in Lu_et_al")
  }
  # ---------------------------------------------------------------------------.
  Lu_et_al <- list2env(Lu_et_al)
}

# =============================================================================.
# Orlando_et_al
# -----------------------------------------------------------------------------.
if(any(c("Orlando_et_al", "ALL") %in% load_dataset)) {
  data_path <- paste0(base_path, "GSE60104_Orlando_et_al_2014/")
  Orlando_et_al <- list(
    name = "Orlando_et_al",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_hg38.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    ROI = readRDS(paste0(data_path, "_RDATA_/ROI.LST.rdata")),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  if(! with(Orlando_et_al, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in Orlando_et_al")
  }
  # ---------------------------------------------------------------------------.
  Orlando_et_al <- list2env(Orlando_et_al)
}

# CONFIGURATIONS ###############################################################

# =============================================================================.
# ESC_BRD
# -----------------------------------------------------------------------------.
wks_name <- "ESC_BRD"
CFG_DTS[[wks_name]] <- list(
  abd.lst = c("H3K27me3", "Suz12"),
  abd.prm = c("ref_roi", "ncl"),
  abd.chk = c("Input"),
  viz = list(
    c("ESC_Input_0_R1", "ESC_X---X_100_R1", "ESC_X---X_10_R1", "ESC_X---X_0_R1"),
    c("ESC_Input_0_R2", "ESC_X---X_100_R2", "ESC_X---X_10_R2", "ESC_X---X_0_R2"),

    c("ESC_X---X_100_R1", "ESC_X---X_10_R1"),
    c("ESC_X---X_100_R2", "ESC_X---X_10_R2"),

    c("ESC_X---X_10_R1",  "ESC_X---X_10_R2")
  ),
  layout  = matrix(1:9, 3, 3, byrow = T),
  nxp     = 0,
  rex     = "Input",
  ref_roi = c(H3K27me3 = "GN_ENS", Suz12 = "CGI_UCSC"),
  bdt     = c(0.6, 0.05),
  ncl     = c(H3K27me3 = 2, Suz12 = 1)
)

# =============================================================================.
# NSC_K27M
# -----------------------------------------------------------------------------.
wks_name <- "NSC_K27M"
CFG_DTS[[wks_name]] <- list(
  abd.lst = c("H3K27me3"),
  abd.prm = c(),
  abd.chk = c("Input"),
  viz = list(
    c("H3WT_DMSO_Input_R1",    "H3WT_DMSO_X---X_R1", "K27M_DMSO_X---X_R1", "H3WT_EPZ_X---X_R1"),
    c("H3WT_DMSO_Input_R2",    "H3WT_DMSO_X---X_R2", "K27M_DMSO_X---X_R2", "H3WT_EPZ_X---X_R2"),

    c("H3WT_DMSO_X---X_R1", "K27M_DMSO_X---X_R1"),
    c("H3WT_DMSO_X---X_R2", "K27M_DMSO_X---X_R2"),

    c("K27M_DMSO_X---X_R1", "K27M_DMSO_X---X_R2")
  ),
  layout  = matrix(1:9, 3, 3, byrow = T),
  nxp     = 0,
  rex     = "Input",
  ref_roi = "GN_ENS",
  bdt     = c(0.5, 0.05),
  ncl     = 1
)

# =============================================================================.
# E14_EPZ
# -----------------------------------------------------------------------------.
wks_name <- "E14_EPZ"
CFG_DTS[[wks_name]] <- list(
  abd.lst = c("H3K27me3", "H3K27me2", "H3K27me1", "Suz12"),
  abd.prm = c("ref_roi", "bdt", "ncl"),
  abd.chk = c("Input"),
  viz = list(
    c("X---X_E14", "X---X_7d", "X---X_4h", "X---X_8h"),
    c("X---X_E14", "X---X_16h", "X---X_24h", "X---X_48h"),
    c("X---X_E14", "X---X_96h"),
    c("ESC_Input_0_R1", "X---X_7d")
  ),
  layout  = matrix(1:9, 3, 3, byrow = T),
  nxp     = 0,
  rex     = "Input|(_7d$)",
  ref_roi = c(H3K27me3 = "GN_ENS", H3K27me2 = "CGI_UCSC", H3K27me1 = "CGI_UCSC", Suz12 = "CGI_UCSC"),
  bdt     = list(H3K27me3 = c(0.5, 0.05), H3K27me2 = c(0.5, 0.10), H3K27me1 = c(0.5, 0.10), Suz12 = c(0.5, 0.05)),
  ncl     = c(H3K27me3 = 2, H3K27me2 = 2, H3K27me1 = 2, Suz12 = 1)
)

# =============================================================================.
# Lu_et_al
# -----------------------------------------------------------------------------.
wks_name <- "Lu_et_al"
CFG_DTS[[wks_name]] <- list(
  abd.lst = c("H3K36me3", "H3K27me3"),
  abd.prm = c(),
  abd.chk = c("Input"),
  viz = list(
    c("H3WT_X---X", "K36M_X---X"),
    c("H3WT_Input", "H3WT_X---X", "K36M_X---X")
  ),
  layout  = matrix(1:9, 3, 3, byrow = T),
  nxp     = 0,
  rex     = "Input",
  ref_roi = "GN",
  bdt     = c(0.6, 0.05),
  ncl     = 2
)

# =============================================================================.
# Orlando_et_al
# -----------------------------------------------------------------------------.
wks_name <- "Orlando_et_al"
CFG_DTS[[wks_name]] <- list(
  abd.lst = c("H3K79me2", "H3K4me3"),
  abd.prm = c("ref_roi", "ncl"),
  abd.chk = c("Input"),
  viz = list(
    c("Input_0_R1", "X---X_100_R1", "X---X_25_R1", "X---X_0_R1"),
    c("Input_0_R2", "X---X_100_R2", "X---X_25_R2", "X---X_0_R2"),

    c("X---X_100_R1", "X---X_25_R1"),
    c("X---X_100_R2", "X---X_25_R2"),

    c("X---X_25_R1",  "X---X_25_R2")
  ),
  layout  = matrix(1:9, 3, 3, byrow = T),
  nxp     = 0,
  rex     = "Input",
  ref_roi = c(H3K79me2 = "GN", H3K4me3 = "GN"),
  bdt     = c(0.5, 0.05),
  ncl     = c(H3K79me2 = 1, H3K4me3 = 2)
)

# PROCESSING ###################################################################

for(dts in target_dataset) { # //////////////////////////////////////////////////

  # Dataset to be analysed
  wks <- get(dts)
  cfg <- CFG_DTS[[dts]]

  for(antibody in cfg$abd.lst) { # /////////////////////////////////////////////

    # Analysis settings
    cfg <- CFG_DTS[[dts]]
    cfg <- MakeSettingsForAntibody(cfg, antibody, wks$ann)

    # cfg$ref_roi <- "GN_ENS"
    message(dts, " | ", antibody, " | ", cfg$ref_roi)

    # Raw counts
    cnt <- with(cfg, wks$CNT[[ref_roi]][, xps])

    # Name of all experiments and controls
    experiments <- colnames(cnt)
    controls    <- experiments[grepl(cfg$rex, experiments, perl = T)]

    # BRD estimation
    brd <- with(
      cfg, BRD(
        cnt, controls = controls, smobs = smobs,
        dither = dither, zscore = zscore, knn = knn,
        bdt = bdt, ncl = ncl, mincs = mincs
      )
    )
    message("status: ", brd$status)
    print(brd$populations)

    # Results ##################################################################

    fname <- paste0(wks$data_path, cfg$ref_roi, "_", dts, "_", antibody)
    fname <- paste0(fname, ifelse(smobs, "_smobs", ""))
    if(save.results) {
      # Reshape dataset to match target experiments
      res <- with(
        wks, list(
          name        = name,
          antibody    = antibody,
          experiments = experiments,
          controls    = controls,
          CFG         = cfg,
          CNT         = ExtractColumns(CNT, experiments),
          ann         = ann[match(experiments, ann$lt_id), ],
          BRD         = c(list(ref_roi = cfg$ref_roi), brd)
        )
      )
      rownames(res$ann) <- NULL
      # Save RData
      res <- list2env(res)
      saveRDS(res, paste0(fname, ".rdata"))
      rm(res)
    }

    # Controls #################################################################

    # Figures ##################################################################

    fname <- paste0(
      base_path, "BRD_Figures/", cfg$ref_roi, "_", dts, "_", antibody
    )
    fname <- paste0(fname, ifelse(smobs, "_smobs", ""))

    # =========================================================================.
    # CDaDaDR, clustering, background ranking
    # -------------------------------------------------------------------------.
    png(
      paste0(fname, "_CDaDaDR.png"),
      width = 2 * 3, height = 2 * 3.45, units = "in", res = 300
    )
    layout(matrix(1:4, 2, 2, byrow = F))
    margins <- c(5, 4, 3, 1)
    if(! with.axes) margins <- margins / 10
    par(mar = margins, pch = 20, cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.8)

    PlotBRD(brd, with.axes = with.axes, with.legend = with.legend, res = 300)

    dev.off()

    # stop("Here")

    # =========================================================================.
    # Scaling factors
    # -------------------------------------------------------------------------.
    png(
      paste0(fname, "_Scaling.png"),
      width = dim(cfg$layout)[2] * 3, height = dim(cfg$layout)[1] * 3.45,
      units = "in", res = 200
    )

    layout(cfg$layout)
    par(
      mar = c(5, 4, 3, 1), pch = 20, cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.8
    )

    cnt <- cnt[brd$nonzero, ]

    for(cmp in cfg$viz) {
      for(lbl in cmp[-1]) {

        idx <- match(c(cmp[1], lbl), experiments)
        d <- knn_density(brd$log2counts[, idx], k = knn)
        o <- order(d)
        plot_samples(
          log2(cnt)[o, ], idx, col = colorize(d[o], mode = "rank")
        )
        plot_samples(
          log2(cnt)[brd$bg_members, ], idx,
          col = rgb(1, 0.9, 0), alpha = 0.15, add = T
        )
        plot_sources_2D(
          list(brd$bg_theta), components = idx,
          clr = rgb(1, 0, 0), lwd = 1
        )
        abline(
          a = diff(brd$bg_theta$mu[idx]), b = 1,
          col = rgb(1, 0.5, 0), lwd = 1.5
        )
      }
    }
    dev.off()
  }
  rm(cmp, lbl, idx, d, o)
}

# =============================================================================.
stop("finished")
# =============================================================================.

# PLAYGROUND ###################################################################

library(Tightrope)

# =============================================================================.
#
# -----------------------------------------------------------------------------.
knuth.welford.start <- function() {
  list(i=0, delta=0, Xmean=0, Vn=0)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
knuth.welford.iterate <- function(Xi, KW) {
  KW$i <- KW$i + 1
  KW$delta <- Xi - KW$Xmean
  KW$Xmean <- KW$Xmean + KW$delta / KW$i
  KW$Vn  <- KW$Vn + KW$delta * (Xi - KW$Xmean)
  KW
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
knuth.welford.end <- function(KW) {
  KW$Vn/(KW$i-1)
  KW
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
global_edge_list <- function(g, i) {
  e <- as_edgelist(g)
  e[, 1] <- i[e[, 1]]
  e[, 2] <- i[e[, 2]]
  e
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
increment_edges <- function(g, q, i) {
  g <- g + edges(t(global_edge_list(q, i)))
  w <- E(g)$weight
  w[is.na(w)] <- 1
  E(g)$weight <- w + count_multiple(g) - 1
  simplify(g)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
density_mode <- function(x, ...) {
  d <- density(x, ...)
  d$x[which.max(d$y)]
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
plot_sources_1D <- function(
  x, theta = NULL, bins = 200, res = bins, xlim = NULL, clr = NULL, ...
) {

  n <- length(x)

  bins <- bins - 1
  r <- diff(range(x), na.rm = T)
  brk <- 0:bins/bins * (1.2 * r) - 0.1 * r + min(x, na.rm = T)
  h.x <- hist(x, breaks = brk, plot = F)

  h.x$counts <- h.x$counts / n
  chk <- h.x$counts > 0
  h.x$counts[chk] <- with(h.x, counts[chk] / sum(diff(breaks[c(chk, T)]) * counts[chk]))


  if(is.null(xlim)) xlim = range(x)
  ylim <- range(h.x$counts)
  ylim <- min(ylim) + c(0, 1.1 * diff(ylim))
  empty.plot(xlim = xlim, ylim = ylim, yaxs = 'i', ...)

  rect(h.x$breaks[-(n+1)], 0, h.x$breaks[-1], h.x$counts, col = "black", border = "black")

  ns <- length(theta) # number of sources
  if(ns > 0) {
    if(is.null(clr)) clr <- rgb(1, 0, 0)
    if(length(clr) == 1) clr = rep(clr, ns)

    b <- SX2Y(0:res, x)
    u <- (b[-1] + b[-(res+1)])/2
    for(i in 1:ns) {
      mu    <- theta[[i]][[2]]
      sigma <- theta[[i]][[3]]
      p <- dnorm(u, mean = mu, sd = sigma)
      p <- theta[[i]][[1]] * p / sum(diff(b) * p)
      lines(u, p / 5, col = clr[i])
    }
  }
}

# =============================================================================.
# Split x in 2 random partitions
# =============================================================================.
# x <- - sodgk(0:100)
# y <- x * sodgk(0:100)
# layout(matrix(1:4, 2, 2, byrow = T))
# plot(x)
# plot(y)
# =============================================================================.
layout(matrix(1:4, 2, 2, byrow = T))
n <- 100
k <- 10
x <- SquareGrid(k = n, d = 2)$x # , subset = 5001:5100
x <- matrix(runif(10000, -1, 1), 5000, 2)
z <-  Wave2D(x, f = 3 , r = rbind(c(1, Inf))) + 2 # rep(3, nrow(x)) #
nn <- get.knn(x, k)
h <- knn_curvature(z, k, i = nn$nn.index, d = nn$nn.dist, smoothing = T)
plot(x, pch = 20, cex = 0.5, col = colorize(z, colors = "WB"))
plot(x[, 1], z, ylim = range(z, h), pch = 20, cex = 0.5)
points(x[, 1], h, pch = 20, cex = 0.5, col = "red")
clr.prm <- defineColors(c(min(h), 0, max(h)), colors = c(rgb(0, 0.5, 1), grey(0.8), rgb(1, 0.5, 0)))
clr <- makeColors(h, parameters = clr.prm)
plot(x, pch = 20, cex = 0.5, col = clr)
# =============================================================================.
n <- 2000
x <- SquareGrid(k = 2, d = 2)$x
x <- ClonalGaussian(n = n, mu = x, sigma = 1/30)$x
x <- rbind(
  matrix(runif(2 * n, min(x), max(x)), n),
  matrix(runif(2 * n, min(x), 0), n),
  matrix(runif(2 * n, 0, max(x)), n),
  x
)
n <- nrow(x)

layout(matrix(1:4, 2, 2, byrow = T))
k <- 100
nn <- get.knn(x, k)
p <- knn_density(x, k, i = nn$nn.index, d = nn$nn.dist, smoothing = T)
h <- knn_curvature(p, k, i = nn$nn.index, d = nn$nn.dist, smoothing = T)

clr.prm <- defineColors(c(min(h), 0, max(h)), colors = c(rgb(0, 0.5, 1), grey(0.8), rgb(1, 0.5, 0)))
clr <- makeColors(h, parameters = clr.prm)
plot(x, pch = 20, cex = 0.5, col = colorize(p, mode = "01", col = "ry"))
plot(x, pch = 20, cex = 0.5, col = colorize(p, mode = "01", col = "ry"))
plot(p, h, pch = 20, cex = 0.5, col = grey(0, 0.2))

Q <-  knuth.welford.start()

G <- graph.empty(n)
G <- set_edge_attr(G, "weight",   value = 0)

layout(matrix(1:4, 2, 2, byrow = T))
rnit <- 10
d <- rep(0, n)
pb <- txtProgressBar(min = 1, max = rnit, char = "|", style= 3)
for(r in 1:rnit) {
  chk <- sample(c(T, F), n, replace = T)
  # First half
  i <- which(chk)
  d[i] <- knn_density(x[i, ], k / 2, smoothing = T)
  qs <- QuickShift(x[i, ], d[i])
  G <- increment_edges(G, qs, i)
  # Second half
  i <- which(! chk)
  d[i] <- knn_density(x[i, ], k / 2, smoothing = T)
  qs <- QuickShift(x[i, ], d[i])
  G <- increment_edges(G, qs, i)
  # Online computation of mean+variance
  Q <- knuth.welford.iterate(2 * d, Q)

  setTxtProgressBar(pb, r)
}
close(pb)

Q <- knuth.welford.end(Q)

e <- as_edgelist(G)
d <- sqrt(rowSums((x[e[, 1], ] - x[e[, 2], ])^2))
w <- E(G)$weight
clr <- colorize(w, mode = "rank")
ColorChannel(clr, "a") <- 0.2
# clr[d > 0.4] <- rgb(1 ,0, 0)
# G <- G - E(G)[d > 0.4]

layout(matrix(1:4, 2, 2, byrow = T))
plot(x, pch = 20, cex = 0.5, col = colorize(p, mode = "rank", col = "ry"))
# plot(x, pch = 20, cex = 0.5, col = colorize(Q$Xmean, mode = "rank", col = "ry"))
# plot(x, pch = 20, cex = 0.5, col = colorize(Q$Vn, mode = "rank", col = "ry"))
# PlotQuickShift(x, G, col = clr)
v <- sqrt(Q$Vn)
# theta <- mv_gmm(v, ns = 2)$theta
# theta <- list(list(1, mu = density_mode(v), sigma = mad(v)))
plot_sources_1D(v, theta)
abline(v = density_mode(v), col = "red")
plot(x[, 1], p, pch = 20, cex = 0.5, col = grey(0, 0.5))
z <- p + rnorm(length(p), sd = density_mode(v))
plot(x[, 1], z, pch = 20, cex = 0.5, col = grey(0, 0.5))
# -----------------------------------------------------------------------------.
# plot(density(log10(d)))
# plot(d, w, log = "xy")
# -----------------------------------------------------------------------------.
# plot(p, Q$Xmean, log = "xy")
# o <- order(Q$Vn)
# plot(Q$Vn[o], cumsum(rankstat(Q$Vn)[o]), log = "x")
# -----------------------------------------------------------------------------.
# g <- QuickShift(x, p)
# e <- as_edgelist(g)
# v.a <- e[, 1]
# v.b <- e[, 2]
# d <- E(g)$distance
# z <- p[v.b] - p[v.a]
#
# q <- data.frame(
#   matrix(
#     F, n, 5, dimnames = list(
#       NULL, c("root", "node", "parent", "distance", "gradient")
#     )
#   )
# )
# q[v.a, 3] <- v.b
# q[v.a, 4] <- d
# q[v.a, 5] <- z
# q$node[v.b] <- T
# q$root <- q$parent == 0
# -----------------------------------------------------------------------------.
# layout(matrix(1:4, 2, 2, byrow = T))
# plot(x, pch = 20, cex = 0.5, col = colorize(p, mode = "rank", col = "ry"))
# with(q, plot(log10(cbind(d, z)), pch = 20, cex = 0.5))
# clr <- colorize(d, mode = "rank")
# clr[log10(d) > -0.5] <- rgb(1 ,0, 0)
# PlotQuickShift(x, g, col = clr)
# =============================================================================.
# g <- QuickShift(x, p)
# e <- as_edgelist(g)
# v.a <- e[, 1]
# v.b <- e[, 2]
# q <- matrix(0, n, 4)
# q[v.b, 1] <- 1
# q[v.a, 2] <- v.b
#
# blur <- rep(0, n)
# pb <- txtProgressBar(min = 1, max = rnit, char = "|", style= 3)
# for(r in 1:rnit) {
#   i <- which(sample(c(T, F), n, replace = T))
#   g <- QuickShift(x[i, ], p[i])
#   e <- as_edgelist(g)
#   v.a <- i[e[, 1]]
#   v.b <- i[e[, 2]]
#   q[v.a, 2] <- q[v.a, 2] + (q[v.a, 1] ==  v.b)
#   q[v.a, 3] <- q[v.a, 3] + 1
#   setTxtProgressBar(pb, r)
# }
# close(pb)
# # hist(q[, 2] / q[, 3], breaks = 30, col = "black")
#
# layout(matrix(1:4, 2, 2, byrow = T))
# plot(x, pch = 20, cex = 0.5, col = grey(0, 0.1))
# plot(x, pch = 20, cex = 0.5, col = colorize(p, mode = "rank", col = "ry"))
# plot(x[q[, 2] / q[, 3] > 0.4, ], pch = 20, cex = 0.5, col = grey(0, 0.5))
#
# for(r in 1:4) {
#   i <- which(sample(c(T, F), n, replace = T))
#   g <- QuickShift(x[i, ], p[i])
#   e <- as_edgelist(g)
#   v.a <- i[e[, 1]]
#   v.b <- i[e[, 2]]
#   q[v.a, 2] <- q[v.a, 2] + (q[v.a, 1] ==  v.b)
#   q[v.a, 3] <- q[v.a, 3] + 1
#   setTxtProgressBar(pb, r)
# }
#
# g <- QuickShift(x, p)

# d <- knn.dist(x, k = 1)
# plot(density(d))

# g <- QuickShift(x, p)
# d <- log10(E(g)$distance)
#
# f <- density(d)
# u <- f$x
# y <- f$y / max(f$y)
#
# f <- mv_mle(d)
# f$sigma <- mad(d)
# z <- with(f, dnorm(u, mu, sigma))
# z <- z / max(z)
#
# plot(u, y, type = "l")
# lines(u, z, col = "red")
#
# p <- pnorm(abs(d - f$mu), 0, f$sigma, lower.tail = F)
#
# PlotQuickShift(x, g, col = ifelse(p < 1E-2, rgb(1, 0, 0, 0.2), grey(0, 0.2)))
