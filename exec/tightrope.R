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
# ESC_BRD, NSC_K27M, E14_EPZ, Lu_et_al
load_dataset   <- "E14_EPZ" # "ALL"
target_dataset <- "E14_EPZ" # c("ESC_BRD", "NSC_K27M", "E14_EPZ", "Lu_et_al")
# -----------------------------------------------------------------------------.
# H3K27me3, H3K27me2, H3K27me1, Suz12, H3K36me3
# antibody <- "H3K27me2"
# =============================================================================.
# Glopal settings
base_path <- "/media/benjamin/USB16GB/LT_WORKS/" # mistral
base_path <- "/Volumes/USB16GB/LT_WORKS/"        # iMac
# -----------------------------------------------------------------------------.
save.results <- F
with.legend  <- T
with.axes    <- T
# -----------------------------------------------------------------------------.
smobs   <- T
dither  <- 5
npc     <- 2
zscore  <- T
knn     <- 300
rare    <- 0.01
method  <- "ica"
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

  data_path <- "/Volumes/USB16GB/LT_WORKS/NextSeq_ESC_BRD_R1_TEST/"
  ESC_BRD$R1 <- list(
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  data_path <- "/Volumes/USB16GB/LT_WORKS/NextSeq_ESC_BRD_R2_TEST/"
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
    data_path = "/Volumes/USB16GB/LT_WORKS/NextSeq_ESC_BRD_RDATA/",
    ann = rbind(
      ESC_BRD$R1$ann,
      ESC_BRD$R2$ann
    ),
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
  data_path <- "/Volumes/USB16GB/LT_WORKS/NextSeq_K27M_spike_in/"
  NSC_K27M <- list(
    name = "NSC_K27M",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
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
  data_path <- "/Volumes/USB16GB/LT_WORKS/NextSeq_E14EPZ/"
  E14_EPZ <- list(
    name = "E14_EPZ",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
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
  data_path <- "/Volumes/USB16GB/LT_WORKS/GSE63195_Lu_et_al_2016/"
  Lu_et_al <- list(
    name = "Lu_et_al",
    data_path = paste0(data_path, "_RDATA_/"),
    ann = read.delim(
      paste0(data_path, "_LittleThumb_/datasets/bowtie2_mm10.txt"),
      comment.char = "#", stringsAsFactors = F
    ),
    CNT = readRDS(paste0(data_path, "_RDATA_/ROI.CNT.rdata"))
  )
  rm(data_path)

  if(! with(Lu_et_al, all(ann$lt_id == colnames(CNT$GN)))) {
    stop("Inconsistent annotation of samples in Lu_et_al")
  }
  # ---------------------------------------------------------------------------.
  Lu_et_al <- list2env(Lu_et_al)
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
        cnt, controls = controls, smobs = smobs, dither = dither,
        npc = npc, zscore = zscore, knn = knn, rare = rare, method = method,
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
