## ----hidden_dependencies, include=FALSE----------------------------------
library(knitr)
library(data.table)
library(Barbouille)
library(Tightrope)

data("EGA91_human") 
data("CGI_hg38") 
data("hg38_blacklist") 

data("Orlando_METADATA") 
data("Orlando_COUNTS") 
data("Orlando_SPIKEIN")

data("EGA91_mouse") 
data("CGI_mm10")
data("mm10_blacklist") 

data("ESC_BRD_METADATA") 
data("ESC_BRD_COUNTS")
data("ESC_BRD_SPIKEIN")

# data("NSC_K27M_METADATA") 
# data("NSC_K27M_COUNTS")


## ------------------------------------------------------------------------
library(Tightrope)

META   <- Orlando_METADATA 
COUNTS <- Orlando_COUNTS

META   <- ESC_BRD_METADATA 
COUNTS <- ESC_BRD_COUNTS

## ----show_dataset, echo=FALSE, out.width="50%"---------------------------
# kable(META, caption = "ChIP-seq experiments")

## ----fig.height=3.1, fig.width=3.1---------------------------------------
CGI <- resize(CGI_hg38, width = 4000, fix = "center", use.names = F)
blacklist <- hg38_blacklist

CGI <- resize(CGI_mm10, width = 4000, fix = "center", use.names = F)
blacklist <- mm10_blacklist

chip <- META[antibody == "H3K79me2"]$sample_id
ctrl <- META[antibody == "Input"]$sample_id

chip <- META[antibody == "H3K27me3" & replicate == "R2"]$sample_id
ctrl <- META[antibody == "Input" & replicate == "R2"]$sample_id

ref <- COUNTS$GNU[, c(chip, ctrl)]
chk <- FiniteValues(log2(ref))
brd <- BRD(ref[chk, ], controls = ctrl, ncl = 2, bdt = c(0.2, 0.03))

chk <- countOverlaps(CGI, blacklist) == 0
cnt <- COUNTS$CGI_UCSC[chk, c(chip, ctrl)]
l2c <- log2(DitherCounts(cnt))

nrm <- t(t(l2c) + brd$normfactors)

## ----fig.height=6.5, out.width="30%"-------------------------------------
# Adjust graphic options
par(
  pch = 20, mar = c(4.5, 4.5, 1, 1), cex = 1.5 # cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
)

# Show the BRD results
PlotBRD(brd)

## ----count_distribution, eval=FALSE--------------------------------------
#  # Create a color mapping function to represent read count distributions
#  cmf <- function(k) colorize(k, mode = "01", col = "Wry")
#  
#  # Adjust plot margins to show the name of experimental conditions
#  par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)
#  
#  # Show read count distributions after spike-in and BRD normalizations
#  r <- SideBySideDensity(
#    l2c, nx = ncol(cnt) * 25, ny = 150, method = "ash", spacing = 0.4,
#    mapper = cmf, las = 2, main = "spike-in"
#  )
#  r <- SideBySideDensity(
#    nrm, nx = ncol(cnt) * 25, ny = 150, method = "ash", spacing = 0.4,
#    mapper = cmf, las = 2, main = "BRD"
#  )

## ----show_count_distribution, echo=FALSE, fig.height=5.5, out.width="40%", ref.label='count_distribution'----
# Create a color mapping function to represent read count distributions 
cmf <- function(k) colorize(k, mode = "01", col = "Wry")

# Adjust plot margins to show the name of experimental conditions
par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1) 

# Show read count distributions after spike-in and BRD normalizations
r <- SideBySideDensity(
  l2c, nx = ncol(cnt) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "spike-in"
)
r <- SideBySideDensity(
  nrm, nx = ncol(cnt) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "BRD"
)

