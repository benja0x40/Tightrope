# =============================================================================.
# Load Tightrope package (this will produce numerous messages)
# -----------------------------------------------------------------------------.
library(Tightrope)

# =============================================================================.
# Create a color mapping function to represent read count distributions
# -----------------------------------------------------------------------------.
cmf <- function(k) colorize(k, mode = "01", col = "Wry")

# =============================================================================.
# Load genome annotations
# -----------------------------------------------------------------------------.
# Human genome annotations (Ensembl GRCh38.p10 release 91)
data("EGA91_human")    # Gene features from Ensembl
data("CGI_hg38")       # CpG Islands from UCSC
data("hg38_blacklist") # Blacklisted regions from ENCODE

# Mouse genome annotations (Ensembl GRCm38.p5 release 91)
data("EGA91_mouse")    # Gene features from Ensembl
data("CGI_mm10")       # CpG Islands from UCSC
data("mm10_blacklist") # Blacklisted regions from ENCODE

# =============================================================================.
# Data from Orlando et al. 2014
# -----------------------------------------------------------------------------.
data("Orlando_METADATA") # Annotation of experimental conditions
data("Orlando_COUNTS")   # Precomputed read count matrixes
data("Orlando_SPIKEIN")  # Number of Drosophila reads (spike-in)

# =============================================================================.
# H3K79me2 (Orlando et al. 2014)
# -----------------------------------------------------------------------------.
message("H3K79me2 (Orlando et al. 2014)")

chip <- Orlando_METADATA[antibody == "H3K79me2" & replicate == "2"]$sample_id
ctrl <- Orlando_METADATA[antibody == "Input" & replicate == "2"]$sample_id

cat("List of ChIP samples:", paste(chip, collapse = ", "), sep = "\n")
cat("List of Input samples:", paste(ctrl, collapse = ", "), sep = "\n")

# Use read counts over genes as reference for BRD normalization
ref <- Orlando_COUNTS$GNU[, c(chip, ctrl)]

# Search for background candidates and estimate normalization factors
brd <- BRD(ref, controls = ctrl, ncl = 1, bdt = c(0.2, 0.03))

# Plot BRD control graphs
par(pch = 20, mar = c(4.5, 4.5, 1, 1), cex = 1.5) # Adjust graphic options
PlotBRD(brd)                                      # Show the BRD results

# Extract read counts over TSS
cnt <- Orlando_COUNTS$TSS[, c(chip, ctrl)]

# Apply dithering and log2 transformation
l2c <- log2(DitherCounts(cnt))

# Apply BRD normalization factors
nrm.brd <- t(t(l2c) + brd$normfactors[c(chip, ctrl)])

# Apply spike-in normalization
spk <- Orlando_METADATA$Drosophila.Unique.Reads
names(spk) <- Orlando_METADATA$sample_id
spk <- spk[c(chip, ctrl)]             # Extract Drosophila read counts
spk <-  log2(mean(spk) / spk)         # Compute normalization factors
nrm.spk <- t(t(l2c) + spk)            # Apply normalization factors

# Adjust plot margins to show the name of experimental conditions
par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)

# Show read count distributions after spike-in and BRD normalizations
layout(matrix(1:2, 2, 1, byrow = T)) # Split plotting device (2 panels)
r <- SideBySideDensity(
  nrm.spk, nx = ncol(nrm.spk) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "spike-in"
)
r <- SideBySideDensity(
  nrm.brd, nx = ncol(nrm.brd) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "BRD"
)

# =============================================================================.
# H3K4me3 (Orlando et al. 2014)
# -----------------------------------------------------------------------------.
message("H3K4me3 (Orlando et al. 2014)")

chip <- Orlando_METADATA[antibody == "H3K4me3"]$sample_id
ctrl <- Orlando_METADATA[antibody == "Input"]$sample_id

# Use read counts over TSS as reference for BRD normalization
ref <- Orlando_COUNTS$TSS[, c(chip, ctrl)]

# Search for background candidates and estimate normalization factors
brd <- BRD(ref, controls = ctrl, ncl = 2, bdt = c(0.2, 0.03))

# Plot BRD control graphs
par(pch = 20, mar = c(4.5, 4.5, 1, 1), cex = 1.5) # Adjust graphic options
PlotBRD(brd)                                      # Show the BRD results

# Extract read counts over TSS
cnt <- Orlando_COUNTS$TSS[, c(chip, ctrl)]

# Apply dithering and log2 transformation
l2c <- log2(DitherCounts(cnt))

# Apply BRD normalization factors
nrm.brd <- t(t(l2c) + brd$normfactors[c(chip, ctrl)])

# Apply spike-in normalization
spk <- Orlando_SPIKEIN[c(chip, ctrl)] # Extract Drosophila read counts
spk <-  log2(mean(spk) / spk)         # Compute normalization factors
nrm.spk <- t(t(l2c) + spk)            # Apply normalization factors

# Adjust plot margins to show the name of experimental conditions
par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)

# Show read count distributions after spike-in and BRD normalizations
layout(matrix(1:2, 2, 1, byrow = T)) # Split plotting device (2 panels)
r <- SideBySideDensity(
  nrm.spk, nx = ncol(nrm.spk) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "spike-in"
)
r <- SideBySideDensity(
  nrm.brd, nx = ncol(nrm.brd) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "BRD"
)

# =============================================================================.
# H3K27me3 (ESC)
# -----------------------------------------------------------------------------.
message("H3K27me3 (ESC)")

data("ESC_BRD_METADATA") # Annotation of experimental conditions
data("ESC_BRD_COUNTS")   # Precomputed read count matrixes
data("ESC_BRD_SPIKEIN")  # Number of Drosophila reads (spike-in)

chip <- ESC_BRD_METADATA[antibody == "H3K27me3" & replicate == "R2"]$sample_id
ctrl <- ESC_BRD_METADATA[antibody == "Input" & replicate == "R2"]$sample_id

# Use read counts over gene units as reference for BRD normalization
ref <- ESC_BRD_COUNTS$GNU[, c(chip, ctrl)]

# Search for background candidates and estimate normalization factors
brd <- BRD(ref, controls = ctrl, ncl = 2, bdt = c(0.2, 0.03))

# Plot BRD control graphs
par(pch = 20, mar = c(4.5, 4.5, 1, 1), cex = 1.5) # Adjust graphic options
PlotBRD(brd)                                      # Show the BRD results

# Extract read counts over CGIs
cnt <- ESC_BRD_COUNTS$CGI_UCSC[, c(chip, ctrl)]

# Apply dithering and log2 transformation
l2c <- log2(DitherCounts(cnt))

# Apply BRD normalization factors
nrm.brd <- t(t(l2c) + brd$normfactors[c(chip, ctrl)])

# Apply spike-in normalization
spk <- ESC_BRD_SPIKEIN[c(chip, ctrl)] # Extract Drosophila read counts
spk <-  log2(mean(spk) / spk)         # Compute normalization factors
nrm.spk <- t(t(l2c) + spk)            # Apply normalization factors

# Adjust plot margins to show the name of experimental conditions
par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)

# Show read count distributions after spike-in and BRD normalizations
layout(matrix(1:2, 2, 1, byrow = T)) # Split plotting device (2 panels)
r <- SideBySideDensity(
  nrm.spk, nx = ncol(nrm.spk) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "spike-in"
)
r <- SideBySideDensity(
  nrm.brd, nx = ncol(nrm.brd) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "BRD"
)

# =============================================================================.
# H3K27me3 (NSC)
# -----------------------------------------------------------------------------.
message("H3K27me3 (NSC)")

data("NSC_K27M_METADATA") # Annotation of experimental conditions
data("NSC_K27M_COUNTS")   # Precomputed read count matrixes
data("NSC_K27M_SPIKEIN")  # Number of Drosophila reads (spike-in)

chip <- NSC_K27M_METADATA[antibody == "H3K27me3"]$sample_id
ctrl <- NSC_K27M_METADATA[antibody == "Input"]$sample_id

# Use read counts over gene units as reference for BRD normalization
ref <- NSC_K27M_COUNTS$GNU[, c(chip, ctrl)]

# Search for background candidates and estimate normalization factors
brd <- BRD(ref, controls = ctrl, ncl = 1, bdt = c(0.2, 0.03))

# Plot BRD control graphs
par(pch = 20, mar = c(4.5, 4.5, 1, 1), cex = 1.5) # Adjust graphic options
PlotBRD(brd)                                      # Show the BRD results

# Extract read counts over CGIs
cnt <- NSC_K27M_COUNTS$CGI_UCSC[, c(chip, ctrl)]

# Apply dithering and log2 transformation
l2c <- log2(DitherCounts(cnt))

# Apply BRD normalization factors
nrm.brd <- t(t(l2c) + brd$normfactors[c(chip, ctrl)])

# Apply spike-in normalization
spk <- NSC_K27M_SPIKEIN[c(chip, ctrl)] # Extract Drosophila read counts
spk <-  log2(mean(spk) / spk)          # Compute normalization factors
nrm.spk <- t(t(l2c) + spk)             # Apply normalization factors

# Adjust plot margins to show the name of experimental conditions
par(mar = c(12, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)

# Show read count distributions after spike-in and BRD normalizations
layout(matrix(1:2, 2, 1, byrow = T)) # Split plotting device (2 panels)
r <- SideBySideDensity(
  nrm.spk, nx = ncol(nrm.spk) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "spike-in"
)
r <- SideBySideDensity(
  nrm.brd, nx = ncol(nrm.brd) * 25, ny = 150, method = "ash", spacing = 0.4,
  mapper = cmf, las = 2, main = "BRD"
)
