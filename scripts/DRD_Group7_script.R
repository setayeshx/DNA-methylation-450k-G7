##############################################################################
## GROUP 7 - DNA Methylation Analysis Pipeline
## Course: DNA/RNA Dynamics - University of Bologna (2025-2026)
##
## Group 7 assignments:
##   Step 3  - Address:          54740449   (AddressB of cg21802055; Type I, Red)
##   Step 5  - detP threshold:   0.05
##   Step 7  - Normalization:    preprocessQuantile
##   Step 9  - Statistical test: t-test
##
## NOTE: This is the packaged pipeline. Set the working directory below to the
## repo root, and make sure Input_data/ and Illumina450Manifest_clean.RData
## are present. Verify against your finalised local copy before relying on it.
##############################################################################

##############################################################################
## STEP 0: Environment and libraries
##############################################################################
rm(list = ls())
setwd("~/your_working_directory")          # <-- set to the repo root
suppressMessages(library(minfi))           # import, QC, normalization, beta/M values
library(gplots)                            # heatmap.2
library(qqman)                             # Manhattan plot

##############################################################################
## STEP 1: Load raw data -> RGChannelSet (RGset)
##############################################################################
list.files("Input_data/")
SampleSheet <- read.csv("Input_data/SampleSheet_Report_II.csv", header = TRUE)

baseDir <- "Input_data"
targets <- read.metharray.sheet(baseDir)

RGset <- read.metharray.exp(targets = targets)
RGset

##############################################################################
## STEP 2: Extract Red and Green fluorescence
##############################################################################
Red   <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))
dim(Red); dim(Green)

##############################################################################
## STEP 3: Address 54740449 - fluorescence and probe type
##############################################################################
Red[rownames(Red)     == "54740449", ]
Green[rownames(Green) == "54740449", ]

load("Illumina450Manifest_clean.RData")
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == "54740449", ]
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID == "54740449", ]
# -> AddressB of cg21802055; Type I probe, Red channel.

##############################################################################
## STEP 4: Create MSet.raw and quality control
##############################################################################
MSet.raw <- preprocessRaw(RGset)
MSet.raw

qc <- getQC(MSet.raw)
plotQC(qc)                                 # -> outputs/qc_plot.png
# Conservative cutoff 10.5 flags 5/8 samples; all retained (detP failure < 0.6%).

controlStripPlot(RGset, controls = "NEGATIVE")   # -> outputs/controlStripPlot.png

##############################################################################
## STEP 5: Detection p-value filtering (threshold 0.05)
##############################################################################
detP <- detectionP(RGset)
failed <- detP > 0.05
colMeans(failed)                           # fraction of failed positions per sample
summary(failed)

##############################################################################
## STEP 6: Beta and M values + group-wise distributions
##############################################################################
beta <- getBeta(MSet.raw)
M    <- getM(MSet.raw)

## Defensive row-ordering: align the sample sheet to the columns of the matrices
## via match() so the CTRL/DIS labelling can never silently mis-map.
pheno <- SampleSheet
ord   <- match(colnames(beta), pheno$Basename)        # adjust key column if needed
stopifnot(!any(is.na(ord)))
pheno <- pheno[ord, ]
grp   <- factor(pheno$Group)                          # CTRL / DIS

beta_CTRL <- beta[, grp == "CTRL"]
beta_DIS  <- beta[, grp == "DIS"]
M_CTRL    <- M[,    grp == "CTRL"]
M_DIS     <- M[,    grp == "DIS"]

mean_beta_CTRL <- apply(beta_CTRL, 1, mean, na.rm = TRUE)
mean_beta_DIS  <- apply(beta_DIS,  1, mean, na.rm = TRUE)
mean_M_CTRL    <- apply(M_CTRL,    1, mean, na.rm = TRUE)
mean_M_DIS     <- apply(M_DIS,     1, mean, na.rm = TRUE)

## -> outputs/beta_m_values.png : density of Beta and M, CTRL vs DIS

##############################################################################
## STEP 7: Normalization (preprocessQuantile) + raw vs normalised comparison
##############################################################################
preprocessQuantile_results <- preprocessQuantile(RGset)
beta_norm <- getBeta(preprocessQuantile_results)

## Raw vs normalised comparison (mean / SD / density / boxplot)
## -> outputs/Raw_normalised_beta.png

## OPTIONAL: group-coloured boxplot of normalised beta values.
## Appropriate here because it lets you eyeball per-sample distribution shifts
## by group after normalisation; colours map to CTRL/DIS.
group_colors <- ifelse(grp == "CTRL", "orange", "purple")
boxplot(beta_norm, col = group_colors, names = pheno$Group,
        las = 2, main = "Normalised beta by sample (group-coloured)")

##############################################################################
## STEP 8: Principal Component Analysis
##############################################################################
pca <- prcomp(t(na.omit(beta_norm)), scale. = TRUE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
print(round(var_explained[1:5] * 100, 2))            # % variance, PC1-PC5

## -> outputs/PCA_groups.png, outputs/PCA_sex.png, outputs/scree_plot.png
## Result: CTRL and DIS do not separate globally (localised differences).

##############################################################################
## STEP 9: Differential methylation - t-test
##############################################################################
My_ttest <- function(x) {
  t.test(x[grp == "CTRL"], x[grp == "DIS"])$p.value
}
pValues <- apply(beta_norm, 1, My_ttest)

stopifnot(length(pValues) == nrow(beta_norm))        # hard guard

sum(pValues < 0.05, na.rm = TRUE)                    # 16,107 nominal DMPs

##############################################################################
## STEP 10: Multiple-testing correction (BH, Bonferroni)
##############################################################################
p_BH   <- p.adjust(pValues, method = "BH")
p_Bonf <- p.adjust(pValues, method = "bonferroni")

sum(p_BH   < 0.05, na.rm = TRUE)                     # 0
sum(p_Bonf < 0.05, na.rm = TRUE)                     # 0

## -> outputs/Histogram_pvalues.png, outputs/Boxplot_corrections.png

##############################################################################
## STEP 11: Manhattan and Volcano plots
##############################################################################
delta <- mean_beta_DIS - mean_beta_CTRL              # effect size (delta beta)

## Build a results table; annotate with chromosome + position from the manifest.
res <- data.frame(
  probe = rownames(beta_norm),
  p     = pValues,
  delta = delta[rownames(beta_norm)]
)

## NA filter BEFORE manhattan() — required, or the call errors on missing rows.
res_man <- merge(res, Illumina450Manifest_clean[, c("IlmnID", "CHR", "MAPINFO")],
                 by.x = "probe", by.y = "IlmnID")
res_man <- res_man[!is.na(res_man$p) &
                   !is.na(res_man$CHR) &
                   !is.na(res_man$MAPINFO), ]
res_man$CHR <- as.numeric(as.character(res_man$CHR))
res_man <- res_man[!is.na(res_man$CHR), ]

manhattan(res_man, chr = "CHR", bp = "MAPINFO", p = "p", snp = "probe")
## -> outputs/manhattan_plot.png

## Volcano: delta beta vs -log10(p)
plot(res$delta, -log10(res$p), pch = 16, cex = 0.4,
     xlab = "delta beta", ylab = "-log10(p)", main = "Volcano plot")
abline(v = c(-0.1, 0.1), lty = 2)
## -> outputs/volcano_plot.png  (candidates with |delta beta| > 0.1)

##############################################################################
## STEP 12: Heatmap of top 100 probes (hierarchical clustering)
##############################################################################
top100 <- names(sort(pValues))[1:100]
mat    <- as.matrix(beta_norm[top100, ])

colside <- ifelse(grp == "CTRL", "orange", "purple")

heatmap.2(mat, col = colorRampPalette(c("blue", "white", "red"))(100),
          ColSideColors = colside, trace = "none",
          hclustfun = function(d) hclust(d, method = "complete"),
          main = "Top 100 probes - complete linkage")
## -> outputs/Complete_linkage_heatmap.png
## Repeat with method = "average" and "single" for the other two heatmaps.
## NB: top-100 selected on the group contrast -> separation is illustrative
## (selection bias), not independent evidence of group structure.

##############################################################################
## END OF PIPELINE
##############################################################################
