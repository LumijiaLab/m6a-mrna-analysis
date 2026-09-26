suppressPackageStartupMessages({
  library(data.table)
})


# --- Parameters ---

DIR_MATERIAL <- "material"
DIR_PROCESS  <- "process"

IN_DEG     <- file.path(DIR_PROCESS,  "deg_iDC_hMPV_vs_Ctrl.csv")
IN_M6A     <- file.path(DIR_PROCESS,  "m6a_diff_iDC_hMPV_vs_Ctrl.csv")
IN_REGION  <- file.path(DIR_MATERIAL, "m6a_peak_region_hMPV.csv")

OUT_POINTS <- file.path(DIR_PROCESS, "fig1E_points_hMPV.csv")

THR_FC <- 1              # |log2FC| cutoff for a significant gene
THR_P  <- 0.05
P_COL  <- "PValue"       # nominal P, not FDR; same rule as Figure 1A

REGION_LEVELS <- c("5'UTR", "StartC", "CDS", "StopC", "3'UTR")

DO_SELFCHECK <- TRUE     # print per-region Pearson r to the console


# --- Functions ---

save_csv <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, bom = TRUE, na = "")
  cat(sprintf("    OK  %s  (%d rows x %d columns)\n", basename(path), nrow(dt), ncol(dt)))
}


# --- Main ---

for (f in c(IN_DEG, IN_M6A, IN_REGION))
  if (!file.exists(f))
    stop("Input file not found: ", f, "\n  Run 01_normalize_inputs.R first")

deg    <- fread(IN_DEG,    encoding = "UTF-8")
m6a    <- fread(IN_M6A,    encoding = "UTF-8")
region <- fread(IN_REGION, encoding = "UTF-8")

deg[, x := fifelse(!is.na(get(P_COL)) & get(P_COL) < THR_P & abs(logFC) >= THR_FC,
                   logFC, 0)]
xmap <- deg[!is.na(gene_name) & nzchar(gene_name),
            .(x = x[which.max(abs(x))]), by = .(gene = toupper(gene_name))]
cat(sprintf("    x: %d genes in the DEG table (%d assigned a non-zero value)\n",
            nrow(xmap), xmap[x != 0, .N]))

pk <- m6a[, .(PeakID, Foldchange, Regulation)][
        region[, .(PeakID, gene = toupper(GeneName), region = Location_std)],
        on = "PeakID", nomatch = 0L]
pk[, signed_fc := fifelse(Regulation == "up", Foldchange, -Foldchange)]
gp <- pk[, .(sum_fc = sum(signed_fc), n_peak = .N), by = .(gene, region)]
gp[, y := fifelse(sum_fc == 0, 0, sign(sum_fc) * log2(abs(sum_fc)))]
cat(sprintf("    y: %d region x gene combinations (from %d peaks)\n", nrow(gp), nrow(pk)))

pts <- merge(gp, xmap, by = "gene", all.x = TRUE)
pts[is.na(x), x := 0]
cat(sprintf("    %d points kept\n", nrow(pts)))
cat(sprintf("    x = 0 in %d points (%.1f%%)\n", pts[x == 0, .N], 100 * pts[x == 0, .N] / nrow(pts)))

pts[, region := factor(region, levels = REGION_LEVELS)]
setorder(pts, region, -y)
pts[, gene := toupper(gene)]

cat("    points per region:\n")
print(pts[, .(n = .N, x_nonzero = sum(x != 0)), by = region])

if (DO_SELFCHECK) {
  cat("    Pearson r per region:\n")
  for (rg in REGION_LEVELS) {
    sub <- pts[region == rg]
    cat(sprintf("      %-6s n = %4d   r = %.4f\n", rg, nrow(sub), cor(sub$x, sub$y)))
  }
}

save_csv(pts[, .(gene, region, x, y, sum_fc, n_peak)], OUT_POINTS)
