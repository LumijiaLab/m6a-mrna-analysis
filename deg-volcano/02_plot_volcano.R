suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})


# --- Parameters ---

DIR_PROCESS <- "process"
DIR_RESULTS <- "results"

IN_CSV     <- file.path(DIR_PROCESS, "deg_iDC_hMPV_vs_Ctrl.csv")
OUT_PREFIX <- "fig1A_volcano_iDC"

THR_FC <- 1
THR_P  <- 0.05
P_COL  <- "PValue"       # nominal P, not FDR

LABEL_GENES <- c("IFIT2", "IFIT3", "IFIH1", "CD274", "IFNB1", "NFKB1", "CD80", "CD63", "DDX58", "IL6", "IFITM1", "TLR8", "TLR7", "IFNL1", "CD86", "IRF7", "IFITM2", "IFNL3", "IFITM3", "TBK1", "HLA-DOA1", "IFNL2", "IRF3", "MAVS", "HLA-B", "HLA-C", "TNF")

NO_LABEL <- FALSE        # TRUE = no gene labels
FIG_W    <- 6.0          # inches
FIG_H    <- 5.4
FIG_DPI  <- 300
COL_UP   <- "#D62728"    # up: red
COL_DOWN <- "#2E86C1"    # down: blue
COL_NS   <- "grey78"     # not significant: grey
PT_ALPHA <- 0.45
PT_SIZE  <- 0.55


# --- Functions ---

classify <- function(dt, fc_thr = THR_FC, p_thr = THR_P, p_col = P_COL) {
  pv <- dt[[p_col]]
  dt[, sig := fifelse(!is.na(pv) & pv < p_thr & logFC >=  fc_thr, "up",
              fifelse(!is.na(pv) & pv < p_thr & logFC <= -fc_thr, "down", "ns"))]
  dt[]
}

count_sig <- function(dt) {
  c(total = nrow(dt),
    up    = dt[sig == "up",   .N],
    down  = dt[sig == "down", .N],
    ns    = dt[sig == "ns",   .N])
}

save_fig <- function(p, name, w = FIG_W, h = FIG_H) {
  dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(DIR_RESULTS, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(DIR_RESULTS, paste0(name, ".png")), p, width = w, height = h,
         dpi = FIG_DPI, bg = "white")
  cat(sprintf("    FIG %s.pdf / .png\n", name))
}

plot_volcano <- function(dt) {
  d <- copy(dt)
  d[, sig := factor(sig, levels = c("up", "ns", "down"))]
  d[, neglog10p := -log10(d[[P_COL]])]
  n <- count_sig(dt)

  p <- ggplot(d, aes(x = logFC, y = neglog10p)) +
    geom_point(aes(colour = sig), alpha = PT_ALPHA, size = PT_SIZE) +
    scale_colour_manual(
      values = c(up = COL_UP, ns = COL_NS, down = COL_DOWN),
      breaks = c("up", "ns", "down"),
      labels = c(sprintf("Sig_Up (%d)", n["up"]),
                 sprintf("NoDiff (%d)", n["ns"]),
                 sprintf("Sig_Down (%d)", n["down"])),
      name = NULL) +
    geom_vline(xintercept = c(-THR_FC, THR_FC), linetype = 2,
               colour = "grey35", linewidth = 0.3) +
    geom_hline(yintercept = -log10(THR_P), linetype = 2,
               colour = "grey35", linewidth = 0.3) +
    labs(x = expression(log[2]*"(FC)"),
         y = expression(-log[10]*"(pVal)")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7),
          panel.grid.minor = element_blank())

  if (!NO_LABEL && "gene_name" %in% names(d)) {
    lab <- d[!is.na(gene_name) & gene_name %in% LABEL_GENES & sig != "ns"]
    lab <- unique(lab, by = "gene_name")
    if (nrow(lab) > 0)
      p <- p + geom_text_repel(data = lab, aes(label = gene_name), size = 2.4,
                               max.overlaps = 30, segment.size = 0.25,
                               min.segment.length = 0, box.padding = 0.3,
                               show.legend = FALSE)
  }
  p
}


# --- Main ---

if (!file.exists(IN_CSV))
  stop("Input csv not found: ", IN_CSV, "\n  Run 01_normalize_deg.R first")
deg <- classify(fread(IN_CSV, encoding = "UTF-8"))
cnt <- count_sig(deg)
cat(sprintf("    loaded %d genes: Sig_Up %d | NoDiff %d | Sig_Down %d\n",
            cnt["total"], cnt["up"], cnt["ns"], cnt["down"]))

save_fig(plot_volcano(deg), OUT_PREFIX)
