suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# --- Parameters ---

DIR_PROCESS <- "process"
DIR_RESULTS <- "results"

IN_POINTS  <- file.path(DIR_PROCESS, "fig1E_points_hMPV.csv")
OUT_PREFIX <- "fig1E_scatter_hMPV"
OUT_STATS  <- file.path(DIR_RESULTS, "fig1E_stats_hMPV.csv")

REGION_LEVELS <- c("5'UTR", "StartC", "CDS", "StopC", "3'UTR")
REGION_COLS   <- c("5'UTR" = "#4D9221", "StartC" = "#F1A340", "CDS" = "#4393C3",
                   "StopC" = "#B2182B", "3'UTR" = "#762A83")

FIG_W    <- 7.5          # inches
FIG_H    <- 5.4
FIG_DPI  <- 300
PT_ALPHA <- 0.35
PT_SIZE  <- 0.7
LINE_W   <- 0.7


# --- Functions ---

fit_regions <- function(dt) {
  rbindlist(lapply(REGION_LEVELS, function(rg) {
    z <- dt[region == rg]
    if (nrow(z) < 3) return(data.table(region = rg, n = nrow(z)))
    ct <- suppressWarnings(cor.test(z$x, z$y))
    ft <- lm(y ~ x, data = z)
    data.table(region = rg, n = nrow(z),
               r = unname(ct$estimate), p = ct$p.value,
               slope = unname(coef(ft)[2]), intercept = unname(coef(ft)[1]))
  }), fill = TRUE)
}

region_label <- function(fit) {
  z <- fit[match(REGION_LEVELS, region)]
  setNames(sprintf("%s   r = %.3f, p %s", z$region, z$r,
                   ifelse(z$p < 1e-3, "< 0.001", sprintf("= %.3f", z$p))),
           z$region)
}

save_fig <- function(p, name, w = FIG_W, h = FIG_H) {
  dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(DIR_RESULTS, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(DIR_RESULTS, paste0(name, ".png")), p, width = w, height = h,
         dpi = FIG_DPI, bg = "white")
  cat(sprintf("    FIG %s.pdf / .png\n", name))
}

save_csv <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, bom = TRUE, na = "")
  cat(sprintf("    OK  %s  (%d rows x %d columns)\n", basename(path), nrow(dt), ncol(dt)))
}

plot_scatter <- function(dt, fit) {
  d <- copy(dt)
  d[, region := factor(region, levels = REGION_LEVELS)]

  ggplot(d, aes(x = x, y = y, colour = region)) +
    geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.25) +
    geom_point(alpha = PT_ALPHA, size = PT_SIZE) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = LINE_W) +
    scale_colour_manual(values = REGION_COLS, breaks = REGION_LEVELS,
                        labels = region_label(fit), name = NULL) +
    labs(x = expression("mRNA differential expression [" * log[2] * "(FC)]"),
         y = expression("m6A differential methylation [" * log[2] * "(FC)]")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7.5),
          panel.grid.minor = element_blank())
}


# --- Main ---

if (!file.exists(IN_POINTS))
  stop("Input file not found: ", IN_POINTS, "\n  Run 02_pair_points.R first")
pts <- fread(IN_POINTS, encoding = "UTF-8")
cat(sprintf("    loaded %d points in %d regions\n", nrow(pts), uniqueN(pts$region)))

fit <- fit_regions(pts)
print(fit[, .(region, n, r = round(r, 4), p = signif(p, 3),
              slope = round(slope, 4), intercept = round(intercept, 4))],
      row.names = FALSE)

save_fig(plot_scatter(pts, fit), OUT_PREFIX)
save_csv(fit[, .(region, n, r, p, slope, intercept)], OUT_STATS)
