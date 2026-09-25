suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})


# --- Parameters ---

DIR_MATERIAL <- "material"
DIR_PROCESS  <- "process"

IN_XLSX    <- file.path(DIR_MATERIAL, "All Comparisons (mRNA).xlsx")
SHEET_UP   <- "up_iDC-hMPV_vs_iDC-Ctrl"
SHEET_DOWN <- "down_iDC-hMPV_vs_iDC-Ctrl"
OUT_CSV    <- file.path(DIR_PROCESS, "deg_iDC_hMPV_vs_Ctrl.csv")

DEG_COLS <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR", "regulation",
              "hMPV1_raw", "hMPV2_raw", "Ctrl1_raw", "Ctrl2_raw",
              "hMPV1_norm", "hMPV2_norm", "Ctrl1_norm", "Ctrl2_norm",
              "gene_name", "gene_biotype", "GeneID", "Synonyms", "dbXrefs",
              "description", "Pathway", "GO_CC", "GO_MF", "GO_BP")

HEADER_KEY      <- "gene_id"         # value marking the real header row
GENE_ID_PATTERN <- "^ENSG|^MSTRG"    # valid gene ID prefixes


# --- Functions ---

coerce_types <- function(d) {
  for (j in seq_along(d)) {
    x <- d[[j]]
    if (is.character(x)) {
      y <- suppressWarnings(as.numeric(x))
      bad <- is.na(y) & !is.na(x) & nzchar(x)
      if (!any(bad)) d[[j]] <- y
    }
  }
  d
}

read_deg_sheet <- function(path, sheet, cols = DEG_COLS, key = HEADER_KEY) {
  raw <- as.data.table(suppressMessages(
    read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal")))
  hit <- which(as.character(raw[[1]]) == key)
  if (length(hit) == 0)
    stop(sprintf("[%s / %s] header row not found (first column should equal \"%s\")",
                 basename(path), sheet, key))
  d <- raw[(hit[1] + 1L):nrow(raw)]
  if (ncol(d) != length(cols))
    stop(sprintf("[%s / %s] column count mismatch: %d found, %d expected",
                 basename(path), sheet, ncol(d), length(cols)))
  setnames(d, cols)
  coerce_types(d)
}

save_csv <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, bom = TRUE, na = "")
  cat(sprintf("    OK  %s  (%d rows x %d columns)\n", basename(path), nrow(dt), ncol(dt)))
}


# --- Main ---

if (!file.exists(IN_XLSX)) stop("Input file not found: ", IN_XLSX)

avail <- excel_sheets(IN_XLSX)
parts <- list()
for (sn in c(SHEET_UP, SHEET_DOWN)) {
  if (!sn %in% avail) { cat("    ! missing sheet: ", sn, "\n"); next }
  d <- copy(read_deg_sheet(IN_XLSX, sn))     # copy() first to avoid shallow-copy warning
  d[, direction := ifelse(sn == SHEET_UP, "up", "down")]
  cat(sprintf("    %-30s %6d rows\n", sn, nrow(d)))
  parts[[length(parts) + 1L]] <- d
}
if (length(parts) == 0) stop("Neither sheet is available; check SHEET_UP / SHEET_DOWN")

deg <- rbindlist(parts, fill = TRUE)

stopifnot("invalid gene_id (header row likely mislocated)" =
            all(grepl(GENE_ID_PATTERN, deg$gene_id)))
cat(sprintf("    merged: %d genes (up %d + down %d)\n",
            nrow(deg), deg[direction == "up", .N], deg[direction == "down", .N]))

save_csv(deg, OUT_CSV)
