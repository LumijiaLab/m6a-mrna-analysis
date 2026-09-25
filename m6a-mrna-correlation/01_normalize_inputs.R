suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})


# --- Parameters ---

DIR_MATERIAL <- "material"
DIR_PROCESS  <- "process"

IN_DEG     <- file.path(DIR_MATERIAL, "All Comparisons (mRNA).xlsx")
SHEETS_DEG <- c("up_iDC-hMPV_vs_iDC-Ctrl", "down_iDC-hMPV_vs_iDC-Ctrl")
OUT_DEG    <- file.path(DIR_PROCESS, "deg_iDC_hMPV_vs_Ctrl.csv")

IN_M6A     <- file.path(DIR_MATERIAL, "Differentially_methylated_sites.mRNA.xlsx")
SHEETS_M6A <- c("up.iDC-hMPV_vs_iDC-Ctrl", "down.iDC-hMPV_vs_iDC-Ctrl")
OUT_M6A    <- file.path(DIR_PROCESS, "m6a_diff_iDC_hMPV_vs_Ctrl.csv")

DEG_COLS <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR", "regulation",
              "hMPV1_raw", "hMPV2_raw", "Ctrl1_raw", "Ctrl2_raw",
              "hMPV1_norm", "hMPV2_norm", "Ctrl1_norm", "Ctrl2_norm",
              "gene_name", "gene_biotype", "GeneID", "Synonyms", "dbXrefs",
              "description", "Pathway", "GO_CC", "GO_MF", "GO_BP")

M6A_COLS <- c("chrom", "txStart", "txEnd", "PeakID", "score", "Peak_length",
              "transcript_id", "GeneName", "Foldchange", "P_value", "FDR",
              "Regulation", "GeneID", "Synonyms", "dbXrefs", "description",
              "Pathway", "GO_CC", "GO_MF", "GO_BP")

HEADER_KEY_DEG <- "gene_id"       # value marking the real header row
HEADER_KEY_M6A <- "chrom"

DEG_ID_PATTERN <- "^ENSG|^MSTRG"  # valid first-column prefixes
M6A_ID_PATTERN <- "^chr"


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

read_sheet <- function(path, sheet, cols, key) {
  raw <- as.data.table(suppressMessages(
    read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal")))
  hit <- which(as.character(raw[[1]]) == key)   # header sits below a file-specific description block
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

merge_sheets <- function(path, sheets, cols, key) {
  avail <- excel_sheets(path)
  parts <- list()
  for (sn in sheets) {
    if (!sn %in% avail) { cat("    ! missing sheet: ", sn, "\n"); next }
    d <- copy(read_sheet(path, sn, cols, key))   # copy() first to avoid shallow-copy warning
    d[, direction := ifelse(grepl("^up", sn), "up", "down")]
    cat(sprintf("    %-32s %6d rows\n", sn, nrow(d)))
    parts[[length(parts) + 1L]] <- d
  }
  if (length(parts) == 0) stop("None of the listed sheets is available: ", basename(path))
  out <- rbindlist(parts, fill = TRUE)
  attr(out, "n_up")   <- sum(out$direction == "up")
  attr(out, "n_down") <- sum(out$direction == "down")
  out
}

save_csv <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, bom = TRUE, na = "")
  cat(sprintf("    OK  %s  (%d rows x %d columns)\n", basename(path), nrow(dt), ncol(dt)))
}


# --- Main ---

if (!file.exists(IN_DEG)) stop("Input file not found: ", IN_DEG)
deg <- merge_sheets(IN_DEG, SHEETS_DEG, DEG_COLS, HEADER_KEY_DEG)
stopifnot("invalid gene_id (header row likely mislocated)" =
            all(grepl(DEG_ID_PATTERN, deg$gene_id)))
cat(sprintf("    merged: %d genes (up %d + down %d)\n",
            nrow(deg), attr(deg, "n_up"), attr(deg, "n_down")))
save_csv(deg, OUT_DEG)

if (!file.exists(IN_M6A)) stop("Input file not found: ", IN_M6A)
m6a <- merge_sheets(IN_M6A, SHEETS_M6A, M6A_COLS, HEADER_KEY_M6A)
stopifnot("invalid chrom (header row likely mislocated)" =
            all(grepl(M6A_ID_PATTERN, m6a$chrom)))
cat(sprintf("    merged: %d peaks (up %d + down %d)\n",
            nrow(m6a), attr(m6a, "n_up"), attr(m6a, "n_down")))
save_csv(m6a, OUT_M6A)
