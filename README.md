# m6a-mrna-analysis

This repository contains two sets of R scripts that draw a differential expression volcano plot and an m6A–mRNA correlation scatter plot from sequencing result tables.

## Structure

```
deg-volcano/              differential expression volcano plot
m6a-mrna-correlation/     m6A methylation vs. mRNA expression correlation scatter plot
```

## Requirements

- R 4.5.1
- Packages: `readxl`, `data.table`, `ggplot2`, `ggrepel`

```r
install.packages(c("readxl", "data.table", "ggplot2", "ggrepel"))
```

Run the scripts in numeric order, for example:

```
Rscript 01_normalize_deg.R
Rscript 02_plot_volcano.R
```

## deg-volcano: differential expression volcano plot

| Script | Purpose |
|---|---|
| `01_normalize_deg.R` | Reads the up- and down-regulated sheets of the result workbook, locates the header row, unifies column names and data types, and merges the sheets into a single normalized csv |
| `02_plot_volcano.R` | Classifies genes as up-regulated, down-regulated or not significant, draws the volcano plot, and reports the class counts in the legend |

## m6a-mrna-correlation: m6A–mRNA correlation scatter plot

| Script | Purpose |
|---|---|
| `01_normalize_inputs.R` | Normalizes the two input workbooks (differential expression, differential methylation) |
| `02_pair_points.R` | Recomputes the x and y coordinates of each "gene × region" combination |
| `03_plot_scatter.R` | Draws scatter points and regression lines per region, and outputs the figures and the statistics table |
