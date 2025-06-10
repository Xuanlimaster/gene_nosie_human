## Load R packages quietly
suppressPackageStartupMessages({
  library(Seurat)
  library(scran)
  library(scater)
  library(Matrix)
  library(batchelor)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(ggpubr)
  library(biomaRt)
  library(SingleCellExperiment)
  library(BASiCS)
  library(harmony)
  library(arrow)
  library(tidyr)
})

dir <- "human_adult_liver"
# Load human scRNA-seq data from MatrixMarket format (.mtx)
#   - 'cells' file contains cell barcodes (colnames)
#   - 'features' file contains gene identifiers (rownames)
counts <- ReadMtx(
  mtx = file.path(dir, "E-MTAB-10553.aggregated_counts.mtx.gz"),
  cells = file.path(dir, "E-MTAB-10553.aggregated_counts.mtx_cols.gz"),
  features = file.path(dir, "E-MTAB-10553.aggregated_counts.mtx_rows.gz")
)

# Process metadata from TSV and SDRF files:
#   - Extract first 9 columns of cell metadata
#   - Add BioSD_SAMPLE ID by parsing 'id' column
#   - Merge with SDRF file to annotate batch/source information
metadata <- fread(file.path(dir, "E-MTAB-10553.cell_metadata.tsv")) %>%
  .[, 1:9] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path(dir, "E-MTAB-10553.sdrf.txt")) %>%
      dplyr::select(where(~ n_distinct(.) > 1)) %>% # Remove constant columns
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        ENA_SAMPLE = "Comment.ENA_SAMPLE.",
        Source = "Source.Name",
        BatchInfo = "Source.Name" # BatchInfo duplicates Source.Name
      ) %>%
      distinct(BioSD_SAMPLE, .keep_all = TRUE), 
    by = "BioSD_SAMPLE"
  ) %>%
  column_to_rownames("id") # Set cell IDs as rownames for SingleCellExperiment

# Create SingleCellExperiment object:
#   - Ensure counts and metadata are matched by cell IDs
#   - Store raw counts in 'assays' slot
#   - Store cell metadata in 'colData' slot
counts <- counts[, match(rownames(metadata), colnames(counts))] # Align cells
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)