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

# Fetch mitochondrial genes from Ensembl
#   - Attributes: ENSEMBL gene IDs and chromosome names
#   - Filter: Mitochondrial chromosome (MT)
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mt_genes <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters = "chromosome_name",
  values = "MT",  # mitochondrial chromosome
  mart = ensembl
)
# Extract ENSEMBL IDs and validate against count matrix
mt_ensembl_ids <- mt_genes$ensembl_gene_id
valid_mt_genes <- intersect(mt_ensembl_ids,
                            rownames(counts)) # Retain actually expressed genes
# Add mitochondrial QC metrics to SCE object
#   - Calculate mitochondrial gene percentage per cell (Mito)
sce <- addPerCellQC(sce, subsets = list(Mito = valid_mt_genes))

# Create QC output directory if not exists
#   - dir: Parent directory path
#   - qc_path: Subdirectory for storing QC plots ("QC plots")
qc_path <- file.path(dir, "QC plots")
if (!dir.exists(qc_path)) {
  dir.create(qc_path, recursive = TRUE)
}

# Plot 1: Library Size vs Detected Genes
#   - x: Total UMI counts (log10-scaled)
#   - y: Number of detected genes per cell
#   - color: Batch information (technical replicates)
#   - Visualizes potential low-quality cells (bottom-left quadrant)
p1 <- plotColData(sce, 
                  x = "sum",                    # Total UMI counts per cell
                  y = "detected",               # Number of detected genes
                  colour_by = "BatchInfo") +
  labs(title = "Library Size vs Genes Detected",
       x = "Total UMI Counts (×10³)",
       y = "Number of Detected Genes",
       color = "Donor") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "none" # Hide legend for cleaner visualization
  ) +
  scale_x_continuous(labels = function(x) x/1000) + # Convert counts to thousands
  scale_color_viridis_d() # Color-blind friendly palette
# Save Plot 1
p1_path <- file.path(qc_path, "BatchEffect_UMI_vs_Genes.png")
ggsave(p1_path, p1, width = 10, height = 6, dpi = 600)

# Plot 2: Mitochondrial Gene Percentage
#   - x: Total UMI counts (log10-scaled)
#   - y: Percentage of mitochondrial reads (MT%)
#   - Red dashed line: Typical QC threshold (10%)
#   - High MT% indicates cell stress/apoptosis[6](@ref)
p2 <- plotColData(sce, 
                  x = "sum",                   # Total UMI counts
                  y = "subsets_Mito_percent",  # MT gene percentage
                  colour_by = "BatchInfo") +
  labs(title = "Mitochondrial Gene Percentage",
       x = "Total UMI Counts (×10³)",
       y = "MT Gene % of Total Counts",
       color = "Donor") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "none",
  ) +
  scale_x_continuous(labels = function(x) x/1000) +
  scale_color_viridis_d() +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") # QC threshold
# Save Plot 2
p2_path <- file.path(qc_path, "BatchEffect_UMI_vs_MTpercent.png")
ggsave(p2_path, p2, width = 10, height = 6, dpi = 600)

# Plot 3: Combined QC Metrics
#   - x: Total UMI counts (log10-scaled)
#   - y: Number of detected genes per cell
#   - color: Mitochondrial gene percentage (MT%)
#     - Blue: Low MT% (<5%) - healthy cells
#     - Grey: Intermediate MT% (5-10%) - borderline cells
#     - Red: High MT% (>10%) - potentially stressed/dying cells
#   - Visualizes the relationship between library size, gene detection, and mitochondrial content
#   - Helps identify low-quality cells in bottom-left quadrant (low UMI + low genes + high MT%)
p3 <- plotColData(sce, 
                  x = "sum",        # Total UMI counts per cell
                  y = "detected",   # Number of detected genes
                  colour_by = "subsets_Mito_percent") +  # Mitochondrial percentage
  labs(title = "Library Size vs Genes Detected",
       subtitle = "Colored by Mitochondrial Gene Percentage (MT%)",
       x = "Total UMI Counts (×10³)",
       y = "Number of Detected Genes",
       color = "MT Gene %"
  ) +
  scale_color_gradient2(
    low = "blue", 
    mid = "grey", 
    high = "red", 
    midpoint = 10
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold")
  )
# Save Plot 3
p3_path <- file.path(qc_path, "QC_Filtering_Matrix.png")
ggsave(p3_path, p3, width = 10, height = 6, dpi = 600)

# Quality Control Filtering (Based on Visualization Plots)
qc.lib <- sce$sum > 500                   # Keep cells with >500 total UMIs
qc.nexprs <- sce$detected > 1000          # Keep cells detecting >1000 genes
qc.mito <- sce$subsets_Mito_percent < 10  # MT% < 10%
keep_genes <- rowSums(counts(sce) > 0) > 20  # Expressed in at least 20 cells
# Filter SingleCellExperiment object
sce <- sce[keep_genes, qc.lib & qc.nexprs & qc.mito]

sce <- logNormCounts(sce)         # Data normalization
dec <- modelGeneVar(sce)          # Gene variance modeling - decomposes technical and biological variation
hvg <- getTopHVGs(dec, n = 2000)  # Selects top 2000 most biologically variable genes
# Dimensional reduction (PCA -> UMAP)
sce <- runPCA(sce,
              subset_row = hvg,   # Use only HVGs for PCA
              exprs_values = "logcounts",
              BSPARAM = BiocSingular::IrlbaParam())  # Fast approximate PCA
umap <- calculateUMAP(sce, 
                      dimred = "PCA",
                      n_neighbors = 15,
                      min_dist = 0.1,
                      metric = "cosine") # Suitable for sparse scRNA-seq data
reducedDim(sce, "UMAP") <- umap

# Batch effect visualization - UMAP colored by batch
p4 <- plotReducedDim(sce, dimred = "UMAP",
                     colour_by = "BatchInfo",
                     point_size = 1.5) +
  ggtitle("UMAP by Batch") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# Save Plot 4
p4_path <- file.path(qc_path, "QC_Filtering_Matrix.png")
ggsave(p4_path, p4, width = 10, height = 6, dpi = 600)

# Batch Effect Quantification
#   - Calculates mean Euclidean distance of cells to their batch centroid in UMAP space
#   - Higher values indicate more dispersed batches, suggesting stronger batch effects
batch_score <- reducedDim(sce, "UMAP") %>% 
  as.data.frame() %>%
  mutate(Batch = sce$BatchInfo) %>% 
  group_by(Batch) %>% 
  summarise(
    Distance_to_centroid = mean(stats::dist(cbind(UMAP1, UMAP2))),
    .groups = "drop"
  )
# Interpretation:
#   - Scores range from 3.41 (donor3_S1) to 9.65 (donor4_S1)
#   - Similar scores for technical replicates (e.g. donor2_S1/S2 = 7.39/7.39)
#   - Large differences suggest batch effects needing correction
print(batch_score)

# Run BASiCS MCMC
chain <- BASiCS_MCMC(
  sce, 
  N = 30000,  # Total iterations
  Thin = 15,  # Keep every 15th sample
  Burn = 15000,  # Initial samples discarded
  Regression = TRUE, # Account for technical noise via regression
  WithSpikes = FALSE, # No spike-in controls available
  PriorParam = BASiCS_PriorParam(
    sce,
    PriorMu = "EmpiricalBayes"), # Empirical Bayes priors for mean expression
  RunName = "human_adult_liver", # Unique identifier for output
  Threads = 22, # Parallel processing
  StoreChains = TRUE, # Save complete MCMC chains
  StoreDir = dir
)
# Load saved MCMC chain (RDS format)
chain <- readRDS(file.path(dir, "chain_human_adult_liver.Rds"))

# Create diagnostic plots directory if not exists
#   - dir: Parent directory path
#   - mcmc_cd_path: Subdirectory for storing diagnostic plots (MCMC convergence diagnostics)
mcmc_cd_path <- file.path(dir, "MCMC convergence diagnostics")
if (!dir.exists(mcmc_cd_path)) {
  dir.create(mcmc_cd_path, recursive = TRUE)
}

# Geweke Diagnostic Plot
#   - test whether the mean of the first 10% of chain differs from last 50%
#   - red dashed lines at ±3 indicate significant divergence (p<0.01)
p5 <- BASiCS_DiagPlot(
  chain,
  Param = "mu",                       # Check convergence for gene means
  Measure = "geweke") +               # Geweke's convergence diagnostic
  theme(legend.position = "bottom") +
  geom_hline(yintercept = c(-3, 3),
             linetype = "dashed",     # Threshold lines for significance
             color = "red")
# Save Plot 5
p5_path <- file.path(mcmc_cd_path, "Geweke.png")
ggsave(p5_path, plot = p5, width = 10, height = 5, dpi = 600)

# Effective Sample Size (ESS) Plot
#   - measure how many independent samples the chain contains
#   - higher ESS indicates better mixing
p6 <- BASiCS_DiagPlot(
  chain,
  Param = "mu",
  Measure = "ess") +                  # Effective sample size
  theme(legend.position = "bottom")
# Save Plot 6
p6_path <- file.path(mcmc_cd_path, "ESS.png")
ggsave(p6_path, plot = p6, width = 10, height = 5, dpi = 600)

# Combined Diagnostic Plot
#   - visualise both Geweke and ESS diagnostics side-by-side
p7 <- geweke + ess + plot_annotation(tag_levels = "A") # Add panel labels
p7_path <- file.path(mcmc_cd_path, "Geweke_ESS.png")
ggsave(p7_path, plot = p7, width = 10, height = 5, dpi = 600)