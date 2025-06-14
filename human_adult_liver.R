## Load R packages quietly
suppressPackageStartupMessages({
  library(Seurat)
  library(scran)
  library(scater)
  library(batchelor)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(ggpubr)
  library(biomaRt)
  library(SingleCellExperiment)
  library(BASiCS)
  library(arrow)
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
      dplyr::select(
        where(~ n_distinct(.) > 1), # Remove constant columns
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment.ENA_SAMPLE." # Comment.ENA_SAMPLE. as BatchInfo
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

sce_dup <- sce
sce_dup <- logNormCounts(sce_dup)         # Data normalization
dec <- modelGeneVar(sce_dup)          # Gene variance modeling - decomposes technical and biological variation
hvg <- getTopHVGs(dec, n = 2000)  # Selects top 2000 most biologically variable genes
# Dimensional reduction (PCA -> UMAP)
sce_dup <- runPCA(sce_dup,
              subset_row = hvg,   # Use only HVGs for PCA
              exprs_values = "logcounts",
              BSPARAM = BiocSingular::IrlbaParam())  # Fast approximate PCA
umap <- calculateUMAP(sce_dup, 
                      dimred = "PCA",
                      n_neighbors = 15,
                      min_dist = 0.1,
                      metric = "cosine") # Suitable for sparse scRNA-seq data
reducedDim(sce_dup, "UMAP") <- umap

# Batch effect visualization - UMAP colored by batch
p4 <- plotReducedDim(sce_dup, dimred = "UMAP",
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
batch_score <- reducedDim(sce_dup, "UMAP") %>% 
  as.data.frame() %>%
  mutate(Batch = sce_dup$BatchInfo) %>% 
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
p7 <- p5 + p6 + plot_annotation(tag_levels = "A") # Add panel labels
p7_path <- file.path(mcmc_cd_path, "Geweke_ESS.png")
ggsave(p7_path, plot = p7, width = 10, height = 5, dpi = 600)

# Extract MCMC summary statistics from the BASiCS chain object
mcmc_summary <- Summary(chain)
# Create data frame of key gene-level parameters:
#   - mean: Baseline expression (μ)
#   - overdisp: Biological overdispersion (δ)
#   - residualoverdisp: Residual overdispersion (ε)
mcmc_parameter <- data.frame(
  mean = displaySummaryBASiCS(mcmc_summary, Param = "mu")[, 1],
  overdisp = displaySummaryBASiCS(mcmc_summary, Param = "delta")[, 1],
  residualoverdisp = displaySummaryBASiCS(mcmc_summary, Param = "epsilon")[, 1]
)

# Detect Highly Variable Genes (HVGs)
hvg <- BASiCS_DetectHVG(
  chain,
  EpsilonThreshold = log(2)) # 2-fold change in residual variation
# Detect Lowly Variable Genes (LVGs)
lvg <- BASiCS_DetectLVG(
  chain,
  EpsilonThreshold = -log(2))
# Merge HVG and LVG results into unified table
vg_table <- merge(
  as.data.frame(lvg, Filter = FALSE), # LVG results (include all genes)
  as.data.frame(hvg, Filter = FALSE), # HVG results (include all genes)
  by = c("GeneName", "GeneIndex", "Mu", "Delta", "Epsilon"), # Merge keys
  suffixes = c("LVG", "HVG")          # Suffixes for overlapping columns
)
# Genes are exclusive as HVG/LVG/Neither
vg_table$VG <- "Not HVG or LVG"     # Default classification
vg_table$VG[vg_table$HVG] <- "HVG"  # Mark HVGs (ε > log(2))
vg_table$VG[vg_table$LVG] <- "LVG"  # Mark LVGs (ε < -log(2))
# Integrate results with gene metadata in SingleCellExperiment object
rowData(sce) <- merge(
  data.frame(ensembl_gene_id = rownames(sce), rowData(sce)),
  vg_table,
  by.x = "ensembl_gene_id", by.y = "GeneName",
  sort = FALSE
)

dir_targetscan <- "TargetScan"
# Read TargetScan's predicted target gene information file
targets_info <- read.delim(
  file.path(dir_targetscan, "Predicted_Targets_Info.default_predictions.txt"), 
  stringsAsFactors = FALSE
)
# Read TargetScan's context score predictions 
conserved_scores <- read.delim(
  file.path(dir_targetscan, "Predicted_Targets_Context_Scores.default_predictions.txt"),
  stringsAsFactors = FALSE
)

# Extraction of liver samples and expression matrix (miRNA TissueAtlas 2025)
dir_atlas <- "miRNATissueAtlas"
# Read sample large expression matrices
obs <- read_parquet(file.path(dir_atlas,"expression.parquet"))
# Read sample metadata containing tissue annotations
meta <- read_parquet(file.path(dir_atlas,"metadata.parquet"))
# Merge expression data with metadata using unique sample IDs
human_mirna <-  merge(obs, meta, by = "__index_level_0__")
# Filter for liver samples and select miRNA expression columns
liver_mirna <- human_mirna %>%
  filter(Tissue == "liver") %>%
  dplyr::select(2:(ncol(.)-6))
# Calculation of average miRNA expression in liver
#   - compute mean expression across all liver samples for each miRNA
miRNA_means <- liver_mirna %>% 
  colMeans(na.rm = TRUE) %>% 
  enframe(name = "miRNA", value = "mean") %>%
  arrange(desc(mean))
# Extraction of top 10% miRNA names
#   - select all miRNA names
all_miRNAs <- miRNA_means$miRNA

# Acquisition of validated liver miRNAs from miRTarBase
dir_mirtarbase <- "miRTarBase"
mirtarbase <- read.csv(file.path(dir_mirtarbase, "hsa_MTI.csv"), stringsAsFactors = FALSE) %>% 
  filter(miRNA %in% all_miRNAs) %>% distinct(miRTarBase.ID, .keep_all = TRUE)
validated_liver_mirna <- unique(mirtarbase$miRNA) %>% sub("^hsa-", "", .)

# Extraction of all target genes for liver miRNAs from TargetScan
liver_targets <- targets_info %>%
  filter(Species.ID == 9606) %>%  # Human
  mutate(miR.Family = strsplit(miR.Family, "/")) %>% 
  unnest(miR.Family) %>% 
  # Screen for matching miRNA families
  filter(miR.Family %in% validated_liver_mirna) %>% 
  mutate(Gene.ID = sub("\\..*", "", Gene.ID)) %>% 
  # Join with conserved scores data
  left_join(
    conserved_scores %>%
      filter(Gene.Tax.ID == 9606) %>% # Human
      mutate(
        Gene.ID = sub("\\..*", "", Gene.ID),
        miRNA = sub("^.*-miR-", "miR-", miRNA) # Remove prefix
      ),
    by = c("Gene.ID", "miR.Family" = "miRNA"),
    relationship = "many-to-many"
  ) %>%
  # Filter based on TargetScan scores
  filter(
    context...score < 0,
    context...score.percentile >90
  )

# Generate gene-miRNA interaction summary table
#   - Input: liver_targets (dataframe containing Gene.ID and miR.Family columns)
#   - Output: tibble with aggregated miRNA targeting information per gene
gene_mirna_summary <- liver_targets %>%
  # Select only relevant columns to minimise memory usage
  dplyr::select(Gene.ID, miR.Family) %>%
  # Calculate the number of target sites for each gene
  group_by(Gene.ID) %>%
  summarise(
    target.number = n(),
    miRNA = paste(unique(miR.Family), collapse = "/")  # Merge all target miRNAs
  ) %>%
  ungroup() # Remove grouping structure for downstream operations

# Label target genes in single-cell data
rowData(sce)$miRNA_target <- ifelse(
  rownames(sce) %in% gene_mirna_summary$Gene.ID,
  "Target",
  "Non-target"
)

# Analyse target gene distribution
gene_stats <- as.data.frame(rowData(sce)) %>%
  mutate(
    gene_class = factor(VG, levels = c("HVG", "LVG", "Not HVG or LVG")),
    is_target = miRNA_target == "Target"
  )

# Create nosie comparison plots directory if not exists
#   - dir: Parent directory path
#   - noise_path: Subdirectory for storing nosie comparison plots (Gene expression nosie comparison)
noise_path <- file.path(dir, "Gene expression nosie comparison")
if (!dir.exists(noise_path)) {
  dir.create(noise_path, recursive = TRUE)
}
# Chi-square test for target gene proportions
count_data <- as.data.frame(rowData(sce)) %>% 
  filter(VG %in% c("HVG", "LVG")) %>% 
  group_by(VG, miRNA_target) %>% 
  summarise(n = n()) %>% 
  tidyr::pivot_wider(names_from = miRNA_target, values_from = n, values_fill = 0)
chisq_test <- chisq.test(as.matrix(count_data[, c("Target", "Non-target")]))
p_label_chisq <- ifelse(chisq_test$p.value < 0.001,
                        "P < 0.001",
                        sprintf("P = %.3f", chisq_test$p.value))
x2_label <- sprintf("X² = %.2f", chisq_test$statistic)

# Target gene proportion plot
p8 <- ggplot(gene_stats %>% filter(gene_class %in% c("HVG", "LVG")), 
             aes(x = gene_class, fill = is_target)) +
  geom_bar(position = "fill", width = 0.6, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("skyblue", "salmon"),
                    labels = c("Non-target", "Target"),
                    name = "") + 
  labs(x = "Gene Category", y = "Proportion", 
       title = "Liver miRNA Targets Distribution",
       subtitle = paste("Chi-square test:", x2_label, ",", p_label_chisq)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "top")
# Save Plot 8
p8_path <- file.path(noise_path, "chi.png")
ggsave(p8_path, plot = p8, width = 10, height = 5, dpi = 600)

# Wilcoxon test for expression noise
noise_comparison <- wilcox.test(
  Delta ~ miRNA_target, 
  data = as.data.frame(rowData(sce))
)
p_label <- ifelse(noise_comparison$p.value < 0.001, 
                  "P < 0.001", 
                  sprintf("P = %.3f", noise_comparison$p.value))
w_label <- sprintf("W = %.0f", noise_comparison$statistic)

# Expression noise comparison boxplot
p9 <- ggplot(rowData(sce), aes(x = miRNA_target, y = Epsilon, fill = miRNA_target)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 21,
    outlier.size = 2,
    outlier.color = "black",
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  labs(x = NULL, y = "Residual Over-dispersion",
       title = "Expression Noise Comparison",
       subtitle = paste("Wilcoxon test:", w_label, ",", p_label)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "none")
# Save Plot 9
p9_path <- file.path(noise_path, "wilcoxon.png")
ggsave(p9_path, plot = p9, width = 10, height = 5, dpi = 600)

# Combined plot
p10 <- p8 + p9 +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')
# Save Plot 10
p10_path <- file.path(noise_path, "wilcoxon_chi.png")
ggsave(p10_path, plot = p10, width = 10, height = 5, dpi = 600)