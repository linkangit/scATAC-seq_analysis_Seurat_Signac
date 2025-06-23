# Load required libraries for single-cell ATAC-seq and RNA-seq analysis
library(Signac)                # For scATAC-seq analysis
library(Seurat)                # For integrated single-cell analysis
library(EnsDb.Hsapiens.v75)    # Human gene annotation database
library(tidyverse)             # Data manipulation and visualization
library(SingleR)               # Automated cell type annotation

# -------------------------------
# 1. Data Import and Setup
# -------------------------------

# Preview the fragment file (first 10 rows) to inspect format
frag_preview <- read.delim('atac_v1_pbmc_10k_fragments.tsv.gz', header = FALSE, nrows = 10)
head(frag_preview)

# Load scATAC-seq count matrix (10x Genomics format)
atac_counts <- Read10X_h5('atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
atac_counts[1:10, 1:10]  # Check a subset of the count matrix

# Construct a ChromatinAssay object, specifying peak format and QC thresholds
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),                              # Genomic coordinate separator
  fragments = "atac_v1_pbmc_10k_fragments.tsv.gz",# Path to fragment file
  min.cells = 10,                                 # Minimum cells per feature
  min.features = 200                              # Minimum features per cell
)

str(chrom_assay)  # Inspect the ChromatinAssay structure

# Import cell metadata (e.g., QC metrics, barcodes)
cell_metadata <- read.csv('atac_v1_pbmc_10k_singlecell.csv', header = TRUE, row.names = 1)
View(cell_metadata)

# Initialize Seurat object for ATAC data, using the ChromatinAssay and metadata
pbmc_atac <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = cell_metadata,
  assay = 'ATAC'
)

str(pbmc_atac)  # Confirm object structure

# -------------------------------
# 2. Gene Annotation
# -------------------------------

# Retrieve gene annotations from EnsDb and convert to UCSC-style chromosome names
gene_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(gene_annotations) <- paste0('chr', seqlevels(gene_annotations))

# Attach gene annotations to the ATAC assay for downstream analyses
Annotation(pbmc_atac) <- gene_annotations

# -------------------------------
# 3. Quality Control (QC)
# -------------------------------

# Compute nucleosome signal per cell (proxy for chromatin structure)
pbmc_atac <- NucleosomeSignal(pbmc_atac)

# Calculate TSS (Transcription Start Site) enrichment per cell (indicator of data quality)
pbmc_atac <- TSSEnrichment(object = pbmc_atac, fast = FALSE)

# Add additional QC metrics: blacklist ratio and fraction of reads in peaks
pbmc_atac$blacklist_ratio <- pbmc_atac$blacklist_region_fragments / pbmc_atac$peak_region_fragments
pbmc_atac$pct_reads_in_peaks <- pbmc_atac$peak_region_fragments / pbmc_atac$passed_filters * 100

View(pbmc_atac@meta.data)  # Inspect all computed QC metrics

# -------------------------------
# 4. QC Visualization
# -------------------------------

# Visualize relationships between key QC metrics
qc_plot1 <- DensityScatter(pbmc_atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
qc_plot2 <- DensityScatter(pbmc_atac, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
qc_plot1 | qc_plot2

# Violin plots for major QC features
VlnPlot(
  object = pbmc_atac,
  features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 6
)

# -------------------------------
# 5. Filtering Low-Quality Cells
# -------------------------------

# Subset cells based on stringent QC thresholds
pbmc_atac <- subset(
  x = pbmc_atac,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

# -------------------------------
# 6. Normalization & Dimensionality Reduction
# -------------------------------

# Normalize ATAC data using TF-IDF (term frequency-inverse document frequency)
pbmc_atac <- RunTFIDF(pbmc_atac)

# Identify top variable features (peaks)
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = 'q0')

# Perform SVD (Singular Value Decomposition) for linear dimensionality reduction
pbmc_atac <- RunSVD(pbmc_atac)

# Assess correlation between sequencing depth and principal components
DepthCor(pbmc_atac)

# -------------------------------
# 7. Clustering & Visualization
# -------------------------------

# Non-linear dimensionality reduction (UMAP) on LSI components 2-30
pbmc_atac <- RunUMAP(pbmc_atac, reduction = 'lsi', dims = 2:30)

# Construct nearest neighbor graph and cluster cells
pbmc_atac <- FindNeighbors(pbmc_atac, reduction = 'lsi', dims = 2:30)
pbmc_atac <- FindClusters(pbmc_atac, algorithm = 3)

# Visualize clusters in UMAP space
DimPlot(pbmc_atac, label = TRUE) + NoLegend()

# -------------------------------
# 8. Gene Activity Matrix
# -------------------------------

# Compute gene activity scores (proxy for gene expression from ATAC peaks)
gene_activity_matrix <- GeneActivity(pbmc_atac)
gene_activity_matrix[1:10, 1:10]

# Add gene activity as a new assay and normalize
pbmc_atac[['RNA']] <- CreateAssayObject(counts = gene_activity_matrix)
pbmc_atac <- NormalizeData(
  object = pbmc_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc_atac$nCount_RNA)
)

# Set RNA as default assay for marker gene visualization
DefaultAssay(pbmc_atac) <- 'RNA'

# Visualize canonical marker gene activity across clusters
FeaturePlot(
  pbmc_atac,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# -------------------------------
# 9. Load and Visualize scRNA-seq Reference
# -------------------------------

# Load pre-processed scRNA-seq data (should contain cell type annotations)
pbmc_rna <- readRDS('pbmc_10k_v3.rds')
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

View(pbmc_rna@meta.data)

# Compare UMAPs of ATAC and RNA datasets before integration
p_umap_atac <- DimPlot(pbmc_atac, reduction = 'umap') + NoLegend() + ggtitle('scATAC-Seq')
p_umap_rna  <- DimPlot(pbmc_rna, reduction = 'umap', group.by = 'celltype', repel = TRUE, label = TRUE) + ggtitle('scRNA-Seq') + NoLegend()
p_umap_atac | p_umap_rna

# -------------------------------
# 10. Integration: Label Transfer
# -------------------------------

# Identify anchors between scRNA-seq reference and scATAC-seq query
transfer_anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc_atac,
  reduction = 'pcaproject'  # Faster than CCA for large datasets
)

# Transfer cell type labels from RNA to ATAC data
predicted_labels <- TransferData(
  anchorset = transfer_anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc_atac[['lsi']],
  dims = 2:30
)

head(predicted_labels)

# Add predicted cell type labels to ATAC metadata
pbmc_atac <- AddMetaData(pbmc_atac, metadata = predicted_labels)
View(pbmc_atac@meta.data)

# Visualize transferred cell type labels on ATAC UMAP
plot_atac_labels <- DimPlot(
  pbmc_atac, reduction = 'umap', group.by = 'predicted.id',
  label = TRUE, repel = TRUE
) + NoLegend() + ggtitle('scATAC-Seq')

plot_rna_labels <- DimPlot(
  pbmc_rna, reduction = 'umap', group.by = 'celltype',
  label = TRUE, repel = TRUE
) + NoLegend() + ggtitle('scRNA-Seq')

plot_atac_labels | plot_rna_labels

# -------------------------------
# 11. Differential Accessibility Analysis
# -------------------------------

# Set cell identities to predicted cell types for comparison
Idents(pbmc_atac) <- pbmc_atac$predicted.id

# Switch back to ATAC assay for peak-level analysis
DefaultAssay(pbmc_atac) <- 'ATAC'

# Identify differentially accessible peaks between two cell types
da_peaks <- FindMarkers(
  object = pbmc_atac,
  ident.1 = 'CD4 Naive',
  ident.2 = 'CD14+ Monocytes',
  test.use = 'LR',                # Logistic regression
  latent.vars = 'nCount_ATAC'     # Adjust for sequencing depth
)

head(da_peaks)

# Visualize accessibility of top differential peak
vln_da <- VlnPlot(
  pbmc_atac, features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c('CD4 Naive', 'CD14+ Monocytes')
)

feat_da <- FeaturePlot(
  pbmc_atac, features = rownames(da_peaks)[1], pt.size = 0.1
)

vln_da | feat_da

# Calculate and sort log2 fold changes between groups
fold_changes <- FoldChange(
  object = pbmc_atac,
  ident.1 = 'CD4 Naive',
  ident.2 = 'CD14+ Monocytes'
)
fold_changes <- fold_changes[order(fold_changes$avg_log2FC, decreasing = TRUE),]
head(fold_changes)

# -------------------------------
# 12. Genomic Region Visualization
# -------------------------------

# Set plotting order for cell types
levels(pbmc_atac) <- unique(pbmc_atac$predicted.id)

# Visualize coverage at a top differential peak (with extended flanking regions)
CoveragePlot(
  object = pbmc_atac,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

# Optionally, launch an interactive coverage browser for a gene of interest
CoverageBrowser(pbmc_atac, region = 'CD8A')
