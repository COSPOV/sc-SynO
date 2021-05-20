# Single nuclei transcriptomics of the Allen Brain atlas data set of mus musculus (https://celltypes.brain-map.org/)
# The underlying analysis is part of the manuscript entitled 
# "Automated annotation of rare-cell types from single-cell RNA-sequencing data through synthetic oversampling"
# Data anaylsis and visalizations were mainly generated with the Seurat R package (https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html)
# Analysis for the iRhythmics project (https://irhythmics.med.uni-rostock.de/)
# Creator: Markus Wolfien (markus.wolfien@gmail.com) or (markus.wolfien@uni-rostock.de)

# Shortcut do save current graphs within the display device as pdf 
dopdf <- function(filename = "dopdf.pdf", pdf.width = 16, pdf.height = 20) {
  dev.copy2pdf(file=filename, width = pdf.width, height = pdf.height)}

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2) # Needed for plotting
library(cowplot)
library(harmony)
library(data.table)
library(DropletUtils)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(magrittr)
library("biomaRt") # Used for transcript annotation
library("psych")   # needed for the dfOrder function
library(data.table)

# Download data from Allen Brain atlas with e.g., wget or curl
# Mice https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv

raw_counts<-read.table(file=paste0("matrix.csv"),sep=",", header = TRUE, row.names = 1, as.is = TRUE)
# If the single matrix file is too large for an R integration use $split -d -l 300000 matrix.csv brain_part_
# to split up the input as it was done in our case to use only the first 300.000 cells

raw_counts<-read.table(file=paste0("brain_part_00.csv"),sep=",", header = TRUE, row.names = 1, as.is = TRUE) 

# Seurat needs cells in colums and features in rows
# transposing of data needed
t_raw_counts <- transpose(raw_counts)
# get row and colnames in order
  colnames(t_raw_counts) <- rownames(raw_counts)
  rownames(t_raw_counts) <- colnames(raw_counts)
 col_raw <- colnames(raw_counts) # save the columns for the other data batches

head(raw_counts)

# Amount of raw counts is too large for a single seurat object, splitting of input data needed ...

########## Murine dataset #####################
brain_1 <- t_raw_counts[,1:29999]
brain_2 <- t_raw_counts[,30000:59999]
brain_3 <- t_raw_counts[,60000:89999]
brain_4 <- t_raw_counts[,90000:119999]
brain_5 <- t_raw_counts[,120000:149999]
brain_6 <- t_raw_counts[,150000:179999]
brain_7 <- t_raw_counts[,180000:209999]
brain_8 <- t_raw_counts[,210000:239999]
brain_9 <- t_raw_counts[,240000:269999]
brain_10 <- t_raw_counts[,270000:299999]

# Use individual Seurat objects and check for a consistent amount of included cells
m.brain_data_1 <- CreateSeuratObject(counts = brain_1, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_1
m.brain_data_2 <- CreateSeuratObject(counts = brain_2, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_2
m.brain_data_3 <- CreateSeuratObject(counts = brain_3, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_3
m.brain_data_4 <- CreateSeuratObject(counts = brain_4, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_4
m.brain_data_5 <- CreateSeuratObject(counts = brain_5, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_5
m.brain_data_6 <- CreateSeuratObject(counts = brain_6, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_6
m.brain_data_7 <- CreateSeuratObject(counts = brain_7, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_7
m.brain_data_8 <- CreateSeuratObject(counts = brain_8, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_8
m.brain_data_9 <- CreateSeuratObject(counts = brain_9, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_9
m.brain_data_10 <- CreateSeuratObject(counts = brain_10, min.cells = 3, min.features = 200, project = "brain_scRNAseq")
m.brain_data_10

# Detect and integrate mitochondrial genes (mt-DNA)
mito.genes_1 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_1), value = TRUE)
mito.genes_2 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_2), value = TRUE)
mito.genes_3 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_3), value = TRUE)
mito.genes_4 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_4), value = TRUE)
mito.genes_5 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_5), value = TRUE)
mito.genes_6 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_6), value = TRUE)
mito.genes_7 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_7), value = TRUE)
mito.genes_8 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_8), value = TRUE)
mito.genes_9 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_9), value = TRUE)
mito.genes_10 <- grep(pattern = "^mt-", x = rownames(x = m.brain_data_10), value = TRUE)


m.brain_data_1[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_1, pattern = "^mt-")
m.brain_data_2[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_2, pattern = "^mt-")
m.brain_data_3[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_3, pattern = "^mt-")
m.brain_data_4[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_4, pattern = "^mt-")
m.brain_data_5[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_5, pattern = "^mt-")
m.brain_data_6[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_6, pattern = "^mt-")
m.brain_data_7[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_7, pattern = "^mt-")
m.brain_data_8[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_8, pattern = "^mt-")
m.brain_data_9[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_9, pattern = "^mt-")
m.brain_data_10[["percent.mt"]] <- PercentageFeatureSet(m.brain_data_10, pattern = "^mt-")

# Filter low quality reads and potential cell doublets
m.brain_data_1 <- subset(m.brain_data_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_2 <- subset(m.brain_data_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_3 <- subset(m.brain_data_3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_4 <- subset(m.brain_data_4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_5 <- subset(m.brain_data_5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_6 <- subset(m.brain_data_6, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_7 <- subset(m.brain_data_7, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_8 <- subset(m.brain_data_8, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_9 <- subset(m.brain_data_9, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
m.brain_data_10 <- subset(m.brain_data_10, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize and transform data
m.brain_data_1 <- NormalizeData(m.brain_data_1) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_2 <- NormalizeData(m.brain_data_2) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_3 <- NormalizeData(m.brain_data_3) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_4 <- NormalizeData(m.brain_data_4) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_5 <- NormalizeData(m.brain_data_5) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_6 <- NormalizeData(m.brain_data_6) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_7 <- NormalizeData(m.brain_data_7) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_8 <- NormalizeData(m.brain_data_8) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_9 <- NormalizeData(m.brain_data_9) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)
m.brain_data_10 <- NormalizeData(m.brain_data_10) %>% FindVariableFeatures() %>% SCTransform() %>% RunPCA(verbose = FALSE)

# Prepare datasets for integration 
features <- SelectIntegrationFeatures(object.list = c(m.brain_data_1, m.brain_data_2, m.brain_data_3, m.brain_data_4, m.brain_data_5, m.brain_data_6, m.brain_data_7, m.brain_data_8,m.brain_data_9, m.brain_data_10), nfeatures = 3000)
brain.list <- PrepSCTIntegration(object.list = c(m.brain_data_1, m.brain_data_2, m.brain_data_3, m.brain_data_4, m.brain_data_5, m.brain_data_6, m.brain_data_7, m.brain_data_8,m.brain_data_9, m.brain_data_10), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = brain.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30)

# Integrate SCT normalized data via anchors
brain.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# Don't run ScaleData here after integration
brain.integrated <- RunPCA(brain.integrated, verbose = FALSE)
brain.integrated <- RunUMAP(brain.integrated, reduction = "pca", dims = 1:30)

pdf('m.brain_data_umap.pdf')
DimPlot(brain.integrated, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

# Cluster the cells
brain.integrated <- FindNeighbors(object = brain.integrated, dims = 1:30)
brain.integrated <- FindClusters(object = brain.integrated, resolution = 0.4)
pdf('m.brain_clusters.pdf')
DimPlot(brain.integrated, reduction = "umap", raster = FALSE, label = TRUE)
dev.off()

# Cluster 30 was picked as a cluster of interest - find markers for this cluster - use different test for compuation and use the joint list of markers as scSynO identification input
# In general more tests are recommended but are not a requirement
cluster30_markers_roc <- FindMarkers(brain.integrated, ident.1 = "30", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster30_markers_lr <- FindMarkers(brain.integrated, ident.1 = "30", logfc.threshold = 0.25, test.use = "LR", only.pos = TRUE)
cluster30_markers_t <- FindMarkers(brain.integrated, ident.1 = "30", logfc.threshold = 0.25, test.use = "t", only.pos = TRUE)

# Export markers as .csv file
write.csv(cluster30_markers_roc, file="cluster30_markers_roc.csv")
write.csv(cluster30_markers_lr, file="cluster30_markers_lr.csv") 
write.csv(cluster30_markers_t, file="cluster30_markers_t.csv")  

# Export the specific cluster 30 and the rest of the cells of the other clusters as a training input for scSynO
# Either use all transcripts as used here or directly use a reduced expression set of the identified markers only
cell_expression_30 <- as.matrix(GetAssayData(brain.integrated, slot = "data")[, WhichCells(brain.integrated, ident = "30")])
cell_expression <- as.matrix(GetAssayData(brain.integrated, slot = "data")[, WhichCells(brain.integrated, ident = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33"))])

# Export the data as .csv files
setwd("/media/md0/mw581/Projects/David/iRhythmics/Sc_heart_nuc_190910/Loras_expression")
write.csv(cell_expression, file="cell_expression_all_brain_300k.csv")
write.csv(cell_expression_30, file="cell_expression_30_brain_300k.csv")

sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS:   /usr/local/lib64/R/lib/libRblas.so
# LAPACK: /usr/local/lib64/R/lib/libRlapack.so

# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] magrittr_2.0.1     data.table_1.14.0  harmony_1.0        Rcpp_1.0.6
# [5] cowplot_1.1.1      ggplot2_3.3.3      Matrix_1.3-2       dplyr_1.0.5
# [9] SeuratObject_4.0.0 Seurat_4.0.1

# loaded via a namespace (and not attached):
# [1] nlme_3.1-152          matrixStats_0.58.0    spatstat.sparse_2.0-0
# [4] RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2
# [7] sctransform_0.3.2     tools_4.0.5           utf8_1.2.1
#[10] R6_2.5.0              irlba_2.3.3           rpart_4.1-15
#[13] KernSmooth_2.23-18    uwot_0.1.10           mgcv_1.8-34
#[16] lazyeval_0.2.2        colorspace_2.0-0      withr_2.4.2
#[19] tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.5
#[22] plotly_4.9.3          scales_1.1.1          lmtest_0.9-38
#[25] spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.4-3
#[28] goftest_1.2-2         stringr_1.4.0         digest_0.6.27
#[31] spatstat.utils_2.1-0  pkgconfig_2.0.3       htmltools_0.5.1.1
#[34] parallelly_1.24.0     fastmap_1.1.0         htmlwidgets_1.5.3
#[37] rlang_0.4.10          shiny_1.6.0           generics_0.1.0
#[40] zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2
#[43] patchwork_1.1.1       munsell_0.5.0         fansi_0.4.2
#[46] abind_1.4-5           reticulate_1.19       lifecycle_1.0.0
#[49] stringi_1.5.3         MASS_7.3-53.1         Rtsne_0.15
#[52] plyr_1.8.6            grid_4.0.5            parallel_4.0.5
#[55] listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1
#[58] crayon_1.4.1          miniUI_0.1.1.1        deldir_0.2-10
#[61] lattice_0.20-41       splines_4.0.5         tensor_1.5
#[64] pillar_1.6.0          igraph_1.2.6          spatstat.geom_2.1-0
#[67] future.apply_1.7.0    reshape2_1.4.4        codetools_0.2-18
#[70] leiden_0.3.7          glue_1.4.2            vctrs_0.3.7
#[73] png_0.1-7             httpuv_1.5.5          gtable_0.3.0
#[76] RANN_2.6.1            purrr_0.3.4           spatstat.core_2.1-2
#[79] polyclip_1.10-0       tidyr_1.1.3           scattermore_0.7
#[82] future_1.21.0         mime_0.10             xtable_1.8-4
#[85] later_1.1.0.1         survival_3.2-10       viridisLite_0.4.0
#[88] tibble_3.1.1          cluster_2.1.1         globals_0.14.0
#[91] fitdistrplus_1.1-3    ellipsis_0.3.1        ROCR_1.0-11
