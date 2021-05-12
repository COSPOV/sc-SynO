# Single-cell RNA-seq analysis of single-cell and nuclei data of adult murine whole hearts
# The underlying analysis is part of the manuscript entitled 
# "Automated annotation of rare-cell types from single-cell RNA-sequencing data through synthetic oversampling"
# Data anaylsis and visalizations were mainly generated with the Seurat R package (https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html)
# Analysis for the iRhythmics project (https://irhythmics.med.uni-rostock.de/)
# Creator: Markus Wolfien (markus.wolfien@gmail.com) or (markus.wolfien@uni-rostock.de)

# Define function to save current graphs within the display device as pdf (just for your convenience)
dopdf <- function(filename = "dopdf.pdf", pdf.width = 16, pdf.height = 20) {
  dev.copy2pdf(file=filename, width = pdf.width, height = pdf.height)}

# Please prepare the .fastq datasets with the provided UNIX script to distinguish between spliced and unspliced transcripts
# Preprocessing protocol from https://github.com/BUStools/getting_started/blob/master/kallisto_bus_mouse_nuclei_tutorial.ipynb

# Install/Load required packages
library(dplyr)
library(Seurat) # Single cell analysis package
library(ggplot2) # Needed for plotting
library(devtools) # To be able to install tools from github
install_github("immunogenomics/harmony")
library(harmony) # https://github.com/immunogenomics/harmony
install_github("velocyto-team/velocyto.R")
library(velocyto.R) # RNA velocity package
devtools::install_github("BUStools/BUSpaRse")
library(BUSpaRse) # Bustools extraction package
library(Rcpp)
library(data.table)
library(DropletUtils)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(magrittr)
library("biomaRt") # Used for transcript annotation
library("psych")   # needed for the dfOrder function

# Define your workspace location of the processed data

setwd(".../bus_output_FztDU") 	# Single-nuclei raw data of the Fzt:DU mice strain processed with the kallisto bustools pipeline as described in Wolfien et al. 2020 (https://doi.org/10.3390/cells9020318) script also available in FairdomHub
setwd(".../bus_output_bl6")  	# Single-nuclei raw data of the BL6 mice strain processed with the kallisto bustools pipeline as described Wolfien et al. 2020 (https://doi.org/10.3390/cells9020318) script also available in FairdomHub
data_dir_vid <- "../cellranger_vidal/outs/filtered_feature_bc_matrix"	# Single-nuclei data from Vidal et al. 2019 (https://insight.jci.org/articles/view/131092) obtained as raw files and processed with CellRanger (v.6.0.0.1) on default settings
data_dir_lin <- "../_csv_lin"	# Single-cell data from Linscheid et al. 2019 (https://www.nature.com/articles/s41467-019-10709-9) obtained as normalized count files in .csv format 


# For the first two datasets, we have to specifically integrate the spliced and unspliced matrices that have been generated as standard output formats via kallisto and bustools to be read into R via BusParse
# Load and save them one after each other

c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = "./spliced",
                                                spliced_name = "spliced",
                                                unspliced_dir = "./unspliced",
                                                unspliced_name = "unspliced")

sum(unspliced@x) / (sum(unspliced@x) + sum(spliced@x))

dim(spliced)
dim(unspliced)

tot_count <- Matrix::colSums(spliced)
summary(tot_count)

bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)

# Can only plot barcodes with both spliced and unspliced counts
bcs_inter <- intersect(colnames(spliced), colnames(unspliced))
s <- colSums(spliced[,bcs_inter])
u <- colSums(unspliced[,bcs_inter])
# Grid points
sr <- sort(unique(exp(round(log(s)*100)/100)))
ur <- sort(unique(exp(round(log(u)*100)/100)))

bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]

# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 1]

# Cells to use, same amount of cells needed, use the same number for both depending on the lower number
sf_bar <- unspliced@Dimnames[[2]][1:592725] # Fzt:DU
sf_bar <- unspliced@Dimnames[[2]][1:321641] # BL66

sf <- spliced[spliced@Dimnames[[1]] %in% genes_use, spliced@Dimnames[[2]] %in% sf_bar]
dim(sf)

# Now use the matching amount of cells also for the unspliced version and obtain a matching matrix for spliced/unspliced =)

uf_bar <- sf@Dimnames[[2]][1:321085] # fztDU
uf_bar <- sf@Dimnames[[2]][1:193155] # Bl6 193155

uf <- unspliced[unspliced@Dimnames[[1]] %in% genes_use, unspliced@Dimnames[[2]] %in% uf_bar]

# Control again both dims
dim(sf)
dim(uf)

# Create the Seurat object with the two assay types spliced and unspliced (during this analysis only the spliced version will be used)
# Improtant note: For a later combination load each dataset individually and load it into the matrix
seu.fzt <- CreateSeuratObject(sf, assay = "spliced")
seu.fzt[["unspliced"]] <- CreateAssayObject(uf)
seu.fzt$tech <- "Fzt:DU"
seu.fzt

seu.bl6 <- CreateSeuratObject(sf, assay = "spliced")
seu.bl6[["unspliced"]] <- CreateAssayObject(uf)
seu.bl6$tech <- "BL6"
seu.bl6

# Load the Vidal (vid) and Linscheid (lin) datasets as well and create the seurat objects
counts.vid <- Read10X(data.dir = data_dir_vid)
seu.vid <- CreateSeuratObject(counts = counts.vid, min.cells = 3, min.features = 200, project = "Vidal_snRNAseq")
seu.vid

counts.lin <-read.table(file=paste0("../_csv_lin/GSE130710_normdata.csv"),sep=",", header = TRUE, row.names = 1, as.is = TRUE)
seu.lin <- CreateSeuratObject(counts = counts.lin, min.cells = 3, min.features = 200, project = "Linscheid_scRNAseq")
seu.lin 

# Filter low quality reads and potential cell doublets. Also restricts the features to be able to run the computations on a usual desktop PC
seu.fzt <- subset(seu.fzt, subset = nFeature_spliced > 200 & nFeature_spliced < 3000 & nFeature_unspliced > 200 & nFeature_unspliced < 3000)
seu.bl6 <- subset(seu.bl6, subset = nFeature_spliced > 200 & nFeature_spliced < 3000 & nFeature_unspliced > 200 & nFeature_unspliced < 3000)
seu.vid <- subset(seu.vid, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
seu.lin  <- subset(seu.lin , subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

# Processing all datasets in the same manner
seu.fzt <- NormalizeData(seu.fzt, normalization.method = "LogNormalize", scale.factor = 10000)
# Seurat normalization global-scaling normalization method “LogNormalize” that normalizes the feature 
# expression measurements for each cell by the total expression
seu.fzt <- FindVariableFeatures(seu.fzt, selection.method = "vst", nfeatures = 3000)
seu.fzt <- ScaleData(seu.fzt, verbose = FALSE) 
# Note:
# Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. 
# Therefore, the default in `ScaleData` is only to perform scaling on the previously identified variable features (2,000 by default). 
# To do this, omit the `features` argument in the previous function call, i.e.'
seu.fzt <- RunPCA(seu.fzt, npcs = 30, verbose = FALSE) # The dimensionalty of the datasets have been determined before individually by Elbow and Jack-Straw
seu.fzt <- RunHarmony(seu.fzt, group.by.vars = "orig.ident", assay.use = "spliced")
seu.fzt <- RunUMAP(seu.fzt, reduction = "pca", dims = 1:30)
seu.fzt <- RunUMAP(seu.fzt, reduction = "harmony", dims = 1:30)

# The BL6 dataset was analysed in the same manner
seu.bl6 <- NormalizeData(seu.bl6, normalization.method = "LogNormalize", scale.factor = 10000)
seu.bl6 <- FindVariableFeatures(seu.bl6, selection.method = "vst", nfeatures = 3000)
seu.bl6 <- ScaleData(seu.bl6, verbose = FALSE)
seu.bl6 <- RunPCA(seu.bl6, npcs = 30, verbose = FALSE)
seu.bl6 <- RunHarmony(seu.bl6, group.by.vars = "orig.ident", assay.use = "spliced")
seu.bl6 <- RunUMAP(seu.bl6, reduction = "pca", dims = 1:30)
seu.bl6 <- RunUMAP(seu.bl6, reduction = "harmony", dims = 1:30)

seu.vid <- NormalizeData(seu.vid, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
seu.lin <- NormalizeData(seu.lin, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

##### Cluster the cells ######

# Seurat v4 applies a graph-based clustering approach, building upon initial strategies in ([Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8)). 
# Importantly, the *distance metric* which drives the clustering analysis (based on previously identified PCs) remains the same. 
# However, the approach to partioning the cellular distance matrix into clusters has dramatically improved. 
# Their approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [[SNN-Cliq, Xu and Su, Bioinformatics, 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract) and CyTOF data [[PhenoGraph, Levine *et al*., Cell, 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). 
# Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature 
# expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'. 

# As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells 
# based on the shared overlap in their local neighborhoods (Jaccard similarity). 
# This step is performed using the `FindNeighbors` function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

# To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [[SLM, Blondel *et al*., Journal of Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), 
# to iteratively group cells together, with the goal of optimizing the standard modularity function. 
# The `FindClusters` function implements this procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering, 
# with increased values leading to a greater number of clusters. 
# It is considered that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. The clusters can be found using the `Idents` function.

seu.fzt <- FindNeighbors(object = seu.fzt, dims = 1:30)
seu.fzt <- FindClusters(object = seu.fzt, resolution = 0.725)
pdf('fzt_clusters.pdf')
DimPlot(seu.fzt, reduction = "umap", label = TRUE)
dev.off()

seu.bl6 <- FindNeighbors(object = seu.bl6, dims = 1:25)
seu.bl6 <- FindClusters(object = seu.bl6, resolution = 0.725)
pdf('bl6_clusters.pdf')
DimPlot(seu.bl6, reduction = "umap", label = TRUE)
dev.off()

seu.vid <- FindNeighbors(object = seu.vid, dims = 1:30)
seu.vid <- FindClusters(object = seu.vid, resolution = 1.525)
DimPlot(seu.vid, reduction = "umap", label = TRUE)

seu.lin <- FindNeighbors(object = seu.lin, dims = 1:18)
seu.lin <- FindClusters(object = seu.lin, resolution = 0.925)
DimPlot(seu.vid, reduction = "umap", label = TRUE)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
fzt.markers <- FindAllMarkers(seu.fzt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fzt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100_fzt # Top 100 Markers
fzt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_fzt # Top 5 Markers

bl6.markers <- FindAllMarkers(seu.bl6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
bl6.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100_bl6 # Top 100 Markers
bl6.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10_bl6 # Top 10 Markers

vid.markers <- FindAllMarkers(seu.vid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
vid.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100_vid # Top 100 Markers
vid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_vid # Top 5 Markers

lin.markers <- FindAllMarkers(seu.lin, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lin.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100_lin # Top 100 Markers
lin.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_lin # Top 5 Markers

DotPlot(object = seu.fzt, features = rev(unique(top5_seu.fzt$gene)))+ RotatedAxis()

# For further usage in scSynO, we annotated and saved the identified top 100 markers in a .csv file
# They can serve as a refined input for the analysis (Tested for top 20, 50, and 100 transcripts)
# Annotate the murine transcipts
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# Collect gene names with
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

sc_anno_top100 <-  cbind(t2g$ext_gene[t2g$ens_gene %in% top100_fzt$gene],t2g$ens_gene[t2g$ens_gene %in% top100_fzt$gene])
colnames(sc_anno_top100) <- c("Gene","ID")

var = c ("gene","cluster","p_val_adj")

Marker_genes100_seu_har <-top100_fzt[var]
fin100 <- unique(merge(Marker_genes100_seu_har, sc_anno_top100, by.x = "gene", by.y = "ID"))

fin100_sorted <- fin100[order(fin100[, "cluster"]), drop = FALSE]

write.csv(fin100_sorted, file="Marker_genes100_fzt.csv") # Save all 100 Marker genes in a .csv file
top100_unique_seu_har <- unique(unlist(fin100_sorted$gene))
top100_unique_seu_har_anno <- rev(unique(unlist(fin100_sorted$Gene)))
plot1 <- DotPlot(object = seu.fzt, features = top100_unique_seu_har)
plot1 + theme(axis.text.x = element_text(angle = 45, siz = 6,hjust = 1))+ scale_x_discrete(labels= top100_unique_seu_har_anno)

# Use the dopdf function to save the current plot in the current directory
dopdf(filename = "top100markers_per_cluster.pdf", pdf.width = 90, pdf.height = 20)

# Use case 1 utilizes the manually identified cardiac glial cells as an input for the scSynO training step
# Here, these cells refer to cluster 7, which are individually pulled from the data matrix
# The rest of the data is stored in seperared matrix
# GetAssayData allows pulling from a specific slot rather than just data
cell_expression_7 <- as.matrix(GetAssayData(seu.fzt, slot = "data")[, WhichCells(seu.fzt, ident = "7")])
cell_expression <- as.matrix(GetAssayData(seu.fzt, slot = "data")[, WhichCells(seu.fzt, ident = c("0","1","2","3","4","5","6","8"))])

# Export the .csv files for further processing with scSYnO
write.csv(cell_expression_all, file="cell_expression_all.csv")
write.csv(cell_expression_7, file="cell_expression_7.csv")

# Check for glial cell expression markers in the validation datasets
DotPlot(object = seu.bl6, features = rev(unique(top5_bl6$gene)))+ RotatedAxis()
dopdf(filename = ".../glialcells/val1_dot.pdf", pdf.width = 30, pdf.height = 16)

DotPlot(object = seu.vid, features = rev(unique(top5_vid$gene)))+ RotatedAxis()
dopdf(filename = ".../glialcells/val2_dot.pdf", pdf.width = 30, pdf.height = 16)

# Indicate cells within UMAP plot, predicted cell from scSynO can be fed back and highligthed for inspection
top20 <- read.table(".../glialcells/top_20/predict_val1_20_norm.txt", header = FALSE, sep = "") # LoRAS predictions of glial cells for the top 20 features
top50 <- read.table(".../glialcells/top_50/predict_val1_50_norm.txt", header = FALSE, sep = "") # LoRAS predictions of glial cells for the top 50 features 
DimPlot(seu.bl6, label=T, cells.highlight= list(top20$V1, top50$V1), cols.highlight = c("darkblue", "darkred"), cols= "grey")
DimPlot(seu.bl6, label=T, cells.highlight= top20$V1, cols.highlight = "darkblue", cols= "grey")

top20 <- read.table(".../glialcells/top_20/predict_val2_20_norm.txt", header = FALSE, sep = "") # LoRAS predictions of glial cells for the top 20 features
top50 <- read.table(".../glialcells/top_50/predict_val2_50_norm.txt", header = FALSE, sep = "") # LoRAS predictions of glial cells for the top 50 features 
DimPlot(seu.vid, label=T, cells.highlight= list(top20$V1, top50$V1), cols.highlight = c("darkblue", "darkred"), cols= "grey")
DimPlot(seu.vid, label=T, cells.highlight= top50$V1, cols.highlight = "darkblue", cols= "grey")

# Integration of both mice strains for the identification of the proliverative cardiomyocytes (use case 2)
anchors <- FindIntegrationAnchors(object.list = c(seu.fzt,seu.bl6), dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Run again the same analyses for the combined data, please be aware of the changing nuclei amount
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident", assay.use = "integrated")
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)

integrated <- FindNeighbors(object = integrated, dims = 1:30)
integrated <- FindClusters(object = integrated, resolution = 0.525)

# Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. 
# The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. 
# Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. 
# As input to the UMAP and tSNE, it is suggested using the same PCs as input to the clustering analysis.

# Please, note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1 <- DimPlot(integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# find markers for every cluster compared to all remaining cells, report only the positive ones
integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
integrated.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100_integrated # Top 100 Markers per cluster

# Again annotate the identified transcripts (only needed for analyses with kallisto bustools)
sc_anno_top100_integrated <-  cbind(t2g$ext_gene[t2g$ens_gene %in% top100_integrated$gene],t2g$ens_gene[t2g$ens_gene %in% top100_integrated$gene])
colnames(sc_anno_top100_integrated) <- c("Gene","ID")

Marker_genes100_integrated <- sc_anno_top100_integrated[var]
fin100_integrated <- unique(merge(Marker_genes100_integrated, sc_anno_top100_integrated, by.x = "gene", by.y = "ID"))

fin100_sorted_integrated <- fin100_integrated[order(fin100_integrated[, "cluster"]), drop = FALSE]

write.csv(fin100_sorted_integrated, file="Marker_genes100_integrated.csv") # Save all 100 Marker genes in a .csv file
top100_unique_integrated <- unique(unlist(fin100_sorted_integrated$gene))
top100_unique_seu_har_anno_integrated <- rev(unique(unlist(fin100_sorted_integrated$Gene)))
plot1 <- DotPlot(object = integrated, features = top100_unique_integrated)
plot1 + theme(axis.text.x = element_text(angle = 45, siz = 6,hjust = 1))+ scale_x_discrete(labels= top100_unique_seu_har_anno_integrated)

# Indicate cells within umap plot

top20 <- read.table(".../proCM/top 20/predict_val2_21.txt", header = FALSE, sep = "") # scSynO predictions of glial cells for the top 20 features
top50 <- read.table(".../proCM/top 100/predict_val1_100.txt", header = FALSE, sep = "") # scSynO predictions of glial cells for the top 50 features 
DimPlot(seu.vid, label=T, cells.highlight= list(top20$V1, top50$V1), cols.highlight = c("darkblue", "darkred"), cols= "grey")
DimPlot(seu.lin, label=T, cells.highlight= top20$V1, cols.highlight = "darkblue", cols= "grey")


########## Characterizing the clusters ##########
# Average mean number of counts for each cell you can run

cluster.averages <- AverageExpression(integrated)
head(cluster.averages[["integrated"]][, 1:5])
colSums(cluster.averages[["integrated"]])

# Average mean number of counts
summary(colSums(integrated))
summary(colSums(seu.fzt))
summary(colSums(seu.bl6))
summary(colSums(seu.vid))
summary(colSums(seu.lin))

# Store cluster identities in object@meta.data$my.clusters
object_fzt<- Seurat::StashIdent(object = seu.fzt, save.name = "my.clusters")
object_bl6 <- Seurat::StashIdent(object = seu.bl6, save.name = "my.clusters")
object_integrated <- Seurat::StashIdent(object = integrated, save.name = "my.clusters")

object_vid <- Seurat::StashIdent(object = seu.vid, save.name = "my.clusters")
object_lin <- Seurat::StashIdent(object = seu.lin, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object_fzt@meta.data$my.clusters, object_fzt@meta.data$orig.ident)
table(object_bl6@meta.data$my.clusters, object_bl6@meta.data$orig.ident)
table(object_integrated@meta.data$my.clusters, object_integrated@meta.data$orig.ident)

table(object_vid@meta.data$my.clusters, object_vid@meta.data$orig.ident)
table(object_lin@meta.data$my.clusters, object_lin@meta.data$orig.ident)

sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux bullseye/sid

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.13.so

# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base

#other attached packages:
# [1] velocyto.R_0.6     Matrix_1.2-18      harmony_1.0        Rcpp_1.0.5
# [5] devtools_2.3.2     usethis_1.6.3      ggplot2_3.3.2      SeuratObject_4.0.0
# [9] Seurat_4.0.0       dplyr_1.0.2

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.15           colorspace_2.0-0     deldir_0.2-3
#  [4] ellipsis_0.3.1       ggridges_0.5.2       rprojroot_1.3-2
#  [7] fs_1.5.0             spatstat.data_2.1-0  leiden_0.3.3
# [10] listenv_0.8.0        remotes_2.2.0        ggrepel_0.8.2
# [13] fansi_0.4.1          codetools_0.2-18     splines_4.0.3
# [16] polyclip_1.10-0      pkgload_1.1.0        jsonlite_1.7.1
# [19] ica_1.0-2            cluster_2.1.0        png_0.1-7
# [22] uwot_0.1.10          shiny_1.5.0          sctransform_0.3.2
# [25] compiler_4.0.3       httr_1.4.2           backports_1.1.8
# [28] assertthat_0.2.1     fastmap_1.0.1        lazyeval_0.2.2
# [31] cli_2.2.0            later_1.1.0.1        htmltools_0.5.1.1
# [34] prettyunits_1.1.1    tools_4.0.3          igraph_1.2.6
# [37] gtable_0.3.0         glue_1.4.2           RANN_2.6.1
# [40] reshape2_1.4.4       rappdirs_0.3.1       spatstat_1.64-1
# [43] Biobase_2.50.0       scattermore_0.7      vctrs_0.3.6
# [46] nlme_3.1-151         lmtest_0.9-37        stringr_1.4.0
# [49] globals_0.14.0       ps_1.4.0             testthat_2.3.2
# [52] mime_0.9             miniUI_0.1.1.1       lifecycle_0.2.0
# [55] irlba_2.3.3          goftest_1.2-2        future_1.21.0
# [58] MASS_7.3-53          zoo_1.8-8            scales_1.1.1
# [61] pcaMethods_1.80.0    promises_1.1.1       spatstat.utils_2.1-0
# [64] parallel_4.0.3       RColorBrewer_1.1-2   curl_4.3
# [67] memoise_1.1.0        reticulate_1.16      pbapply_1.4-2
# [70] gridExtra_2.3        rpart_4.1-15         stringi_1.5.3
# [73] desc_1.2.0           BiocGenerics_0.36.0  pkgbuild_1.1.0
# [76] rlang_0.4.10         pkgconfig_2.0.3      matrixStats_0.57.0
# [79] lattice_0.20-41      ROCR_1.0-11          purrr_0.3.4
# [82] tensor_1.5           patchwork_1.0.1      htmlwidgets_1.5.2
# [85] cowplot_1.0.0        tidyselect_1.1.0     processx_3.4.5
# [88] parallelly_1.21.0    RcppAnnoy_0.0.18     plyr_1.8.6
# [91] magrittr_2.0.1       R6_2.5.0             generics_0.1.0
# [94] pillar_1.4.7         withr_2.3.0          mgcv_1.8-33
# [97] fitdistrplus_1.1-1   survival_3.2-7       abind_1.4-5
#[100] tibble_3.0.4         future.apply_1.5.0   crayon_1.3.4
#[103] KernSmooth_2.23-18   plotly_4.9.2.1       grid_4.0.3
#[106] data.table_1.13.2    callr_3.5.1          digest_0.6.27
#[109] xtable_1.8-4         tidyr_1.1.2          httpuv_1.5.4
#[112] munsell_0.5.0        viridisLite_0.3.0    sessioninfo_1.1.1
