# sc-SynO
single-cell Synthetic Oversampling
- The reprository contains python codes for implemetation of sc-SynO through the [LoRAS](https://github.com/narek-davtyan/LoRAS) algorithm.
- We provide the analysis results in form of Jupyter notebooks.
- R scripts that are used to pre-process the data are made available.
- Two separate folders show the results for our analysis on Gial cell and Proliferative cardiomyocite detection respectively.
- The raw data files used for our benchmark are publicly accessible via ArrayExpress (E-MTAB-7869, E-MTAB-8751, E-MTAB-8848), the [Allen Brain Atlas](https://celltypes.brain-map.org/), and GEO (GSE130710).
- A preprint of sc-SynO is available at [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.20.427486v1)

## Why sc-SynO?
The research landscape of single-cell and single-nuclei RNA-sequencing is evolving rapidly. In particular, the area for the detection of rare cells was highly facilitated by this technology. However, an automated, unbiased, and accurate annotation of rare subpopulations is challenging. Once rare cells are identified in one dataset, it is usually necessary to generate further specific datasets to enrich the analysis (e.g., with samples from other tissues). From a machine learning perspective, the challenge arises from the fact that rare-cell subpopulations constitute an imbalanced classification problem. 

We here introduce a Machine Learning (ML)-based oversampling method that uses gene expression counts of already identified rare cells as an input to generate synthetic cells to then identify similar (rare) cells in other publicly available experiments. We utilize single-cell synthetic oversampling (sc-SynO), which is based on the Localized Random Affine Shadowsampling (LoRAS) algorithm. The algorithm corrects for the overall imbalance ratio of the minority and majority class. 

## How to use?
sc-SynO can be used within all commonly available single-cell data analysis pipelines. The input files are *count.csv* files from raw or normalized (e.g., log, or SCT) data of the target rare-cell cluster (starting from 3 cells) and reference cells. The training can be done on all processed features or already identified marker genes for the rare-cell type of interest. Identified markers can be obtained elsewhere from in-build algorithms of Seurat, ML-based approaches, databases or domain expertise.

### Identify and obtain rare cells
Here, a brief example to extract cluster information of an R Seurat object for rare-cell types of interest and the reference cells for training the model. In general, all methods are compatible with sc-SynO, since there are only two separted arrays for the rare and reference cells needed.

```r
# GetAssayData allows pulling from a specific slot rather than just data. Here in cluster 7 are our rare, cells and the other cells specify our reference cell set.
cell_expression_rare_cells <- as.matrix(GetAssayData(seurat_object, slot = "data")[, WhichCells(seurat_object, ident = "7")])
cell_expression_reference <- as.matrix(GetAssayData(seurat_object, slot = "data")[, WhichCells(seurat_object, ident = c("0","1","2","3","4","5","6","8"))])

# Export the .csv files for further processing with sc-SynO in Python
write.csv(cell_expression_rare_cells, file="cell_expression_rare_cells.csv")
write.csv(cell_expression_reference, file="cell_expression_reference.csv")
```

### The core of sc-SynO is LoRAS
Localized Randomized Affine Shadowsampling (LoRAS) oversampling technique

#### Installation
The latest version is available on PyPi and installable with the command: `pip install loras`

#### Usage
There is just one method `fit_resample(maj_class_points, min_class_points, k, num_shadow_points, list_sigma_f, num_generated_points, num_aff_comb, random_state=42)`

There are two mandatory inputs:  

- `maj_class_points` : Majority class parent data points, which is a non-empty list containing numpy arrays acting as points (reference cells)
- `min_class_points` : Minority class parent data points, which is a non-empty list containing numpy arrays acting as points (rare cells) 

There are also optional parameters:

- `k` : Number of nearest neighbours to be considered per parent data point (default value: `8 if len(min_class_points)<100 else 30`)
- `num_shadow_points` : Number of generated shadowsamples per parent data point (default value: `ceil(2*num_aff_comb / k)`)
- `list_sigma_f` : List of standard deviations for normal distributions for adding noise to each feature (default value: `[0.005, ... , 0.005]`)
- `num_generated_points` : Number of shadow points to be chosen for a random affine combination (default value: `ceil((len(maj_class_points) + len(min_class_points)) / len(min_class_points))`)
- `num_aff_comb` : Number of generated LoRAS points for each nearest neighbours group (default value: `min_class_points.shape[1]`)   


### Visualization of identified cells in R
After applying sc-SynO in Python it is possible to check and refine the results with the original data. For this one has to integrate the obtained cell IDs with the following commands:

```r
# Indicate cells within the UMAP plot, predicted cells from sc-SynO can be fed back and highligthed for inspection.
top20 <- read.table(".../predict_val_20.txt", header = FALSE, sep = "") # sc-SynO predictions for the top 20 features
top50 <- read.table("...predict_val_50.txt", header = FALSE, sep = "") # sc-SynO predictions for the top 50 features 

# Highlight one set of cells
DimPlot(seurat_object, label=T, cells.highlight= top20$V1, cols.highlight = "darkblue", cols= "grey")
DimPlot(seurat_object, label=T, cells.highlight= top50$V1, cols.highlight = "darkblue", cols= "grey")

# Or highlight two or more sets of cells in the same UMAP plot, e.g., to see the effect of different feature amounts
DimPlot(seurat_object, label=T, cells.highlight= list(top20$V1, top50$V1), cols.highlight = c("darkblue", "darkred"), cols= "grey") 
```


## What are possible use cases?
 Applying sc-SynO on a novel dataset to identify the same rare-cell type is magnitudes less time-consuming than manually curated data processing and annotation of scRNA-Seq data. These facilitated cell enrichment can be used for more in-depth downstream analyses with the cell type of interest, without re-analyzing all the datasets.  Such a scenario can be of high interest for single-cell identification in cancer, hypothesis testing on larger cell sets, cell homology search across tissues, or other individual applications.