# sc-SynO
single-cell Synthetic Oversampling
- The reprository contains python codes for implemetation of sc-SynO through the [LoRAS](https://github.com/narek-davtyan/LoRAS) algorithm.
- We provide the analysis results in form of Jupyter notebooks.
- R scripts that are used to pre-process the data are made available.
- Two separate folders show the results for our analysis on Gial cell and Proliferative cardiomyocite detection respectively.
- The data files used for our benchmark are publicly accessible at via ArrayExpress (E-MTAB-7869, E-MTAB-8751, E-MTAB-8848), the [Allen Brain Atlas](https://celltypes.brain-map.org/), and GEO (GSE130710).
- A preprint of sc-SynO is available at [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.20.427486v1)

## Why sc-SynO?
The research landscape of single-cell and single-nuclei RNA-sequencing is evolving rapidly. In particular, the area for the detection of rare cells was highly facilitated by this technology. However, an automated, unbiased, and accurate annotation of rare subpopulations is challenging. Once rare cells are identified in one dataset, it is usually necessary to generate further specific datasets to enrich the analysis (e.g., with samples from other tissues). From a machine learning perspective, the challenge arises from the fact that rare-cell subpopulations constitute an imbalanced classification problem. 

We here introduce a Machine Learning (ML)-based oversampling method that uses gene expression counts of already identified rare cells as an input to generate synthetic cells to then identify similar (rare) cells in other publicly available experiments. We utilize single-cell synthetic oversampling (sc-SynO), which is based on the Localized Random Affine Shadowsampling (LoRAS) algorithm. The algorithm corrects for the overall imbalance ratio of the minority and majority class. 

## How to use?
sc-SynO can be used within all commonly used single-cell data analysis pipeline. The input files are *count.csv* files from raw or normalized (e.g., log, or SCT) data of the target rare-cell cluster (starting from 3 cells) and reference cells. The training can be done on all processed features or already identified marker genes for the rare-cell type of interest. Identified markers can be obtained elsewhere from in-build algorithms of Seurat, ML-based approaches, databases or domain expertise.

## What are possible use cases?
 Applying sc-SynO on a novel dataset to identify the same rare-cell type is magnitudes less time-consuming than manually curated data processing and annotation of scRNA-Seq data. These facilitated cell enrichment can be used for more in-depth downstream analyses with the cell type of interest, without re-analyzing all the datasets.  Such a scenario can be of high interest for single-cell identification in cancer, hypothesis testing on larger cell sets, cell homology search across tissues, or other individual applications.