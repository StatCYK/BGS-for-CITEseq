# RNA-ADT co-expression network with CBMC CITE-seq data

## Run experiments

### CITE-seq data preprocessing ####
1. You need to download `GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz` and `GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz` at the NCBI Gene Expression Omnibus (GEO; https://www.ncbi.nlm.nih.gov/geo/).

2. Upload the files to the folder `/data/`

3. Run `data preprocessing.R`

### Fit regression model & calculate BICs ###

Run `BIC_calculate.R`

### Base algorithms of Seurat & Scanpy for comparison ###

For Seurat, run `seurat_compare.R`
For Scanpy, run `scanpy.ipynb`

### Analysis with BGS algorithm ###

For markers with fixed sizes, run `BGS_select_fix_size.py`

For markers with all sizes, run `BGS_select.py`

For results visualization, run `Analysis_BGS.ipynb`
