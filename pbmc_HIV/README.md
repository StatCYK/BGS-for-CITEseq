# Panel design for immune cell type identification

## Run experiments

### CITE-seq data preprocessing ####
1. You need to download `pbmc_multimodal.h5seurat` at https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
2. Upload the files to the folder `/data/`

3. Run `marker_data_generate.R`

### Fit logistic regression model & calculate BICs, AUCs, specificity, sensitivity ###

Run `BIC_calculate.R`


### Analysis with BGS algorithm ###

run `BGS_select_fix_size.py`
