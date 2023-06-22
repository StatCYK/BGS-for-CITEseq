# Name: realdata1.R
# Function: Comparing with Seurat & scanpy
# Note: Using new normalization, and not removing y==0
# ------------------------------------------------------------------------------
# rm(list = ls())

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(dbscan)
library(leaps)
library(caret) # PCA
library(jackstraw)
library(paletteer)
# BiocManager::install("MAST")
# BiocManager::install("DESeq2")

setwd("~/Desktop/rs_BGS/codes/")
set.seed(123)
t1 = Sys.time()
# ------- Functions ------------------------------------------------------------
clr_function = function(x) { # For normalizing ADT (x)
  return(log(x = (1+x)/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
}

r_sq <- function(predicted, actual) { # For computing R-squared
  RSS <- sum((actual - predicted)^2)
  TSS <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (RSS / TSS)
  return(R_squared)
}

a_r_sq<- function(predicted, actual, p) {
  RSS <- sum((actual - predicted)^2)
  TSS <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (RSS / TSS)
  n <- length(actual)
  adjusted_R_squared <- 1 - ((1 - R_squared) * (n - 1) / (n - p - 1))
  return(adjusted_R_squared)
}


# ------ read and preprocess ---------------------------------------------------
cbmc.rna = as.sparse(read.csv(file = "../data/cbmc_CITEseq/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))
cbmc.rna = CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc.adt = as.sparse(read.csv(file = "../data/cbmc_CITEseq/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))
cbmc = CreateSeuratObject(counts = cbmc.rna)
adt_assay = CreateAssayObject(counts = cbmc.adt)
cbmc[["ADT"]] = adt_assay

## normalizing ADT (x) ##
ADT.data0 = scale(apply(cbmc@assays$ADT@counts, 2, clr_function))
ADT.data  = data.frame(t(ADT.data0))
dim(ADT.data)

## normalizing RNA (y) ##
DefaultAssay(cbmc) = 'RNA'
cbmc = NormalizeData(cbmc) %>% FindVariableFeatures() %>% ScaleData()

CD14 = cbmc["rna_CD14"]
y = matrix(CD14@assays$RNA@data)

### Saving data ###
data_8000 = cbind(y,ADT.data)
row.names(data_8000) = NULL
# write.csv(data_8000,"../data/data_8000.csv", row.names = F)

### Saving un-normalized data ###
unnormalized_ADT = t(as.matrix(cbmc@assays$ADT@counts))
data_8000_unnormalized = as.data.frame(cbind(y, unnormalized_ADT))
row.names(data_8000_unnormalized) = NULL

# ----- Settings ---------------------------------------------------------------

num_sims = 100

# ----- Main -------------------------------------------------------------------
dat = read.csv("../data/datanew/data_8000_new.csv")
data_8000_unnormalized = read.csv("../data/datanew/data_8000_unnormalized.csv")
dat = dat[,-1]
data_8000_unnormalized = data_8000_unnormalized[,-1]
head(dat)
p = ncol(dat) - 1

# === Main

DE_method_list = c( "wilcox", "bimod", "roc", "t", "LR", "MAST")
SCAN_method_list = c('logreg', 't-test', 'wilcoxon')

method_names = c('BSS', paste0("Seurat_", DE_method_list),
                 SCAN_method_list,
                 "Corr", "Corr_ABS")

num_methods = length(method_names)

num_selected_results = matrix(NA, nrow = num_sims, ncol = num_methods)
MSE_results   = matrix(NA, nrow = num_sims, ncol = num_methods)
R2_results    = matrix(NA, nrow = num_sims, ncol = num_methods)
adjR2_results = matrix(NA, nrow = num_sims, ncol = num_methods)
colnames(MSE_results) = colnames(R2_results) = colnames(adjR2_results) = 
  colnames(num_selected_results) = method_names

truncR2_res_array = array(NA, dim = c(num_sims,num_methods, p))
MSE_res_array     = array(NA, dim = c(num_sims,num_methods, p))
pearson_res_array = array(NA, dim = c(num_sims,num_methods, p))
spearman_res_array= array(NA, dim = c(num_sims,num_methods, p))
Tau_res_array     = array(NA, dim = c(num_sims,num_methods, p))
for(case_sim in 1: num_sims){
  print(paste0("case_sim: ", case_sim," / ",num_sims))
  set.seed(case_sim)
  n_rows = nrow(dat)
  row_indices = sample(1:n_rows)
  sep_index = round(0.75 * n_rows)
  train_indices = row_indices[1:sep_index]
  test_indices = row_indices[(sep_index+1):n_rows]
  
  train = dat[train_indices, ]
  test = dat[test_indices, ]
  
  train_unnormal = data_8000_unnormalized[train_indices,]
  
  
  write.csv(train_indices, file = paste0("../indices/real1_sim_",case_sim ,"_train.csv"))
  write.csv(test_indices, file = paste0("../indices/real1_sim_",case_sim ,"_test.csv"))
  
  method_id = 0
  # ========= Method 1: BSS ========= 
  method_id = method_id + 1
  best_subset = regsubsets(y ~ ., data = train, 
                           nvmax = 13, method = "exhaustive")
  best_summary = summary(best_subset)
  lowest_bic = which.min(best_summary$bic)
  best_model_vars0 = as.character(best_summary$which[lowest_bic, ])
  best_model_vars =colnames(train)[best_model_vars0 == TRUE][-1]
  
  fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
  fit = lm(fmla, data = train)
  y_hat = predict(fit,newdata = test)
  
  num_selected_results[case_sim, method_id] = length(best_model_vars)
  MSE_results[case_sim, method_id]   = mean((y_hat - test$y)^2)
  R2_results[case_sim, method_id]    = r_sq(y_hat, test$y)
  adjR2_results[case_sim, method_id] = a_r_sq(y_hat, test$y, p = length(best_model_vars))
  
  for(trun_num in 1:p){
    
    best_subset = regsubsets(y ~ ., data = train,
                             nvmax = trun_num, method = "exhaustive")
    best_summary = summary(best_subset)
    
    best_model_vars0 <- best_summary$which[(trun_num), , drop = FALSE]
    
    
    # lowest_bic = which.min(best_summary$bic)
    
    # best_model_vars0 = as.character(best_summary$which[lowest_bic, ])
    best_model_vars =colnames(train)[best_model_vars0 == TRUE][-1]
    
    fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
    fit = lm(fmla, data = train)
    y_hat = predict(fit,newdata = test)
    
    MSE_res_array[case_sim, method_id, trun_num] = mean((y_hat - test$y)^2)
    truncR2_res_array[case_sim, method_id, trun_num] = r_sq(y_hat, test$y)
    pearson_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "pearson")
    spearman_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "spearman")
    Tau_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "kendall")
  } # case_trun_num
  
  # ========= Method 2: Seurat_DE methods =========

  for(DE_method in DE_method_list){
    method_id = method_id + 1
    
    train_for_DE = train

    seurat_object = CreateSeuratObject(counts = t(train_for_DE[,-1]))
    my_labels = as.character(as.numeric(train$y==0))

    Idents(object = seurat_object) =  as.factor(my_labels)
    marker_result = FindMarkers(seurat_object, ident.1 = "1", ident.2 = "0",
                                test.use = DE_method) # ("wilcox" is default)
    best_model_vars = rownames(marker_result)

    fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
    fit = lm(fmla, data = train)
    y_hat = predict(fit,newdata = test)

    num_selected_results[case_sim, method_id] = length(best_model_vars)
    MSE_results[case_sim, method_id]   = mean((y_hat - test$y)^2)
    R2_results[case_sim, method_id]    = r_sq(y_hat, test$y)
    adjR2_results[case_sim, method_id] = a_r_sq(y_hat, test$y, p = length(best_model_vars))


    for(trun_num in 1:p){
      num_vars_to_use = min(nrow(marker_result), trun_num)
      best_model_vars = rownames(marker_result)[1:num_vars_to_use]

      fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
      fit = lm(fmla, data = train)
      y_hat = predict(fit, newdata = test)

      MSE_res_array[case_sim, method_id, trun_num] = mean((y_hat - test$y)^2)
      truncR2_res_array[case_sim, method_id, trun_num] = r_sq(y_hat, test$y)
      pearson_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "pearson")
      # spearman_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "spearman")
      # Tau_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "kendall")
    } # trun_num
  }
    # ========= Method 3: SCANPY methods =========
    for(SCAN_method in SCAN_method_list){
      method_id = method_id + 1

      SCAN_file = paste0("../scanpy_results/real1_sim", case_sim, "_",
                         SCAN_method, ".csv")

      for(trun_num in 1:p){


        scan_res = read.csv(SCAN_file)

        num_vars_to_use = min(nrow(scan_res), trun_num)
        best_model_vars0 = scan_res[1:num_vars_to_use, 1]
        best_model_vars = colnames(train_for_DE[,-1])[colnames(train_for_DE[,-1]) %in% best_model_vars0]


        fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
        fit = lm(fmla, data = train)
        y_hat = predict(fit, newdata = test)

        MSE_res_array[case_sim, method_id, trun_num] = mean((y_hat - test$y)^2)
        truncR2_res_array[case_sim, method_id, trun_num] = r_sq(y_hat, test$y)
        pearson_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "pearson")
        # spearman_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "spearman")
        # Tau_res_array[case_sim, method_id, trun_num] = cor(y_hat, test$y, method = "kendall")
      } # trun_num

  } # SCAN_method
  
  
 
  
}# case_sim


save.image(file = "plot_results_realdata1.rda")

# load("plot_results_realdata1.rda")

# # ===== Plotting results =====

{
  # plotting_metric = "Pearson correlation"
  plotting_metric = "MSE"
  
  if(plotting_metric == "Pearson correlation"){ylim = c(0.4, 0.65)}
  if(plotting_metric == "MSE"){ylim = c(0.054, 0.07)}
  
par(mfrow = c(1, 1))



SCAN_method_list = c('logreg', 't-test', 'wilcoxon', 't-test_var' )
ploting_methods = c("DESeq2", "wilcox", "bimod", "roc", "t", "poisson", "LR", "MAST")
ploting_methods = c('BSS', paste0("Seurat_", ploting_methods),
                    paste0("SCANPY_", SCAN_method_list), "Corr", "Corr_ABS")


random_per_sd = 0.0
lwd = 2
# === Reading quantum results

quantum_res = read.csv("GBS_select 1.csv")
quantum_res = quantum_res[,-1]
colMeans(quantum_res)

dim(quantum_res)


# 
if(plotting_metric == "Pearson correlation"){current_array = (pearson_res_array)}
if(plotting_metric == "MSE"){current_array = (MSE_res_array)}

dim(current_array)



colors =  paletteer_d("tvthemes::kimPossible")
col_set = c('red', colors[-4])
col_set[11] = 'blue'

mean_array <- apply(current_array, c(2, 3), mean)
sd_array <- apply(current_array, c(2, 3), sd)
rownames(mean_array) = ploting_methods
rownames(sd_array) = ploting_methods
colnames(mean_array) = paste0('num', 1:13)

x_range = 1:12
col_id = 1
par(mar = c(5, 5, 5, 5))

# mean_array[10,2] = mean_array[1,2]
## plot BSS
random_per = random_per_sd
plot(x_range + random_per, mean_array[1,x_range], 
     col = col_set[col_id], 
     pch = col_id,
     ylim = ylim,
     xlim = range(x_range), #range(x_range),
     lwd = lwd,
     ylab = plotting_metric, xlab = "Selected number of variables",
     cex.lab = 1.5,
     cex.axis = 1.5, type = "b")


## plot BGS
col_id = 1
random_per = random_per + random_per_sd
BGS_mean = colMeans(quantum_res)
BGS_sd = apply(quantum_res, 2, sd)
lines(x_range + random_per, BGS_mean[x_range],
      pch = col_id,
      col = col_set[col_id],
      lwd = lwd, type = "n")


## plot Seurat
for(case_method in ploting_methods[c(3,4,5,6,8,9)]){
  col_id = col_id + 1
  random_per = random_per + random_per_sd
  lines(x_range + random_per, mean_array[case_method,x_range], 
        pch = col_id,
        col = col_set[col_id], 
        lwd = lwd, type = "b")
}

## plot Others
for(case_method in ploting_methods[c(10,11,12)]){
  
  col_id = col_id + 1
  random_per = random_per + random_per_sd
  lines(x_range + random_per, mean_array[case_method, x_range], 
        pch = col_id,
        col = col_set[col_id], 
        lwd = lwd, type = "n")
}



# fin
actual_methods = ploting_methods[c(1,
                                   3,4,5,6,8,9,
                                   10,11,12,14)]
# legend('topright',legend = actual_methods[1:6],
#        col=col_set[1:6], lty = 1,
#        pch = 1:6,
#        cex = 1.5, lwd = lwd, horiz = F)
# #
# legend('bottomright',legend = actual_methods[7:10],
#        col=col_set[7:11], lty = 1,
#        pch = 7:11,
#        cex = 1.5, lwd = lwd, horiz = F)

}


print(Sys.time() - t1)













