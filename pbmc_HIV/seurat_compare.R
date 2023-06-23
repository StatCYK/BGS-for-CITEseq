# Name: seurat_compare.R
# Function: Compare with Seurat and scanpy
# ------------------------------------------------------------------------------
# rm(list = ls())

library(ggplot2)
library(scales)
library(leaps)
library(MASS)
library(ggplot2)
library(cowplot)
library(dplyr)

library(foreach)
library(doParallel)

set.seed(123)
t1 = Sys.time()
# ------------------------------------------------------------------------------


RNA_namelist = c("B intermediate kappa", "B intermediate lambda", "CD4 TEM_3",
                 "B memory kappa", "B memory lambda", "B naive kappa",
                 "B naive lambda", "CD4 CTL", "CD4 Naive", "CD4 TCM_1",
                 "CD4 TCM_2", "CD4 TCM_3", "CD4 TEM_1", "CD8 Naive", "CD8 TCM_1",
                 "CD8 TCM_2", "CD8 TCM_3", "CD8 TEM_1", "CD8 TEM_2", "Treg Memory",
                 "Treg Naive", "CD8 TEM_4", "CD8 TEM_5", "CD8 TEM_6", "CD14 Mono",
                 "CD16 Mono", "cDC2_1", "cDC2_2", "gdT_1", "gdT_2", "gdT_3",
                 "gdT_4", "MAIT", "NK Proliferating", "NK_1", "NK_2", "NK_3",
                 "NK_4", "NK_CD56bright", "pDC", "Platelet")

ttt = 0

cl <- makeCluster(8)
registerDoParallel(cl)
# for(RNA_name in RNA_namelist){
foreach(RNA_name=RNA_namelist) %dopar% {
  library(Seurat)
dat0 = read.csv(paste0("../data/data/", RNA_name, "_160000.csv"))
dim(dat0)
dat = dat0[,2:12]
colnames(dat)
dim(dat)
dat = dat[,c(11,1:10)]
colnames(dat)
dim(dat)

p = ncol(dat) - 1

FP_mat_list <- vector("list", p)
FN_mat_list <- vector("list", p)
AUC_mat_list <- vector("list", p)

for(trun_num in 1:(p-1)){
FP_mat_list[[trun_num]] <- read.csv(paste0("../data/result_BIC/results_", RNA_name, "/num", trun_num, "_test_error1.csv"))
FN_mat_list[[trun_num]] <- read.csv(paste0("../data/result_BIC/results_", RNA_name, "/num", trun_num, "_test_error2.csv"))
AUC_mat_list[[trun_num]] <- read.csv(paste0("../data/result_BIC/results_", RNA_name, "/num", trun_num, "_test_auc.csv"))
}

# ----- Settings ---------------------------------------------------------------
num_sims = 30

# ----- Main -------------------------------------------------------------------

DE_method_list = c("wilcox", "bimod", "roc", "t", "LR", "MAST")


method_names = c('BSS', paste0("Seurat_", DE_method_list))

SCAN_method_list = c('logreg', 't-test', 'wilcoxon')

method_names = c(method_names, paste0("Scanpy_", SCAN_method_list))

num_methods = length(method_names)

MSE_res_array = array(NA, dim = c(num_sims,num_methods, p-1))
FP_res_array  = array(NA, dim = c(num_sims,num_methods, p-1))
FN_res_array  = array(NA, dim = c(num_sims,num_methods, p-1))
AUC_res_array = array(NA, dim = c(num_sims,num_methods, p-1))

dim(MSE_res_array)
for(case_sim in 1: num_sims){
  
  print(paste0("case_sim: ", case_sim," / ",num_sims))
  print(Sys.time() - t1)
  
  # print(case_subset)
  set.seed(123*case_sim)
  n_rows = nrow(dat)
  row_indices <- sample(1:n_rows)
  sep_index <- round(0.75 * n_rows)
  train <- dat[row_indices[1:sep_index], ]
  test <- dat[row_indices[(sep_index+1):n_rows], ]
  dim(train)
  
  # =========== 1 .BSS method
  method_id = 1
  for(trun_num in 1:(p-1)){

    result_name = paste0("../data/result_BIC/results_", RNA_name, "/num",
                         trun_num,"_BIC.csv")
    BIC_mat = read.csv(result_name)
    BIC_colmeans = colMeans(BIC_mat)

    best_model_vars = colnames(BIC_mat)[which.min(BIC_colmeans)[1]]
    result_name_colname = paste(best_model_vars, collapse = "...")

    FP_mat  = FP_mat_list[[trun_num]]
    FN_mat  = FN_mat_list[[trun_num]]
    AUC_mat = AUC_mat_list[[trun_num]]
    
    FP_res_array[case_sim, method_id, trun_num]  = FP_mat[case_sim,result_name_colname]
    FN_res_array[case_sim, method_id, trun_num]  = FN_mat[case_sim,result_name_colname]
    AUC_res_array[case_sim, method_id, trun_num] = AUC_mat[case_sim,result_name_colname]

  }
  
  
  # # =========== 2. Seurat methods
  # 
  train_for_DE = train
  seurat_object = CreateSeuratObject(counts = t(train_for_DE[,-1]))
  my_labels = as.character(train$y)
  Idents(object = seurat_object) =  as.factor(my_labels)

  for(DE_method in DE_method_list){
    # asdasd
    method_id = method_id + 1

    marker_result = FindMarkers(seurat_object, ident.1 = "1", ident.2 = "0",
                                test.use = DE_method) # ("wilcox" is default)

    for(trun_num in 1:(p-1)){
      num_vars_to_use = min(nrow(marker_result), trun_num)
      best_model_vars0 = rownames(marker_result)[1:num_vars_to_use]
      best_model_vars = colnames(train_for_DE[,-1])[colnames(train_for_DE[,-1]) %in% best_model_vars0]

      result_name_colname = paste(best_model_vars, collapse = "...")

      FP_mat  = FP_mat_list[[num_vars_to_use]]
      FN_mat  = FN_mat_list[[num_vars_to_use]]
      AUC_mat = AUC_mat_list[[num_vars_to_use]]
      
      FP_res_array[case_sim, method_id, trun_num]  = FP_mat[case_sim,result_name_colname]
      FN_res_array[case_sim, method_id, trun_num]  = FN_mat[case_sim,result_name_colname]
      AUC_res_array[case_sim, method_id, trun_num] = AUC_mat[case_sim,result_name_colname]
    } # case_trun_num

  } # DE_method
  
  # ========= Method 3: SCANPY methods =========
  for(SCAN_method in SCAN_method_list){
    method_id = method_id + 1
    
    SCAN_file = paste0("../scanpy_results/real2_", RNA_name, "_sim", case_sim,
                       "_", SCAN_method, ".csv")
    scan_res = read.csv(SCAN_file, header = T, skip = 1)
    for(trun_num in 1:(p-1)){
      
      
      num_vars_to_use =  trun_num
      best_model_vars0 = scan_res[1:num_vars_to_use,1]
      best_model_vars = colnames(train[,-1])[colnames(train[,-1]) %in% best_model_vars0]
      
      result_name_colname = paste(best_model_vars, collapse = "...")
      
      FP_mat  = FP_mat_list[[num_vars_to_use]]
      FN_mat  = FN_mat_list[[num_vars_to_use]]
      AUC_mat = AUC_mat_list[[num_vars_to_use]]
      
      FP_res_array[case_sim, method_id, trun_num]  = FP_mat[case_sim,result_name_colname]
      FN_res_array[case_sim, method_id, trun_num]  = FN_mat[case_sim,result_name_colname]
      AUC_res_array[case_sim, method_id, trun_num] = AUC_mat[case_sim,result_name_colname]
    } # trun_num
    
  } # SCAN_method
  
} # case_sim


save(file = paste0("../results/Jun8/AUC_res_array ", RNA_name,".rda"), 
     AUC_res_array, FP_res_array, FN_res_array)

}

stopCluster(cl)
 
par(mfrow = c(1, 1)) 

for(RNA_name in RNA_namelist){


# === Reading benchmark methods results
  # load(paste0("all_objects_real2_", RNA_name, ".rda"))


# === Reading quantum results

quantum_res = read.csv(paste0("../data/result_BGS/results_",
                              RNA_name, "/BGS_select_auc.csv"))
quantum_res = quantum_res[,-1]
dim(quantum_res)

# current_array = (MSE_res_array)
# current_array = (FN_res_array)
current_array = (AUC_res_array)
dim(current_array)


  random_per_sd = 0
  lwd = 1
  ploting_methods = c( "wilcox", "bimod", "roc", "t", "LR", "MAST")
  ploting_methods = c('BSS','BGS', paste0("Seurat_", ploting_methods))


  ploting_methods = method_names


  col_set = rainbow(length(ploting_methods))
  mean_array <- apply(current_array, c(2, 3), mean)
  sd_array <- apply(current_array, c(2, 3), sd)
  rownames(mean_array) = method_names
  rownames(sd_array) = method_names

  x_range = 1:(p-1)
  col_id = 1
  par(mar = c(5, 5, 5, 5))

  yrange = c(0.9,1)

  ## plot BSS
  random_per = rnorm(1,0,random_per_sd)
  plot(x_range + random_per, mean_array[1,x_range], type = "b",
       col = col_set[col_id],
       main = RNA_name,
       pch = col_id,
       ylim = yrange,
       xlim = range(x_range),
       lwd = lwd,
       ylab = "AUC", xlab = "Selected number of variables",
       cex.lab = 1.5,
       cex.axis = 1.5)

  ## plot others
  for(case_method in ploting_methods[-(1)]){
    col_id = col_id + 1
    random_per = rnorm(1,0,random_per_sd)
    lines(x_range + random_per, mean_array[case_method,x_range],
          pch = col_id,
          col = col_set[col_id],
          lwd = lwd, type = "b")
  }


  legend('bottomright',legend = ploting_methods,
         col=col_set, lty = rep(1,length(ploting_methods)),
         pch = 1:length(ploting_methods),
         cex = 1, lwd = lwd)
}

print(Sys.time() - t1)


