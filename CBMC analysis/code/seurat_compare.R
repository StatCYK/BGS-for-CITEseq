# Name: seurat_compare.R
# Function: Under current problem setting, comparing BSS with Seurat & scanpy
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
# BiocManager::install("MAST")
# BiocManager::install("DESeq2")

setwd("../CBMC analysis/code")
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



# ----- Settings ---------------------------------------------------------------

num_sims = 100

# ----- Main -------------------------------------------------------------------
dat = read.csv("../data/data_8000.csv",row.names=1)
data_8000_unnormalized = read.csv("../data/data_8000_unnormalized.csv")
head(dat)
p = ncol(dat) - 1

# === Main

DE_method_list = c("DESeq2", "wilcox", "bimod", "roc", "t", "poisson", "LR",
                   "MAST")

method_names = c('BSS', paste0("Seurat_", DE_method_list),"PCA")
num_methods = length(method_names)

num_selected_results = matrix(NA, nrow = num_sims, ncol = num_methods)
MSE_results   = matrix(NA, nrow = num_sims, ncol = num_methods)
R2_results    = matrix(NA, nrow = num_sims, ncol = num_methods)
adjR2_results = matrix(NA, nrow = num_sims, ncol = num_methods)
colnames(MSE_results) = colnames(R2_results) = colnames(adjR2_results) = 
  colnames(num_selected_results) = method_names

truncR2_res_array = array(NA, dim = c(num_sims,num_methods, p))
MSE_res_array     = array(NA, dim = c(num_sims,num_methods, p))

for(case_sim in 1: num_sims){
  print(paste0("case_sim: ", case_sim," / ",num_sims))
  set.seed(case_sim)
  n_rows = nrow(dat)
  row_indices = sample(1:n_rows)
  sep_index = round(0.75 * n_rows)
  train = dat[row_indices[1:sep_index], ]
  test = dat[row_indices[(sep_index+1):n_rows], ]
  
  train_unnormal = data_8000_unnormalized[row_indices[1:sep_index],]
  
  method_id = 0
  # ========= Method 1: Our BSS ========= 
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
    lowest_bic = which.min(best_summary$bic)
    best_model_vars0 = as.character(best_summary$which[lowest_bic, ])
    best_model_vars =colnames(train)[best_model_vars0 == TRUE][-1]
    
    fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
    fit = lm(fmla, data = train)
    y_hat = predict(fit,newdata = test)
    
    MSE_res_array[case_sim, method_id, trun_num] = mean((y_hat - test$y)^2)
    truncR2_res_array[case_sim, method_id, trun_num] = r_sq(y_hat, test$y)
  } # case_trun_num
  
  # ========= Method 2: Seurat_DE methods ========= 
  DE_method_list = c("DESeq2", "wilcox", "bimod", "roc", "t", "poisson", "LR",
                         "MAST")
  
  for(DE_method in DE_method_list){
    method_id = method_id + 1
    
    if(DE_method %in% c("poisson", "DESeq2")){
      train_for_DE = train_unnormal
    }else{
      train_for_DE = train
    }
    
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
    } # case_trun_num
    
  } # DE_method
  
  
  # ========= Method 3: PCA ========= 
  method_id = method_id + 1
  pca_fit <- prcomp(train[,-1], center = TRUE, scale. = TRUE)
  pca_train = predict(pca_fit, newdata = train[,-1])
  pca_test <- predict(pca_fit, newdata = test[,-1])
  
  marker_result = colnames(pca_test)
  
  pca_traindata = data.frame(cbind(train$y, pca_train))
  pca_testdata = data.frame(cbind(test$y, pca_test))
  colnames(pca_testdata)[1] = colnames(pca_traindata)[1] = 'y'
  
  for(trun_num in 1:p){
    num_vars_to_use = min(length(marker_result), trun_num)
    best_model_vars = marker_result[1:num_vars_to_use]
    
    fmla = as.formula(paste("y ~", paste(best_model_vars, collapse = " + ")))
    fit = lm(fmla, data = pca_traindata)
    y_hat = predict(fit, newdata = pca_testdata)
    
    MSE_res_array[case_sim, method_id, trun_num] = mean((y_hat - test$y)^2)
    truncR2_res_array[case_sim, method_id, trun_num] = r_sq(y_hat, test$y)
  } # case_trun_num
  
}# case_sim


### May 10:  R^2(1:13), line plot
### extract model selected variables

save.image(file = "all_objects_real1.rda")
load("all_objects_real1.rda")
# ===== Plotting results =====

ploting_methods = c("DESeq2", "wilcox", "bimod", "roc", "t", "LR", "MAST")
ploting_methods = c('BSS','BGS', paste0("Seurat_", ploting_methods))

round(colMeans(MSE_results[,ploting_methods]),4)
boxplot(MSE_results[,ploting_methods], 
        xlab = "Methods", ylab = "Testing Set MSE",
        col = rainbow(num_methods))

round(colMeans(num_selected_results[,ploting_methods]),4)
boxplot(num_selected_results[,ploting_methods], 
        xlab = "Methods", ylab = "Number of ADTs selected",
        col = rainbow(num_methods))

colMeans(R2_results)
boxplot(R2_results, main = "Comparison of Groups", 
        xlab = "Groups", ylab = "Test R-squared", col = rainbow(num_methods))

colMeans(adjR2_results)
boxplot(adjR2_results, main = "Comparison of Groups", 
        xlab = "Groups", ylab = "Test adjusted R-squared", col = rainbow(num_methods))


{
ploting_methods = c("DESeq2", "wilcox", "bimod", "roc", "t", "LR", "MAST")
ploting_methods = c('BSS','BGS', paste0("Seurat_", ploting_methods))
# === Reading quantum results
  
quantum_res = read.csv("BGS_select.csv")
quantum_res = quantum_res[,-1]
dim(quantum_res)

current_array = (MSE_res_array[,,2:12])
dim(current_array)

col_set = rainbow(length(ploting_methods) + 1)
mean_array <- apply(current_array, c(2, 3), mean)
sd_array <- apply(current_array, c(2, 3), sd)
rownames(mean_array) = method_names
rownames(sd_array) = method_names

x_range = 1:11
col_id = 1
par(mar = c(5, 5, 5, 5))

## plot BSS
random_per = rnorm(1,0,0.3)
plot(x_range + random_per + 1, mean_array[1,x_range], type = "b", 
     col = col_set[col_id], 
     pch = col_id,
     ylim = c(0.47, 0.56),
     xlim = range(x_range)+1,
     lwd = 3,
     ylab = "MSE", xlab = "Number of variables",
     cex.lab = 1.5,
     cex.axis = 1.5)
arrows(x_range+ random_per+1, mean_array[1,x_range] - sd_array[1,x_range],
       x_range+ random_per+1, mean_array[1,x_range] + sd_array[1,x_range],
       angle = 90, code = 3, length = 0.05, col = col_set[col_id])

## plot BGS
col_id = col_id + 1
random_per = rnorm(1,0,0.3)
BGS_mean = colMeans(quantum_res)
BGS_sd = apply(quantum_res, 2, sd)
lines(x_range + random_per + 1, BGS_mean, 
      pch = 1, #col_id,
      col = col_set[col_id], 
      lwd = 3, type = "b")
arrows(x_range+ random_per+1, BGS_mean[x_range] - BGS_sd[x_range],
       x_range+ random_per+1, BGS_mean[x_range] + BGS_sd[x_range],
       angle = 90, code = 3, length = 0.05, col = col_set[col_id])

## plot others
for(case_method in ploting_methods[-(1:2)]){
  col_id = col_id + 1
  random_per = rnorm(1,0,0.0)
  lines(x_range + random_per + 1, mean_array[case_method,x_range], 
        pch = 1, #col_id,
        col = col_set[col_id], 
        lwd = 3, type = "b")
  arrows(x_range + random_per, mean_array[case_method,x_range] - sd_array[case_method,x_range],
         x_range + random_per, mean_array[case_method,x_range] + sd_array[case_method,x_range],
         angle = 90, code = 3, length = 0.05, col = col_set[col_id])
}


legend('topright',legend = ploting_methods, 
       col=col_set, lty = rep(1,length(ploting_methods)),
       pch = 1:length(ploting_methods),
       cex = 1.5)
}




print(Sys.time() - t1)


