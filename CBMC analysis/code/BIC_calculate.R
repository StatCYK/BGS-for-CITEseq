# Compute all BICs
# ------------------------------------------------------------------------------
rm(list = ls())

library(foreach)
library(doParallel)
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("./CBMC analysis/code")
set.seed(123)
# ----- Settings ---------------------------------------------------------------

num_sims = 100

# ----- Main -------------------------------------------------------------------
dat = read.csv("../data/data_8000.csv",row.names=1)


x_names = colnames(dat)[-1]
num_vars = length(x_names)
cl <- makeCluster(38)
registerDoParallel(cl)
for(case_num_vars in 1:length(x_names)){
  
  print(case_num_vars)
  
  all_subset_matrix = t(combn(x_names, case_num_vars))
  
  num_subsets = nrow(all_subset_matrix)
  
  set_fmlas = rep(NA, num_subsets)
  matrix_BIC       = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
  matrix_test_cor = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
  matrix_test_MSE  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
  matrix_test_spearman  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
  for(case_subset in 1: num_subsets){
    subset = all_subset_matrix[case_subset,]
    aaa<- foreach(case_sim = c(1: num_sims), .combine="rbind")%dopar%{
      set.seed(case_sim)
      n_rows <- nrow(dat)
      row_indices <- sample(1:n_rows)
      sep_index <- round(0.75 * n_rows)
      train <- dat[row_indices[1:sep_index], ]
      test <- dat[row_indices[(sep_index+1):n_rows], ]
      fmla = as.formula(paste("y ~", paste(subset, collapse = " + ")))
      fit = lm(fmla, data = train)
      y_prd = predict(fit, newdata = test)
      
      c(BIC(fit), mean((y_prd - test$y)^2),cor(y_prd , test$y),cor(y_prd , test$y,method = "spearman"))
    }
    matrix_BIC[,case_subset] = aaa[,1]
    matrix_test_MSE[,case_subset] = aaa[,2]
    matrix_test_cor[,case_subset] = aaa[,3]
    matrix_test_spearman[,case_subset] = aaa[,4]
    set_fmlas[case_subset] =  paste(subset, collapse = " + ")
    
  }
  colnames(matrix_BIC)       = set_fmlas
  colnames(matrix_test_cor)   = set_fmlas
  colnames(matrix_test_MSE)  = set_fmlas
  colnames(matrix_test_spearman)   = set_fmlas

  write.csv(matrix_BIC,     paste0("../results_new/","num",case_num_vars,"_BIC.csv"),       row.names = F)
  write.csv(matrix_test_cor, paste0("../results_new/","num",case_num_vars,"_test_cor.csv"),   row.names = F)
  write.csv(matrix_test_MSE, paste0("../results_new/","num",case_num_vars,"_test_MSE.csv"), row.names = F)
  write.csv(matrix_test_spearman, paste0("../results_new/","num",case_num_vars,"_test_spearman.csv"), row.names = F)
}


stopCluster(cl)
































