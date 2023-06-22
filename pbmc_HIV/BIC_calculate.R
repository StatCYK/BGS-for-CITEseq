library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library(pROC)

marker.names <- read.csv("./data/ADTmarkers.csv")
cl <- makeCluster(100)
registerDoParallel(cl)


for(cell.name in unique(marker.names$celltype.l3)){
  print(cell.name)
  dat = read.csv(paste0("./data/",cell.name,"_160000.csv"))
  if (sum(dat$y)>500){
  res_path = paste0("./results_",cell.name)
  if (!dir.exists(res_path)){
    dir.create(res_path)
  }else{
    print("dir exists")
  }
  # the 10 markers find by WNN paper for cell markers are 
  cell.marker <- marker.names$protein[which(marker.names$celltype.l3 == cell.name)]
  cell.marker2 <- c()
  for(i in 1:length(cell.marker)){
      tmp = gsub("`","",cell.marker[i])
      tmp = gsub("-",".",tmp)
      tmp = gsub("/",".",tmp)
      cell.marker2 <- c(cell.marker2, tmp)
  }
  cell.marker2 <- colnames(dat)[2:11]
  num_sims = 100
  for(case_num_vars in 1:length(cell.marker2)){
    print(case_num_vars)
    
    all_subset_matrix = t(combn(cell.marker2, case_num_vars))
    num_subsets = nrow(all_subset_matrix)
    
    set_fmlas = rep(NA, num_subsets)
    matrix_BIC       = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
    matrix_test_accu  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
    matrix_test_error1  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
    matrix_test_error2  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
    matrix_test_auc  = data.frame(matrix(NA, nrow = num_sims, ncol = num_subsets))
    for(case_subset in 1:num_subsets){
      print(case_subset)
      subset.select = all_subset_matrix[case_subset,]
      select.string = paste(subset.select, collapse = " + ")
      dat.subset = dat[c("y",subset.select)]
      fmla = as.formula(paste("y ~",select.string ))
      # downsample to balance the cluster size
      n_rows <- nrow(dat)
      aaa<- foreach(case_sim = c(1:num_sims), .combine="rbind",.packages="pROC") %dopar% {
        set.seed(123*case_sim)
        row_indices <- sample(1:n_rows)
        sep_index <- round(0.75 * n_rows)
        train <- dat.subset[row_indices[1:sep_index], ]
        train.weight <- (0.5*train$y*(1/sum(train$y == 1)-1/sum(train$y == 0))+0.5/sum(train$y == 0))
        train['weight'] = train.weight*nrow(train)
        test <- dat.subset[row_indices[(sep_index+1):n_rows],]
        
        fit = glm(fmla, family = binomial("logit"),data = train,weights = weight)
        y_prd = (predict(fit, newdata = test,type = "response")>0.5)+0
        
        c(BIC(fit),mean((y_prd - test$y)^2),sum((y_prd - test$y)^2*(test$y == 1))/sum(test$y == 1),
          sum((y_prd - test$y)^2*(test$y == 0))/sum(test$y == 0), auc(test$y,predict(fit, newdata = test,type = "response")))
      }
      matrix_BIC[,case_subset] = aaa[,1]
      matrix_test_accu[, case_subset] = aaa[,2]
      matrix_test_error1[, case_subset] = aaa[,3]
      matrix_test_error2[, case_subset] =aaa[,4]
      matrix_test_auc[,case_subset] =aaa[,5]
      
      set_fmlas[case_subset] =  paste(subset.select, collapse = " + ")
      
    }
    colnames(matrix_BIC)       = set_fmlas
    colnames(matrix_test_accu)   = set_fmlas
    colnames(matrix_test_error1)  = set_fmlas
    colnames(matrix_test_error2)  = set_fmlas
    colnames(matrix_test_auc)  = set_fmlas
    
    write.csv(matrix_BIC,     paste0(res_path,"/num",case_num_vars,"_BIC.csv"),       row.names = F)
    write.csv(matrix_test_accu, paste0(res_path,"/num",case_num_vars,"_test_accu.csv"),   row.names = F)
    write.csv(matrix_test_error1, paste0(res_path,"/num",case_num_vars,"_test_error1.csv"),   row.names = F)
    write.csv(matrix_test_error2, paste0(res_path,"/num",case_num_vars,"_test_error2.csv"), row.names = F)
    write.csv(matrix_test_auc, paste0(res_path,"/num",case_num_vars,"_test_auc.csv"), row.names = F)
  }
  }
}
stopCluster(cl)