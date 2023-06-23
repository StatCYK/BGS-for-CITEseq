# Name: sim2_logistic_wtd.R
# Function: simulation for Weighted logistic regression using BSS 
# ==============================================================================
rm(list = ls())

library(leaps)
library(MASS)
library(dplyr)

library(foreach)
library(doParallel)


t0 = Sys.time()
# ==============================================================================
samp_size = 2000
p_set = 3:14
num_sim = 100

result_mat = matrix(NA, nrow = num_sim, ncol = length(p_set))


# cores <- 7 # replace this with the number of cores you want to use
# registerDoParallel(cores)

for(case_p in 1: length(p_set)){
  
  print(Sys.time() - t0)
  # === data generation 
  dimension = p_set[case_p]
  print(dimension)
  
  result_this_p <- matrix(0, nrow = num_sim, ncol = 1)
  
  result_this_p <- foreach(case_sim = 1:num_sim, .combine = 'cbind') %dopar% {
    # for(case_sim in 1:num_sim){

    set.seed(case_sim)
    
    cov_matrix <- matrix(0, nrow = dimension, ncol = dimension)
    
    for(i in 1:dimension){
      for(j in 1:dimension){cov_matrix[i, j] <- 0.7^abs(i - j)}
    }
    
    mean_vector <- rep(0, dimension)
    
    timess = 5
    
    X_mat <- mvrnorm(n = samp_size * timess, mu = mean_vector, Sigma = cov_matrix)
    dim(X_mat)
    
    X_mat = as.matrix(X_mat)
    colnames(X_mat) = paste0("X", 1:dimension)
    
    beta_vec = c(rep(1, floor(dimension/2)),
                 rep(0, ceiling(dimension/2)))
    beta_vec = as.matrix(beta_vec)
    
    p_vec = 1/(1 + exp(1)^(- X_mat %*%  beta_vec))
    
    Y_vec = rbinom(samp_size * timess,1,p_vec)
    
    
    index1 = which(Y_vec==1)
    index0 = which(Y_vec==0)
    
    index_fin = c(index1[1:250], index0[1:1750])
    
    X_mat_fin = X_mat[index_fin,]
    Y_vec_fin =  Y_vec[index_fin]
    
    dat = as.data.frame(X_mat_fin)
    dat$Y = Y_vec_fin
    
    ratio = sum(dat$Y)/samp_size
    print(ratio)
    
    # stop()
    
    temp <- (0.5*dat$Y*(1/sum(dat$Y== 1) - 1/sum(dat$Y == 0)) + 0.5/sum(dat$Y == 0))
    dat_weights = temp*length(dat$Y)

    # === BSS

    # Define predictors and response
    predictors <- paste0("X", 1:dimension)
    response <- 'Y'

    # Create an empty data frame to store results
    results <- data.frame(
      subset = character(0),
      BIC = numeric(0)
    )

    # Loop over all possible subsets
    BIC_set = c()
    subset_set = c()
    for (k in 1:length(predictors)) {

      subsets <- combn(predictors, k)
      for (i in 1:ncol(subsets)) {
        subset <- subsets[, i]
        formula <- as.formula(paste(response, paste(subset, collapse = '+'), sep = '~'))
        model <- glm(formula, data = dat,
                     # weights=dat_weights,
                     family = binomial("logit"))
        BIC <- BIC(model)

        results <- rbind(results, data.frame(subset = paste(subset, collapse = '+'),
                                             BIC = BIC))
      } # subsets
    } # dimensions

    write.csv(results,
              file = paste0("../simulation2_weighted/dim", dimension,"_sim",
                            case_sim, ".csv"))

    # Find the model with the smallest BIC
    best_model <- results[which.min(results$BIC), ]

    res = best_model$subset
    true_res = paste0(paste0("X", (1:dimension)[as.logical(beta_vec)]),
                      collapse = '+')

    ress = as.numeric(res == true_res)

    ress
  } #case_sim

  result_mat[, case_p] = result_this_p
  
} #case_p
colMeans(result_mat)
stopImplicitCluster()
print(Sys.time() - t0)














