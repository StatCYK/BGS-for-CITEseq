# ==============================================================================

library(ggplot2)
library(paletteer)


# ==============================================================================

metric = "AUC"; ya_rg = c(-0.45,0.1)



DE_method_list = c("wilcox", "bimod", "roc", "t", "LR", "MAST")

method_names = c('BSS', "BGS", paste0("Seurat-", DE_method_list))

SCAN_method_list = c('LR', 't', 'wilcox')

the_method_names = c(method_names, paste0("SCANPY-", SCAN_method_list))


RNA_namelist = c("B intermediate kappa", "B intermediate lambda", "CD4 TEM_3", 
                 "B memory kappa", "B memory lambda", "B naive kappa", 
                 "B naive lambda", "CD4 CTL", "CD4 Naive", "CD4 TCM_1", 
                 "CD4 TCM_2", "CD4 TCM_3", "CD4 TEM_1", "CD8 Naive", "CD8 TCM_1",
                 "CD8 TCM_2", "CD8 TCM_3", "CD8 TEM_1", "CD8 TEM_2", "Treg Memory", 
                 "Treg Naive", "CD8 TEM_4", "CD8 TEM_5", "CD8 TEM_6", "CD14 Mono", 
                 "CD16 Mono", "cDC2_1", "cDC2_2", "gdT_1", "gdT_2", "gdT_3", 
                 "gdT_4", "MAIT", "NK Proliferating", "NK_1", "NK_2", "NK_3", 
                 "NK_4", "NK_CD56bright", "pDC", "Platelet")
RNA_namelist1 = c("B intermediate kappa")

pch_li = c(20,1:14,0,15:19)


par(mfrow = c(1,1))
plot_id = 0

all_result_array = array(dim = c(11, 9, 41))

new_array = array(NA, dim = c(3,9, 41))



case_RNA = 0
for(RNA_namee in RNA_namelist){
  case_RNA = case_RNA + 1
  
  load(paste0("../results/Jun14/AUC_res_array ", RNA_namee,".rda"))
  
  if(metric == "AUC"){current_array = (AUC_res_array)}
  if(metric == "specificity"){current_array = (1 - FP_res_array)}
  if(metric == "sensitivity"){current_array = (1 - FN_res_array)}
  #####################################################
  dim(current_array)
  
  random_per_sd = 0
  
  DE_method_list = c("wilcox", "bimod", "roc", "t", "LR", "MAST")
  
  method_names = c('BSS', paste0("Seurat_", DE_method_list))
  
  SCAN_method_list = c('logreg', 't-test', 'wilcoxon') #, 't-test_overestim_var')
  
  method_names = c(method_names, paste0("SCANPY_", SCAN_method_list))
  
  
  ploting_methods = method_names
  
  mean_array <- apply(current_array, c(2, 3), mean)
  
  rownames(mean_array) = method_names
  
  all_result_array[c(1,3:11), , case_RNA] =  mean_array
  
  
  
  #
  if(metric == "AUC"){quantum_res = read.csv(paste0("../real2newdata/select_result/results_",
                                                    RNA_namee, "/GBS_select_auc.csv"))}
  if(metric == "specificity"){quantum_res = 1 - read.csv(paste0("../real2newdata/select_result/results_",
                                                   RNA_namee, "/GBS_select_error1.csv"))}
  if(metric == "sensitivity"){quantum_res = 1 - read.csv(paste0("../real2newdata/select_result/results_",
                                                   RNA_namee, "/GBS_select_error2.csv"))}
  #####################################################
  quantum_res = quantum_res[,-1]
  BGS_mean = colMeans(quantum_res)
  
  all_result_array[2,,case_RNA] = BGS_mean
  #+#+#+#+##+

  
  new_array[1,1:9,case_RNA] =  BGS_mean
  for(case_num in 1:9){
    new_array[2, case_num, case_RNA] =  max(mean_array[2:7,case_num])
    new_array[3, case_num, case_RNA] =  max(mean_array[8:10,case_num])
  }
   

  
} # RNA_namee










dim(new_array)
# all

my_array = all_result_array[-1,,]




dimnames(my_array) <- list(the_method_names[-1], 
                           paste0(1:9),
                           paste0(1:41))
my_df <- as.data.frame.table(my_array, responseName = "value")
names(my_df) <- c( "Method", "panel_size","cell_type", "value")

fontsize = 20


all_wotBSS = all_result_array[-1,,]
dim(all_wotBSS)
ratio_array = array(NA, dim = c(9, 9, 41))

for(case_method in 1:9){
  for(case_num in 1:9){
    for(case_RNA in 1:41){
      ratio_array[case_method,case_num,case_RNA] = 
        all_wotBSS[case_method + 1,case_num,case_RNA] / all_wotBSS[1,case_num,case_RNA]
    }
  }
}
NUMrange = c(2,4,6,8)
ratio_array = log(ratio_array[,NUMrange,])

#====== final ratio plotting
{
# plotting_package = "Seurat"
# plotting_package = "SCANPY"
  plotting_package = "both"
if(plotting_package == "Seurat"){
  my_array = ratio_array[1:6,,]
  pkg_methods = the_method_names[-(1:2)][1:6]
  color_palette <- c("#174D79", "#0000FF","#00CCFF","#A2CAED","#886EB2","#CC00FF")
  y_rg = c(-0.6,0.1)
}
if(plotting_package == "SCANPY"){
  my_array = ratio_array[7:9,,]
  pkg_methods = the_method_names[-(1:2)][7:9]
  color_palette <- c("#33A65CFF","#A2B627FF","#F8B620FF")
}
if(plotting_package == "both"){
  my_array = ratio_array[1:9,,]
  pkg_methods = the_method_names[-(1:2)][1:9]
  color_palette <- c(c("#174D79", "#0000FF","#00CCFF","#A2CAED","#886EB2","#CC00FF"),
                     "#33A65CFF","#A2B627FF","#F8B620FF")
  }
dimnames(my_array) <- list(pkg_methods, 
                           paste0(NUMrange),
                           paste0(1:41))
my_df <- as.data.frame.table(my_array, responseName = "value")
names(my_df) <- c( "Method", "panel_size","cell_type", "value")
fontsize = 20
method_colors <- setNames(color_palette, unique(my_df$Method))
# box - all
ggplot(my_df, aes(x=Method, y=value, fill=Method)) +
  geom_boxplot() +
  # geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~panel_size, nrow = 1, strip.position = "bottom") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = fontsize + 3),      # Increase size of axis labels
        axis.text = element_text(size = fontsize), 
        strip.text.x = element_text(size = fontsize), 
        # legend.position = "none",
        legend.position = "bottom",
        legend.text = element_text(size = fontsize-2),     # Increase size of legend text
        legend.title = element_text(size = fontsize),
        legend.spacing.y = unit(0.25, "cm")) +
  # ylim(c(-0.02,0.02))+
  labs(title = "",
       x = "Panel size",
       y = paste0("log ", metric, " ratio"))+
  scale_fill_manual(values = method_colors,
                    guide = guide_legend(byrow = TRUE,direction = "horizontal"))  +
  coord_cartesian(ylim = y_rg)         
}



