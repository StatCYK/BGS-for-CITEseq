# real1_boxplot.R
# =============================================================================
library(ggplot2)
DE_method_list = c("wilcox", "bimod", "roc", "t", "LR", "MAST")
method_names = c("BGS", paste0("Seurat-", DE_method_list))
SCAN_method_list = c('LR', 't', 'wilcox')
method_names = c(method_names, paste0("SCANPY-", SCAN_method_list))


quantum_res = read.csv("GBS_select 1.csv")
quantum_res = quantum_res[1:100,-1]


plotting_panel = c(2,4,6,8,10)

dim(MSE_res_array)
# 100  15  13
current_array = (MSE_res_array[,c(1,3:6,8:9,10:12),1:12])

current_array[,1,] = as.matrix(quantum_res)

my_array = current_array[,,plotting_panel]

dimnames(my_array) <- list(paste0("sim", 1:num_sims), 
                           method_names, 
                           paste0(1:12)[plotting_panel])
my_df <- as.data.frame.table(my_array, responseName = "value")
names(my_df) <- c("sim", "method", "parameter", "value")

color_palette <- c(c("red","#174D79", "#0000FF","#00CCFF","#A2CAED","#886EB2","#CC00FF"),
                   "#33A65CFF","#A2B627FF","#F8B620FF")
method_colors <- setNames(color_palette, unique(my_df$Method))

# box - all
ggplot(my_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot() +
  facet_wrap(~parameter, nrow = 1, strip.position = "bottom") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = fontsize + 3),      # Increase size of axis labels
        axis.text = element_text(size = fontsize), 
        strip.text.x = element_text(size = fontsize), 
        legend.text = element_text(size = fontsize),     # Increase size of legend text
        legend.title = element_text(size = fontsize))+
  scale_fill_manual(values = method_colors) +
  labs(title = "Boxplots for each method and parameter",
       x = "the number of selected ADTs",
       y = "MSE")
##########ratio
current_array = (MSE_res_array[,c(1,3:6,8:9,10:12),1:12])

current_array[,1,] = as.matrix(quantum_res)

all_wotBGS = current_array[,-1,]
dim(all_wotBGS)
ratio_array = array(NA, dim = dim(all_wotBGS))

for(case_method in 1:9){
  for(case_num in 1:12){
    for(case_sim in 1:num_sims){
      ratio_array[case_sim, case_method ,case_num] = 
        all_wotBGS[case_sim, case_method ,case_num] / as.matrix(quantum_res)[case_sim, case_num]
    }
  }
}
plotting_panel = c(2,4,6,8,10)
ratio_array = log(ratio_array[,,plotting_panel])


#Ratio plot
met_ids = 1:6


my_array = ratio_array[,met_ids,]
dim(my_array )

dimnames(my_array) <- list(paste0("sim", 1:num_sims), 
                           method_names[-1][met_ids], 
                           paste0(1:12)[plotting_panel])
my_df <- as.data.frame.table(my_array, responseName = "value")
names(my_df) <- c("sim", "method", "parameter", "value")

color_palette <- c(c("#174D79", "#0000FF","#00CCFF","#A2CAED","#886EB2","#CC00FF"),
                   "#33A65CFF","#A2B627FF","#F8B620FF")
method_colors <- setNames(color_palette, unique(my_df$Method))

# box - all
ggplot(my_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot() +
  facet_wrap(~parameter, nrow = 1, strip.position = "bottom") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = fontsize + 3),      # Increase size of axis labels
        axis.text = element_text(size = fontsize), 
        strip.text.x = element_text(size = fontsize), 
        legend.text = element_text(size = fontsize),     # Increase size of legend text
        legend.title = element_text(size = fontsize))+
  scale_fill_manual(values = method_colors[met_ids]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  labs(title = "",
       x = "the number of selected ADTs",
       y = "log MSE ratio")

#Ratio plot

met_ids = 7:9

my_array = ratio_array[,met_ids,]
dim(my_array )

dimnames(my_array) <- list(paste0("sim", 1:num_sims), 
                           method_names[-1][met_ids], 
                           paste0(1:12)[plotting_panel])
my_df <- as.data.frame.table(my_array, responseName = "value")
names(my_df) <- c("sim", "method", "parameter", "value")

color_palette <- c(c("#174D79", "#0000FF","#00CCFF","#A2CAED","#886EB2","#CC00FF"),
                   "#33A65CFF","#A2B627FF","#F8B620FF")
method_colors <- setNames(color_palette, unique(my_df$Method))

# box - all
ggplot(my_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot() +
  facet_wrap(~parameter, nrow = 1, strip.position = "bottom") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = fontsize + 3),      # Increase size of axis labels
        axis.text = element_text(size = fontsize), 
        strip.text.x = element_text(size = fontsize), 
        legend.text = element_text(size = fontsize),     # Increase size of legend text
        legend.title = element_text(size = fontsize))+
  scale_fill_manual(values = method_colors[met_ids]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  labs(title = "",
       x = "the number of selected ADTs",
       y = "log MSE ratio")













