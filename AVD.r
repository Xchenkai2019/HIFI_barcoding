# 加载所需的包
library(igraph)
cpf.otu<-otu.relative[grepl("P__Cyanobacteria",otu.relative$Taxonomy) |grepl("P__Planctomycetota",otu.relative$Taxonomy)  |grepl("P__Chytridiomycota",otu.relative$Taxonomy) |grepl("P__Aphelidea",otu.relative$Taxonomy),]

cor_matrix <- cor(t(cpvf.otu[,2:80]), use = "pairwise.complete.obs") # 计算相关系数矩阵
adj_matrix <- ifelse(abs(cor_matrix) > 0.5, 1, 0) # 基于相关系数构建邻接矩阵，阈值设定为0.6
# 转换为 igraph 图对象
network <- graph.adjacency(adj_matrix, mode = "undirected")

# 移除不同比例的物种，并计算稳定性指标
remove_species <- function(network, proportion) {
  num_nodes <- vcount(network) # 物种数量
  num_removed <- round(proportion * num_nodes) # 需要移除的物种数量
  removed_nodes <- sample(1:num_nodes, num_removed, replace = FALSE) # 随机选择需要移除的物种
  reduced_network <- delete.vertices(network, removed_nodes) # 移除物种后的网络
  
  # 计算网络的平均度量和变异度量
  average_degree <- mean(degree(reduced_network)) # 平均度量
  degree_variability <- sd(degree(reduced_network)) # 度量变异度
  
  # 计算AVD指数
  AVD_index <- average_degree / degree_variability
  
  return(AVD_index)
}

# 计算不同比例移除物种后的稳定性指标
proportions <- seq(0, 0.9, by = 0.1) # 移除比例
num_iterations <- 999 # 重复次数
AVD_indices <- matrix(NA, nrow = length(proportions), ncol = num_iterations) # 存储结果

for (i in 1:num_iterations) {
  for (j in 1:length(proportions)) {
    AVD_indices[j, i] <- remove_species(network, proportions[j])
  }
}

# 显示稳定性指标的分布情况
colnames(AVD_indices) <- paste("Iteration", 1:num_iterations)
rownames(AVD_indices) <- paste("Proportion Removed", proportions)
cpf.AVD <- apply(ai,2,sum)/(1*nrow(cpf.otu[,2:80]))
cpf.ai <- abs(cpf.otu[,2:80]-apply(cpf.otu[,2:80], 1, mean))/apply(cpf.otu[,2:80], 1, sd)
print(AVD_indices)