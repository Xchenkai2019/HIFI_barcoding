# ��������İ�
library(igraph)
cpf.otu<-otu.relative[grepl("P__Cyanobacteria",otu.relative$Taxonomy) |grepl("P__Planctomycetota",otu.relative$Taxonomy)  |grepl("P__Chytridiomycota",otu.relative$Taxonomy) |grepl("P__Aphelidea",otu.relative$Taxonomy),]

cor_matrix <- cor(t(cpvf.otu[,2:80]), use = "pairwise.complete.obs") # �������ϵ������
adj_matrix <- ifelse(abs(cor_matrix) > 0.5, 1, 0) # �������ϵ�������ڽӾ�����ֵ�趨Ϊ0.6
# ת��Ϊ igraph ͼ����
network <- graph.adjacency(adj_matrix, mode = "undirected")

# �Ƴ���ͬ���������֣��������ȶ���ָ��
remove_species <- function(network, proportion) {
  num_nodes <- vcount(network) # ��������
  num_removed <- round(proportion * num_nodes) # ��Ҫ�Ƴ�����������
  removed_nodes <- sample(1:num_nodes, num_removed, replace = FALSE) # ���ѡ����Ҫ�Ƴ�������
  reduced_network <- delete.vertices(network, removed_nodes) # �Ƴ����ֺ������
  
  # ���������ƽ�������ͱ������
  average_degree <- mean(degree(reduced_network)) # ƽ������
  degree_variability <- sd(degree(reduced_network)) # ���������
  
  # ����AVDָ��
  AVD_index <- average_degree / degree_variability
  
  return(AVD_index)
}

# ���㲻ͬ�����Ƴ����ֺ���ȶ���ָ��
proportions <- seq(0, 0.9, by = 0.1) # �Ƴ�����
num_iterations <- 999 # �ظ�����
AVD_indices <- matrix(NA, nrow = length(proportions), ncol = num_iterations) # �洢���

for (i in 1:num_iterations) {
  for (j in 1:length(proportions)) {
    AVD_indices[j, i] <- remove_species(network, proportions[j])
  }
}

# ��ʾ�ȶ���ָ��ķֲ����
colnames(AVD_indices) <- paste("Iteration", 1:num_iterations)
rownames(AVD_indices) <- paste("Proportion Removed", proportions)
cpf.AVD <- apply(ai,2,sum)/(1*nrow(cpf.otu[,2:80]))
cpf.ai <- abs(cpf.otu[,2:80]-apply(cpf.otu[,2:80], 1, mean))/apply(cpf.otu[,2:80], 1, sd)
print(AVD_indices)