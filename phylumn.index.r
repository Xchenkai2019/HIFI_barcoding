# 计算每个门的总丰度值
phylum_abundance <- aggregate(Abundance ~ Phylum, data = df, sum)

# 计算每个门所属的界包含的门数量
phylum_kingdom_count <- aggregate(Phylum ~ Kingdom, data = df, FUN = function(x) length(unique(x)))

# 合并门的总丰度值和门所属的界包含的门数量
phylum_info <- merge(phylum_abundance, phylum_kingdom_count, by.x = "Phylum", by.y = "Kingdom")

# 计算门的相对丰度值
phylum_info$Relative_Abundance <- (phylum_info$Abundance * phylum_info$Phylum) / phylum_info$Kingdom