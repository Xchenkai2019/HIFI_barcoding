# ����ÿ���ŵ��ܷ��ֵ
phylum_abundance <- aggregate(Abundance ~ Phylum, data = df, sum)

# ����ÿ���������Ľ������������
phylum_kingdom_count <- aggregate(Phylum ~ Kingdom, data = df, FUN = function(x) length(unique(x)))

# �ϲ��ŵ��ܷ��ֵ���������Ľ������������
phylum_info <- merge(phylum_abundance, phylum_kingdom_count, by.x = "Phylum", by.y = "Kingdom")

# �����ŵ���Է��ֵ
phylum_info$Relative_Abundance <- (phylum_info$Abundance * phylum_info$Phylum) / phylum_info$Kingdom