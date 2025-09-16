sample_info<-read.csv('sample_info.csv',row.name=1)
OTU<-read.csv('OTU.csv',row.names = 1)
taxonomy<-read.csv('Taxonomy.csv',row.name=1)
env_data<-read.csv('env.csv',row.name=1)
dataset <- microtable$new(sample_table = sample_info,otu_table=otu,tax_table = taxonomy)
dataset$tidy_dataset()
alltaxonenv <- trans_env$new(dataset = dataset, add_data = env_data)
alltaxonenv$cal_ordination(method = "RDA", taxa_level = "Taxon")
alltaxonenv$trans_ordination(show_taxa = 14, adjust_arrow_length = TRUE, 
                             max_perc_env = 1.5, max_perc_tax = 1.5, 
                             min_perc_env = 0.2, min_perc_tax = 0.2
                             )
alltaxonenv$plot_ordination(plot_color = "Season",taxa_text_italic = FALSE,taxa_arrow_color ="#B3B3B3",taxa_text_color = "#B3B3B3",color_values = c("#33A02C","#E31A1C","#FF7F00", "#1F38B4"))
alltaxonenv$cal_ordination_anova()
alltaxonenv$cal_ordination_envfit()