#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/pathwaygeneht.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/geneht.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/venndiagrams.R", echo = TRUE)

if (sys.nframe() == 0){
  ###Preprocessing and Setup ----------------------------------------------------------  
  ### Load in Data
  all_df <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/all/combinedpathway_withmore_mainfunction.txt", sep = "\t"))
  unique_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/unique/unique_combineddepathway_with_mainfunction.txt", sep = "\t"))
  common_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/common/common_combineddepathway_with_mainfunction.txt", sep = "\t"))

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  unique_df = unique_df[unique_df$pval <0.05, ]
  common_df = common_df[common_df$pval<0.05, ]

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  unique_df = unique_df[!is.na(unique_df$NES) | !is.na(unique_df$pval), ]
  common_df = common_df[!is.na(common_df$NES) | !is.na(common_df$pval), ]

  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Immune System")

  ###Create Figure 2 ----------------------------------------------------------------------
  venn_plot = list()
  for (i in 0:3){
    venn_plot[[i+1]] = venn_diagram(c(i))
  }
  venn_plot[[5]] = venn_diagram(c(0,1,2,3)) 
  venn = (venn_plot[[1]] | venn_plot[[2]] | venn_plot[[3]] | venn_plot[[4]] | venn_plot[[5]])  
  plot_list = list()
  plot_list[[1]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, dose = "2Gy_vs_0Gy","Unique", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list[[2]] = grid.grabExpr(draw(make_gene_heatmap2(unique_df, "2Gy_vs_0Gy", "Unique", category_names, rownames_fs = 12, rowtitle_fs = 12)))
  plot_list[[3]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, dose = "6Gy_vs_0Gy","Unique", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list[[4]] = grid.grabExpr(draw(make_gene_heatmap2(unique_df, "6Gy_vs_0Gy", "Unique", category_names, rownames_fs = 12)))
  plot_list[[5]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, dose = "2Gy_vs_0Gy","Common", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list[[6]] = grid.grabExpr(draw(make_gene_heatmap2(common_df, "2Gy_vs_0Gy", "Common", category_names, rownames_fs = 12, rowtitle_fs = 12)))
  plot_list[[7]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, dose = "6Gy_vs_0Gy","Common", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list[[8]] = grid.grabExpr(draw(make_gene_heatmap2(common_df, "6Gy_vs_0Gy", "Common", category_names, rownames_fs = 12, rowtitle_fs = 12)))

  figure = venn /
    wrap_plots(plot_list[[1]], plot_list[[2]], ncol = 2) /
    wrap_plots(plot_list[[3]], plot_list[[4]], ncol = 2) /
    wrap_plots(plot_list[[5]], plot_list[[6]], ncol = 2) /
    wrap_plots(plot_list[[7]], plot_list[[8]], ncol = 2) 
  ggsave("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/main.pdf",device = "pdf", plot = figure, width = 8, height = 11.5, units = "in", dpi = 300)
  
  

  

  ###Create Figure 4 ----------------------------------------------------------------------
}
