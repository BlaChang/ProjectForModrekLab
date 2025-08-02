#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

source("scripts/pathwaygeneht.R", echo = TRUE)
source("scripts/geneht.R", echo = TRUE)
source("scripts/venndiagrams.R", echo = TRUE)
source("scripts/atac_motif_heatmap.R", echo = TRUE)
source("scripts/homergoht.R", echo = TRUE)

if (sys.nframe() == 0){
  ###Preprocessing and Setup ----------------------------------------------------------  
  ### Load in Data
  all_df <- data.frame(read.delim("data/all/combinedpathway_with_mainfunction.txt", sep = "\t"))
  unique_df = data.frame(read.delim("data/unique/unique_combinedpathway_with_mainfunction.txt", sep = "\t"))
  common_df = data.frame(read.delim("data/common/common_combinedpathway_with_mainfunction.txt", sep = "\t"))

  all_gene_df = data.frame(read.delim("data/all/rna_marker_concatenated_output.txt", sep = " "))
  unique_gene_df = data.frame(read.delim("data/unique/unique_rna_marker_concatenated_output.txt", sep = " "))
  common_gene_df = data.frame(read.delim("data/common/common_rna_marker_concatenated_output.txt", sep = " "))

  all_gene_df = all_gene_df[all_gene_df$p_val < 0.05, ]

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  unique_df = unique_df[unique_df$pval <0.05, ]
  common_df = common_df[common_df$pval<0.05, ]

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  unique_df = unique_df[!is.na(unique_df$NES) | !is.na(unique_df$pval), ]
  common_df = common_df[!is.na(common_df$NES) | !is.na(common_df$pval), ]

  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Innate Immune System", "Adaptive Immune System", "Cytokine Signalling Immune System")

  ###Create Figure 2 ----------------------------------------------------------------------
  venn_plot = list()
  for (i in 0:3){
    venn_plot[[i+1]] = venn_diagram(c(i))
  }
  venn_plot[[5]] = venn_diagram(c(0,1,2,3)) 
  venn = (venn_plot[[1]] | venn_plot[[2]] | venn_plot[[3]] | venn_plot[[4]] | venn_plot[[5]])  
  pht_list = list()
  ght_list = list()


  pht_gpar = list(rownames_fs = 5, rowtitle_fs = 6, title_fs = 6, colnames_fs = 5, cell_width = unit(1.5, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE)
  pht_list[[1]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, "2Gy_vs_0Gy","Unique", category_names, other = TRUE, gpar = pht_gpar)))
  pht_list[[2]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, "6Gy_vs_0Gy","Unique", category_names, other = TRUE, gpar = pht_gpar)))
  pht_list[[3]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, "2Gy_vs_0Gy","Common", category_names, other = TRUE, gpar = pht_gpar)))
  pht_list[[4]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, "6Gy_vs_0Gy","Common", category_names, other = TRUE, gpar = pht_gpar)))

  ght_gpar = list(rownames_fs = 5, rowtitle_fs = 6, title_fs = 6, colnames_fs = 5, cell_width = unit(1.5, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE)
  ght_list[[1]] = grid.grabExpr(draw(make_gene_heatmap3(unique_df, all_gene_df, "2Gy_vs_0Gy", "Unique", category_names, gpar = ght_gpar)))
  ght_list[[2]] = grid.grabExpr(draw(make_gene_heatmap3(unique_df, all_gene_df, "6Gy_vs_0Gy", "Unique", category_names, gpar = ght_gpar)))
  ght_list[[3]] = grid.grabExpr(draw(make_gene_heatmap3(common_df, all_gene_df, "2Gy_vs_0Gy", "Common", category_names, gpar = ght_gpar)))
  ght_list[[4]] = grid.grabExpr(draw(make_gene_heatmap3(common_df, all_gene_df, "6Gy_vs_0Gy", "Common", category_names, gpar = ght_gpar)))

  figure1 = venn /
    wrap_plots(pht_list[[1]], ght_list[[1]], ncol = 2, widths = c(1.5,1)) /
    wrap_plots(pht_list[[2]], ght_list[[2]], ncol = 2, widths = c(1.5,1)) /
    wrap_plots(pht_list[[3]], ght_list[[3]], ncol = 2, widths = c(1.5,1)) /
    wrap_plots(pht_list[[4]], ght_list[[4]], ncol = 2, widths = c(1.5,1)) + plot_layout(heights = c(0.6,0.5,1,0.8,0.8)) 

  #ggsave("output/main.pdf",device = "pdf", plot = figure, width = 8, height = 11.5, units = "in", dpi = 300)





  ###Create Figure 4 --------------------------------------------------------------------
  atac_goht_list = list()
  unique_motifht_list = list()
  common_motifht_list = list()


  homer_gpar = list(rownames_fs = 4, rowtitle_fs = 6, colnames_fs = 4, title_fs = 6, cell_width = unit(1, "mm"), cell_height = unit(1, "mm"), cluster_rows = TRUE)
  atac_goht_list[[1]] = grid.grabExpr(draw(make_atac_go_heatmap("All", "Up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  atac_goht_list[[2]] = grid.grabExpr(draw(make_atac_go_heatmap("All", "Down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  atac_goht_list[[3]] = grid.grabExpr(draw(make_atac_go_heatmap("All", "Up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  atac_goht_list[[4]] = grid.grabExpr(draw(make_atac_go_heatmap("All", "Down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))

  unique_motifht_list[[1]] = grid.grabExpr(draw(make_motif_heatmap("Unique", "up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  unique_motifht_list[[2]] = grid.grabExpr(draw(make_motif_heatmap("Unique", "down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  unique_motifht_list[[3]] = grid.grabExpr(draw(make_motif_heatmap("Unique", "up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  unique_motifht_list[[4]] = grid.grabExpr(draw(make_motif_heatmap("Unique", "down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))

  common_motifht_list[[1]] = grid.grabExpr(draw(make_motif_heatmap("Common", "up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  common_motifht_list[[2]] = grid.grabExpr(draw(make_motif_heatmap("Common", "down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  common_motifht_list[[3]] = grid.grabExpr(draw(make_motif_heatmap("Common", "up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  common_motifht_list[[4]] = grid.grabExpr(draw(make_motif_heatmap("Common", "down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left"))
  



  figure2 = wrap_plots(atac_goht_list[[1]], unique_motifht_list[[1]], common_motifht_list[[1]], ncol = 3, widths = c(2,1,1)) /
    wrap_plots(atac_goht_list[[2]], unique_motifht_list[[2]], common_motifht_list[[2]], ncol = 3, widths = c(2,1,1)) /
    wrap_plots(atac_goht_list[[3]], unique_motifht_list[[3]], common_motifht_list[[3]], ncol = 3, widths = c(2,1,1)) /
    wrap_plots(atac_goht_list[[4]], unique_motifht_list[[4]], common_motifht_list[[4]], ncol = 3, widths = c(2,1,1)) + plot_layout(heights = c(0.4,1,1,0.8))

  pdf("output/main.pdf", width = 11, height = 8.5)
  print(figure1)
  print(figure2)  # This is very badly formatted, for script putting each plot on a separate page, check scripts/atac.R
  dev.off()


}
