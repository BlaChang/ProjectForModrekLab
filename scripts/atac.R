#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

source("scripts/geneht.R", echo = TRUE)
source("scripts/atac_motif_heatmap.R", echo = TRUE)
source("scripts/homergoht.R", echo = TRUE)
source("scripts/view_df.R", echo = TRUE)


if (sys.nframe() == 0){
  atac_goht_list = list()
  unique_motifht_list = list()
  common_motifht_list = list()


  homer_gpar = list(rownames_fs = 6, rowtitle_fs = 12, colnames_fs = 12, title_fs = 18, cell_width = unit(5, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE)
  atac_goht_list[[1]] = draw(make_atac_go_heatmap("All", "Up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  atac_goht_list[[2]] = draw(make_atac_go_heatmap("All", "Down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  atac_goht_list[[3]] = draw(make_atac_go_heatmap("All", "Up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  atac_goht_list[[4]] = draw(make_atac_go_heatmap("All", "Down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")

  homer_gpar = list(rownames_fs = 12, rowtitle_fs = 12, colnames_fs = 12, title_fs = 18, cell_width = unit(5, "mm"), cell_height = unit(5, "mm"), cluster_rows = TRUE)
  unique_motifht_list[[1]] = draw(make_motif_heatmap("Unique", "up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  unique_motifht_list[[2]] = draw(make_motif_heatmap("Unique", "down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  unique_motifht_list[[3]] = draw(make_motif_heatmap("Unique", "up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  unique_motifht_list[[4]] = draw(make_motif_heatmap("Unique", "down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")

  common_motifht_list[[1]] = draw(make_motif_heatmap("Common", "up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  common_motifht_list[[2]] = draw(make_motif_heatmap("Common", "down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  common_motifht_list[[3]] = draw(make_motif_heatmap("Common", "up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  common_motifht_list[[4]] = draw(make_motif_heatmap("Common", "down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")



  pdf("output/atac.pdf", width = 11, height = 8.5)
  print(atac_goht_list[[1]]) 
  print(unique_motifht_list[[1]])
  print(common_motifht_list[[1]])
  print(atac_goht_list[[2]]) 
  print(unique_motifht_list[[2]]) 
  print(common_motifht_list[[2]])
  print(atac_goht_list[[3]]) 
  print(unique_motifht_list[[3]]) 
  print(common_motifht_list[[3]])
  print(atac_goht_list[[4]]) 
  print(unique_motifht_list[[4]]) 
  print(common_motifht_list[[4]])
  dev.off()
}


