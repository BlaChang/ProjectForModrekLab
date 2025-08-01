#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)


#Basically make Motif 
make_motif_heatmap = function(atac_df, direction, dose, name, rownames_fs, value.var = "count_peaks_with_motif", cluster_rows = TRUE){

  print(glue("{dose} {direction} {name}"))
  filtered_atac_df = atac_df[atac_df$pval < 0.05 & atac_df$status == direction & atac_df$DoseComparison == dose, ]

  casted_df = dcast(filtered_atac_df, motif_gene ~ Cluster, value.var = value.var, fun.aggregate = mean, fill = 0)

  rownames = casted_df$motif_gene
  required_columns = c("0", "1", "2", "3")
  for (col in required_columns) {
    if (!col %in% names(casted_df)) {
      casted_df[[col]] = 0
    }
  }
  mat = as.matrix(casted_df[, c("0", "1", "2", "3")])
  max_val = max(unlist(mat))
  print(max)
  print(head(mat))
  col_fun = colorRamp2(c(0,max_val), c("white", "red"))
  rownames(mat) = rownames
  colnames(mat) = required_columns
  mat[is.na(mat)] = 0


  ht = Heatmap(mat,
               cluster_columns = FALSE,
               cluster_rows = cluster_rows,
               show_heatmap_legend = TRUE,
               show_row_dend= FALSE,
               #width = unit(2, "cm"),
               #height = unit(10, "cm"),
               row_names_gp = gpar(fontsize = rownames_fs),
               row_names_max_width = max_text_width(
                                                    rownames(mat), 
                                                    gp = gpar(fontsize = rownames_fs)
                                                    ),
               column_title = glue("{name} {dose} {direction} Motif Analysis"),
               row_title_rot = 0 ,
               border_gp = gpar(col = "black", lty = 2),
               col = col_fun) 
  return(ht)


}

if (sys.nframe() == 0){
  all_atac_df = data.frame(read.delim("data/all/all_homer_distancezero_concat.txt", sep = " " )) 
  unique_atac_df = data.frame(read.delim("data/unique/unique_homer_distancezero.txt", sep = " " )) 
  common_atac_df = data.frame(read.delim("data/common/common_homer_distancezero.txt", sep = " " ))


  #Possibly later filter out NA rows
  pdf("output/atac_motif_heatmap.pdf", height = 11, width = 8.5)
  draw(make_motif_heatmap(unique_atac_df, "up", "2Gy_vs_0Gy", "Unique", rownames_fs = 12))
  draw(make_motif_heatmap(unique_atac_df, "up", "6Gy_vs_0Gy", "Unique", rownames_fs = 12))
  draw(make_motif_heatmap(unique_atac_df, "down", "2Gy_vs_0Gy", "Unique", rownames_fs = 12))
  draw(make_motif_heatmap(unique_atac_df, "down", "6Gy_vs_0Gy", "Unique", rownames_fs = 12))

  draw(make_motif_heatmap(common_atac_df, "up", "2Gy_vs_0Gy", "Common", rownames_fs = 12))
  draw(make_motif_heatmap(common_atac_df, "up", "6Gy_vs_0Gy", "Common", rownames_fs = 12))
  draw(make_motif_heatmap(common_atac_df, "down", "2Gy_vs_0Gy", "Common", rownames_fs = 12))
  draw(make_motif_heatmap(common_atac_df, "down", "6Gy_vs_0Gy", "Common", rownames_fs = 12))

  dev.off()


}

