#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)
library(readr)

source("scripts/pathwaygeneht.R", echo = TRUE)
source("scripts/geneht.R", echo = TRUE)
source("scripts/venndiagrams.R", echo = TRUE)
source("scripts/view_df.R", echo = TRUE)

df_to_heatmap <- function(df, 
                             dose,
                             name,
                             gpar = list(rownames_fs = 14, rowtitle_fs = 14, title_fs = 8, cell_width = unit(10, "mm"),
                                         cell_height = unit(4.5, "mm"), cluster_rows = TRUE)){
  #plots cluster vs Reactome pathways
  df = df[df$pval < 0.05, ]
  casted_df = dcast(df[df$DoseComparison == dose,], pathway ~  Cluster, value.var = "NES", fill = 0) 
  required_columns = c("0", "1", "2", "3")
  for (col in required_columns) {
    if (!col %in% names(casted_df)) {
      casted_df[[col]] = 0
    }
  }
  rownames = casted_df$pathway
  mat = data.matrix(casted_df[, c("0","1","2","3")])
  rownames(mat) = rownames
  colnames(mat) = c("0", "1", "2", "3")
  mat[is.na(mat)] <- 0
  col_fun = colorRamp2(c(-2,0,2), c("blue", "white", "red"))
  ht = Heatmap(mat,
               cluster_columns = FALSE,
               cluster_rows = gpar$cluster_rows,
               show_heatmap_legend = TRUE,
               width = ncol(mat)* gpar$cell_width, 
               height = nrow(mat)*gpar$cell_height, #1.5
               show_row_dend= FALSE,
               row_names_gp = gpar(fontsize = gpar$rownames_fs),
               row_names_max_width = max_text_width(
                                                    rownames(mat), 
                                                    gp = gpar(fontsize = gpar$rownames_fs)
                                                    ),
               show_row_names = TRUE,
               row_title_gp = gpar(fontsize = gpar$rowtitle_fs),
               column_title = glue("{name} {dose} Enriched Pathways"),
               column_title_gp = gpar(fontsize = gpar$title_fs),
               column_names_gp = gpar(fontsize = gpar$colnames_fs),
               heatmap_legend_param = list(title = "NES"),
               row_title_rot = 0,
               border_gp = gpar(col = "black", lty = 2),
               col = col_fun)
  return(ht)
}

make_homer_ht = function(name,
                         direction,
                         dose, 
                         gpar = list(rownames_fs = 14, rowtitle_fs = 14, title_fs = 8, cell_width = unit(10, "mm"), cell_height = unit(4.5, "mm"), cluster_rows = TRUE)){
  lowercase_name = tolower(name)
  lowercase_direction = tolower(direction)
  df = data.frame(read.delim(glue("data/{lowercase_name}/{lowercase_name}_gobp_homer_{lowercase_direction}_goanalysis.txt"), sep = "\t"))


  ht = df_to_heatmap(df, dose, name = glue("{name} {direction}"), gpar = gpar)
  return(ht)

}
if (sys.nframe() == 0){
  ### Problem Statement is this, create a pathway heatmap with Cluster on the X-axis and Reactome GSEA pathways on the row-axis,
  ###but have this vary over direction(up,down), common, unique, and all, and dosage.  



  homer_gpar = list(rownames_fs = 6, rowtitle_fs = 6, title_fs = 8, cell_width = unit(10, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE)
  pdf("output/homergoht.pdf", width = 11, height = 11)  # Create PDF output file
  draw(make_homer_ht("All", "Up", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  draw(make_homer_ht("All", "Down", "2Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  draw(make_homer_ht("All", "Up", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")
  draw(make_homer_ht("All", "Down", "6Gy_vs_0Gy", gpar = homer_gpar), heatmap_legend_side = "left")

  #draw(make_homer_ht("Common", "Up", "2Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("Common", "Down", "2Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("Common", "Up", "6Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("Common", "Down", "6Gy_vs_0Gy", gpar = homer_gpar))

  #draw(make_homer_ht("unique", "up", "2Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("unique", "down", "2Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("unique", "up", "6Gy_vs_0Gy", gpar = homer_gpar))
  #draw(make_homer_ht("unique", "down", "6Gy_vs_0Gy", gpar = homer_gpar))
  dev.off()







}

