#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)


#Basically make Motif 
#make_motif_heatmap = function(atac_df, direction, dose, name, rownames_fs, value.var = "count_peaks_with_motif", cluster_rows = TRUE){
#
#  print(glue("{dose} {direction} {name}"))
#  filtered_atac_df = atac_df[atac_df$pval < 0.05 & atac_df$status == direction & atac_df$DoseComparison == dose, ]
#
#  casted_df = dcast(filtered_atac_df, motif_gene ~ Cluster, value.var = value.var, fun.aggregate = mean, fill = 0)
#
#  rownames = casted_df$motif_gene
#  required_columns = c("0", "1", "2", "3")
#  for (col in required_columns) {
#    if (!col %in% names(casted_df)) {
#      casted_df[[col]] = 0
#    }
#  }
#  mat = as.matrix(casted_df[, c("0", "1", "2", "3")])
#  max_val = max(unlist(mat))
#  print(max)
#  print(head(mat))
#  col_fun = colorRamp2(c(0,max_val), c("white", "red"))
#  rownames(mat) = rownames
#  colnames(mat) = required_columns
#  mat[is.na(mat)] = 0
#
#
#  ht = Heatmap(mat,
#               cluster_columns = FALSE,
#               cluster_rows = cluster_rows,
#               show_heatmap_legend = TRUE,
#               show_row_dend= FALSE,
#               #width = unit(2, "cm"),
#               #height = unit(10, "cm"),
#               row_names_gp = gpar(fontsize = rownames_fs),
#               row_names_max_width = max_text_width(
#                                                    rownames(mat), 
#                                                    gp = gpar(fontsize = rownames_fs)
#                                                    ),
#               column_title = glue("{name} {dose} {direction} Motif Analysis"),
#               row_title_rot = 0 ,
#               border_gp = gpar(col = "black", lty = 2),
#               col = col_fun) 
#  return(ht)
#
#}


source("scripts/view_df.R", echo = TRUE)
source("scripts/geneht.R", echo = TRUE)

make_motif_heatmap = function(name,
                              direction,
                              dose, 
                              category_names = NULL,
                              other = FALSE,
                              gpar = list(rownames_fs = 14, rowtitle_fs = 14, title_fs = 8, cell_width = unit(10, "mm"), cell_height = unit(4.5, "mm"), cluster_rows = TRUE)){
  ###This function takes the pathways from the ALL go analysis, and then creates a heatmap of the genes of 

  lowercase_name = tolower(name)
  lowercase_direction = tolower(direction)

  homer_pathway_df = data.frame(read.delim(glue("data/all/all_gobp_homer_{lowercase_direction}_goanalysis.txt"), sep = "\t"))
  atac_df = data.frame(read.delim(glue("data/{lowercase_name}/{lowercase_name}_homer_distancezero.txt"), sep = " " ))  
  atac_df = atac_df[atac_df$status == lowercase_direction & atac_df$DoseComparison == dose, ]

  direction = glue("{direction} {name}")
  ht = make_gene_heatmap(homer_pathway_df, atac_df, dose, name = direction, value.var = "count_peaks_with_motif", gpar = gpar)
  return(ht)
}


if (sys.nframe() == 0){
  ### We are making the ATAC Motif heatmap based on the GSEA of the Motif genes
  ### Pathways from GSEA -> Leading Edge Genes -> plot these genes on a heatmap
  ### Options can very over unique/common motif genes, 2gy vs 0gy, and up/down of the ATAC(aka open/closed of the chromatin regions) 

  #Possibly later filter out NA rows

  #category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Innate Immune System", "Adaptive Immune System", "Cytokine Signalling Immune System")

  pdf("output/atac_motif_heatmap.pdf", height = 11, width = 8.5)
  draw(make_motif_heatmap("Unique", "up", "2Gy_vs_0Gy"))
  draw(make_motif_heatmap("Unique", "up", "6Gy_vs_0Gy"))
  draw(make_motif_heatmap("Unique", "down", "2Gy_vs_0Gy"))
  draw(make_motif_heatmap("Unique", "down", "6Gy_vs_0Gy"))

  draw(make_motif_heatmap("Common", "up", "2Gy_vs_0Gy"))
  draw(make_motif_heatmap("Common", "up", "6Gy_vs_0Gy"))
  draw(make_motif_heatmap("Common", "down", "2Gy_vs_0Gy"))
  draw(make_motif_heatmap("Common", "down", "6Gy_vs_0Gy"))
  dev.off()


}

