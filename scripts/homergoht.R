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

make_homer_ht = function(name,
                         direction,
                         dose, 
                         category_names = NULL,
                         other = FALSE,
                         cluster_rows = TRUE,
                         gpar = list(rownames_fs = 14, rowtitle_fs = 14, title_fs = 8, cell_width = unit(10, "mm"), cell_height = unit(4.5, "mm"), cluster_rows = TRUE)){
  lowercase_name = tolower(name)
  lowercase_direction = tolower(direction)
  df = data.frame(read.delim(glue("data/{lowercase_name}/{lowercase_name}_homer_{lowercase_direction}_goanalysis_concat_with_mainfunction.txt"), sep = "\t"))

  ht = create_pathway_heatmap(df, dose, name = glue("{name} {direction}"), category_names = category_names, other = other, gpar = gpar)
  return(ht)

}
if (sys.nframe() == 0){
  ### Problem Statement is this, create a pathway heatmap with Cluster on the X-axis and Reactome GSEA pathways on the row-axis,
  ###but have this vary over direction(up,down), common, unique, and all, and dosage.  



  homer_gpar = list(rownames_fs = 6, rowtitle_fs = 6, title_fs = 8, cell_width = unit(10, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE)
  pdf("output/homergoht.pdf", width = 8.5, height = 11)  # Create PDF output file
  draw(make_homer_ht("All", "Up", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("All", "Down", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("All", "Up", "6Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("All", "Down", "6Gy_vs_0Gy", gpar = homer_gpar))

  draw(make_homer_ht("Common", "Up", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("Common", "Down", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("Common", "Up", "6Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("Common", "Down", "6Gy_vs_0Gy", gpar = homer_gpar))

  draw(make_homer_ht("unique", "up", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("unique", "down", "2Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("unique", "up", "6Gy_vs_0Gy", gpar = homer_gpar))
  draw(make_homer_ht("unique", "down", "6Gy_vs_0Gy", gpar = homer_gpar))
  dev.off()







}

