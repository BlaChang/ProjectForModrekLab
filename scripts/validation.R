#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/pathwaygeneht.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/venndiagrams.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/venndiagrams.R", echo = TRUE)

change_column_name = function(df, old, new){
  colnames(df)[which(names(df) == old)] = new
  return(df)
}

if (sys.nframe() == 0){
  ###Preprocessing and Setup ----------------------------------------------------------  
  ### Load in Data
  all_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/all/combinedpathway_withmore_mainfunction.txt", sep = "\t"))
  hf2354_df = data.frame(read.delim("data/validation/hf2354_C2_CP-REACTOME_pathway_with_mainfunction.txt")) %>% mutate(Cluster =1) %>% mutate(DoseComparison = "Fractionated")
  hf3016_df = data.frame(read.delim("data/validation/hf3016_C2_CP-REACTOME_pathway_with_mainfunction.txt")) %>% mutate(Cluster =1) %>% mutate(DoseComparison = "Fractionated")

  all_gene_df = data.frame(read.delim("data/all/rna_marker_concatenated_output.txt", sep = " "))
  hf2354_gene_df = data.frame(read.delim("data/validation/hf2354_markers_same_direction.txt", sep = "\t"))
  hf3016_gene_df = data.frame(read.delim("data/validation/hf3016_markers_same_direction.txt", sep = "\t"))

  all_gene_path = "data/all/rna_marker_concatenated_output.txt"
  hf2354_gene_path = "data/validation/hf2354_markers_same_direction.txt"
  hf3016_gene_path = "data/validation/hf3016_markers_same_direction.txt"

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  hf2354_df = hf2354_df[hf2354_df$pval <0.05, ]
  hf3016_df = hf3016_df[hf3016_df$pval<0.05, ]

  

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  hf2354_df = hf2354_df[!is.na(hf2354_df$NES) | !is.na(hf2354_df$pval), ]
  hf3016_df = hf3016_df[!is.na(hf3016_df$NES) | !is.na(hf3016_df$pval), ]

  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Immune System")


  all_df_pathways = unique(all_df[all_df$Cluster == 1, ]$common_name)
  hf2354_df_pathways = unique(hf2354_df$common_name)
  hf3016_df_pathways = unique(hf3016_df$common_name)


  ###Get pathways that were enriched in both our data and hf2354 cell line validation data
  hf2354_common_pathways = intersect(all_df_pathways, hf2354_df_pathways)
  hf3016_common_pathways = intersect(all_df_pathways, hf3016_df_pathways)

  common_hf2354_df = hf2354_df[hf2354_df$common_name %in% hf2354_common_pathways, ]
  common_hf3016_df = hf3016_df[hf3016_df$common_name %in% hf3016_common_pathways, ]
  
  hf2354_ht = create_pathway_heatmap(common_hf2354_df, "hf2354", dose = "Fractionated", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)
  hf3016_ht = create_pathway_heatmap(common_hf3016_df, "hf3016", dose = "Fractionated", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)

  hf2354_gene_df = change_column_name(hf2354_gene_df, "avg_log2FC_2", "avg_log2FC")
  hf2354_gene_df = change_column_name(hf2354_gene_df, "p_val_2", "p_val")
  hf3016_gene_df = change_column_name(hf3016_gene_df, "avg_log2FC_2", "avg_log2FC")
  hf3016_gene_df = change_column_name(hf3016_gene_df, "p_val_2", "p_val")
  
  #hf2354_gene_ht = make_gene_heatmap2(common_hf2354_df, "2Gy_vs_0Gy", "hf2354", filter_df_by_dose = FALSE, gene_df = hf2354_gene_df, category_names, rownames_fs = 12, rowtitle_fs = 12)
  print(nrow(hf3016_df))
  hf3016_gene_ht = make_gene_heatmap2(common_hf3016_df, "2Gy_vs_0Gy", "hf3016", filter_df_by_dose = FALSE, gene_df = hf3016_gene_df, category_names = NULL, rownames_fs = 12, rowtitle_fs = 12)

#  hf2354xgsc20_gene_ht = make_gene_heatmap2(common_hf2354_df, "Fractionated", "h2354 compared with GSC20", category_names, gene_path = all_gene_path,rownames_fs = 12, rowtitle_fs = 12)
#  hf3016xgsc20_gene_ht = make_gene_heatmap2(common_hf3016_df, "Fractionated", "h3016 compared with GSC20", category_names, gene_path = all_gene_path,rownames_fs = 12, rowtitle_fs = 12)



  ###Draw Validation Figure ----------
  #hf2354_ht = grid.grabExpr(draw(hf2354_ht))
  #hf3016_ht = grid.grabExpr(draw(hf3016_ht))
#  hf2354_gene_ht = grid.grabExpr(draw(hf2354_gene_ht))
#  hf3016_gene_ht = grid.grabExpr(draw(hf3016_gene_ht))
#
#  hf2354xgsc20_gene_ht = grid.grabExpr(draw(hf2354xgsc20_gene_ht))
#  hf3016xgsc20_gene_ht = grid.grabExpr(draw(hf3016xgsc20_gene_ht))

#  figure = wrap_plots(hf2354_ht, hf2354_gene_ht, hf2354xgsc20_gene_ht, ncol = 3)
#  figure2 = wrap_plots(hf3016_ht, hf3016_gene_ht, hf3016xgsc20_gene_ht, ncol = 3)

  pdf("output/validation.pdf", height = 11, width = 8.5)
  draw(hf2354_ht)
  draw(hf3016_ht)
  #draw(hf2354_gene_ht)
  draw(hf3016_gene_ht)
  dev.off()

}

