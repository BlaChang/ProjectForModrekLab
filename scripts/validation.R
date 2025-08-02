#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)
library(purrr)
library(stringr)

change_column_name = function(df, old, new){
  colnames(df)[which(names(df) == old)] = new
  return(df)
}


source("scripts/pathwayht.R", echo = TRUE)
source("scripts/geneht.R", echo = TRUE)

if (sys.nframe() == 0){
  ###Preprocessing and Setup ----------------------------------------------------------  
  ### Load in Data
  all_df = data.frame(read.delim("data/all/combinedpathway_with_mainfunction.txt", sep = "\t"))
  hf2354_df = data.frame(read.delim("data/validation/hf2354_C2_CP-REACTOME_pathway_with_mainfunction.txt")) %>% mutate(Cluster =1) %>% mutate(DoseComparison = "Fractionated")
  hf3016_df = data.frame(read.delim("data/validation/hf3016_C2_CP-REACTOME_pathway_with_mainfunction.txt")) %>% mutate(Cluster =1) %>% mutate(DoseComparison = "Fractionated")

  all_gene_df = data.frame(read.delim("data/all/rna_marker_concatenated_output.txt", sep = " "))
  hf2354_gene_df = data.frame(read.delim("data/validation/hf2354_markers_same_direction.txt", sep = "\t"))
  hf3016_gene_df = data.frame(read.delim("data/validation/hf3016_markers_same_direction.txt", sep = "\t"))

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  hf2354_df = hf2354_df[hf2354_df$pval <0.05, ]
  hf3016_df = hf3016_df[hf3016_df$pval<0.05, ]

  all_gene_df = all_gene_df[all_gene_df$pval_2 < 0.05, ] %>% filter(!is.na(Cluster))
  hf2354_gene_df = hf2354_gene_df[hf2354_gene_df$p_val_2 <0.05, ]
  hf3016_gene_df = hf3016_gene_df[hf3016_gene_df$p_val_2<0.05, ]

  

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  hf2354_df = hf2354_df[!is.na(hf2354_df$NES) | !is.na(hf2354_df$pval), ]
  hf3016_df = hf3016_df[!is.na(hf3016_df$NES) | !is.na(hf3016_df$pval), ]

  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Innate Immune System", "Adaptive Immune System", "Cytokine Signalling Immune System")


  all_df_pathways = unique(all_df[all_df$Cluster == 1, ]$common_name)
  hf2354_df_pathways = unique(hf2354_df$common_name)
  hf3016_df_pathways = unique(hf3016_df$common_name)

  hf2354_common_pathways = intersect(all_df_pathways, hf2354_df_pathways)
  hf3016_common_pathways = intersect(all_df_pathways, hf3016_df_pathways)
  
  filtered_df = hf3016_df[hf3016_df$common_name %in% hf3016_common_pathways, ]
  gene_list = unique(strsplit(paste(filtered_df[,"leadingEdge"], collapse = ","), ","))[[1]]
  formatted_gene_list = map_chr(gene_list, str_trim)
  print(hf3016_gene_df$gene)
  c = intersect(formatted_gene_list, hf3016_gene_df$gene)

  print(c)
  common_hf2354_df = hf2354_df[hf2354_df$common_name %in% hf2354_common_pathways, ]
  common_hf3016_df = hf3016_df[hf3016_df$common_name %in% hf3016_common_pathways, ]
  
  
  validation_gpar = list(rownames_fs = 5, rowtitle_fs = 6, title_fs = 6, colnames_fs = 5, 
                  cell_width = unit(1.5, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE, show_heatmap_legend = FALSE)
  hf2354_ht = create_pathway_heatmap(common_hf2354_df, "hf2354", dose = "Fractionated", category_names = NULL, other = TRUE, gpar = validation_gpar)
  hf3016_ht = create_pathway_heatmap(common_hf3016_df, "hf3016", dose = "Fractionated", category_names = NULL, other = TRUE, gpar = validation_gpar)

  hf2354_gene_df = change_column_name(hf2354_gene_df, "avg_log2FC_2", "avg_log2FC")   ### Need to rename columns to be compatible with premade functions 
  hf2354_gene_df = change_column_name(hf2354_gene_df, "p_val_2", "p_val")
  hf3016_gene_df = change_column_name(hf3016_gene_df, "avg_log2FC_2", "avg_log2FC")
  hf3016_gene_df = change_column_name(hf3016_gene_df, "p_val_2", "p_val")
  
  
  val_gene_gpar = list(rownames_fs = 5, rowtitle_fs = 6, title_fs = 6, colnames_fs = 5, 
                  cell_width = unit(1.5, "mm"), cell_height = unit(1.5, "mm"), cluster_rows = TRUE, show_heatmap_legend = FALSE)

  hf2354_gene_df2 = hf2354_gene_df[hf2354_gene_df$DoseComparison == "2Gy_vs_0Gy", ]
  hf3016_gene_df2 = hf3016_gene_df[hf3016_gene_df$DoseComparison == "2Gy_vs_0Gy", ]


  hf2354_gene_ht2 = make_gene_heatmap2(hf2354_df, hf2354_gene_df2, "2Gy_vs_0Gy", "hf2354", filter_df_by_dose = FALSE, category_names = NULL, gpar = val_gene_gpar)
  hf3016_gene_ht2 = make_gene_heatmap2(hf3016_df, hf3016_gene_df2, "2Gy_vs_0Gy", "hf3016", filter_df_by_dose = FALSE, category_names = NULL, gpar = val_gene_gpar)

  hf2354_gene_df6 = hf2354_gene_df[hf2354_gene_df$DoseComparison == "6Gy_vs_0Gy", ]
  hf3016_gene_df6 = hf3016_gene_df[hf3016_gene_df$DoseComparison == "6Gy_vs_0Gy", ]
  hf2354_gene_ht6 = make_gene_heatmap2(hf2354_df, hf2354_gene_df6, "6Gy_vs_0Gy", "hf2354", filter_df_by_dose = FALSE, NULL, gpar = val_gene_gpar)
  hf3016_gene_ht6 = make_gene_heatmap2(hf3016_df, hf3016_gene_df6, "6Gy_vs_0Gy", "hf3016", filter_df_by_dose = FALSE, category_names = NULL, gpar = val_gene_gpar)

  pdf("output/validation.pdf", height = 11, width = 8.5)
  print(hf2354_ht)
  print(hf3016_ht)
  print(hf2354_gene_ht2)
  print(hf3016_gene_ht2)
  print(hf2354_gene_ht6)
  print(hf3016_gene_ht6)
  dev.off()

}

