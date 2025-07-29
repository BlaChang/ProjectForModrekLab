#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(purrr)
library(stringr)

make_gene_heatmap2 <- function(df,
                               gene_df,
                               dose,
                               name,
                               category_names = NULL,
                               filter_df_by_dose = TRUE,
                               rownames_fs = 4,
                               rowtitle_fs = 14){
  col_fun = colorRamp2(c(-2,0,2), c("blue", "white", "red"))
  ### Get all Genes that participate in Hypoxia, Aptosis, DNA Repair, and Epigenetics
  if (is.null(category_names)){
    category_names = unique(df$main_function)
    print(category_names)
    print(unique(df$main_function))
  }
  categories = list()
  for (category in category_names){
    pathway_list = unique(df[df$main_function %in% category,]$pathway)
    categories[[category]] <- pathway_list
  }
  ht_list = NULL
  for (category in category_names){
    filtered_df = df[df$pathway %in% categories[[category]],]
    if (filter_df_by_dose){
      filtered_df = filtered_df[filtered_df$DoseComparison == dose, ]
    }
    gene_list = unique(strsplit(paste(filtered_df[,"leadingEdge"], collapse = ","), ","))[[1]]
    formatted_gene_list = map_chr(gene_list, str_trim)
    filtered_gene_df = gene_df[gene_df$gene %in% formatted_gene_list & gene_df$DoseComparison == dose & gene_df$p_val < .05,] 
    if (nrow(filtered_gene_df) == 0){
      print(glue("No genes were found in {category}"))
      next
    } 
    filtered_gene_df = dcast(filtered_gene_df, gene ~  Cluster, value.var = "avg_log2FC", fill = 0)
    rownames = filtered_gene_df[["gene"]]
    required_columns <- c("0", "1", "2", "3")
    for (col in required_columns) {
      if (!col %in% names(filtered_gene_df)) {
        filtered_gene_df[[col]] <- 0
      }
    }
    filtered_gene_df = filtered_gene_df %>% select(any_of(c("0", "1", "2","3")))
    mat = data.matrix(filtered_gene_df)
    rownames(mat) = rownames
    mat[is.na(mat)] <- 0
    #mat = mat[rowSums(mat) != 0 ,]

    if (category == "Extracellular Matrix Formation"){category  <- "ECM Formation"}
    # Rows to highlight
    uniqueRows <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/unique/unique_rna_marker_concatenated_output.txt", sep = " "))$gene
    commonRows <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/common/common_rna_marker_concatenated_output.txt", sep = " "))$gene
    # Set stylings for row names and make our selected rows unique
    unique_row_idx <- which(rownames(mat) %in% uniqueRows)
    common_row_idx <- which(rownames(mat) %in% commonRows)

    fontsizes <- rep(rownames_fs,  nrow(mat))
    fontcolors <- rep('black', nrow(mat))
    fontcolors[unique_row_idx] <- 'red'
    fontcolors[common_row_idx] <- 'blue'
    fontfaces <- rep('plain',nrow(mat))
    fontfaces[unique_row_idx] <- 'bold'
    fontfaces[common_row_idx] <- 'bold'

    # Create text annotation object for displaying row names
    rowAnno <- rowAnnotation(rows = anno_text(rownames(mat), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

    ht_list = Heatmap(mat, 
                      cluster_columns = FALSE,
                      show_heatmap_legend = FALSE,
                      row_names_gp = gpar(fontsize = rownames_fs),
                      row_names_max_width = max_text_width(
                                                           rownames(mat), 
                                                           gp = gpar(fontsize = 12)
                                                           ),
                      show_row_names = FALSE,
                      right_annotation = rowAnno,
                      row_title_gp = gpar(fontsize = rowtitle_fs),
                      row_title_rot = 0 ,
                      column_title = glue("{dose} Differential Expression for {name} Enriched pathways"),
                      row_title = glue({category}),
                      cluster_rows = TRUE,
                      border_gp = gpar(col = "black", lty = 2),
                      col = col_fun)  %v% ht_list 
  }
  return(ht_list)
}

if (sys.nframe() == 0){
  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Immune System")
  ### Load in Data
  all_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/all/combinedpathway_withmore_mainfunction.txt", sep = "\t"))
  unique_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/unique/unique_combineddepathway_with_mainfunction.txt", sep = "\t"))
  common_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/common/common_combineddepathway_with_mainfunction.txt", sep = "\t"))
  
  all_gene_df = data.frame(read.delim("data/all/rna_marker_concatenated_output.txt", sep = " "))
  unique_gene_df <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/unique/unique_rna_marker_concatenated_output.txt", sep = " "))
  common_gene_df <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/common/common_rna_marker_concatenated_output.txt", sep = " "))

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  unique_df = unique_df[unique_df$pval <0.05, ]
  common_df = common_df[common_df$pval<0.05, ]

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  unique_df = unique_df[!is.na(unique_df$NES) | !is.na(unique_df$pval), ]
  common_df = common_df[!is.na(common_df$NES) | !is.na(common_df$pval), ]

  pdf("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/genes_heatmap.pdf", height = 11, width = 15)
  draw(make_gene_heatmap2(all_df, all_gene_df, "2Gy_vs_0Gy", "all", category_names, rownames_fs = 4))
  draw(make_gene_heatmap2(all_df, all_gene_df, "6Gy_vs_0Gy", "all", category_names, rownames_fs = 3))
  draw(make_gene_heatmap2(unique_df, unique_gene_df, "2Gy_vs_0Gy", "unique", category_names, rownames_fs = 4))
  draw(make_gene_heatmap2(unique_df, unique_gene_df, "6Gy_vs_0Gy", "unique", category_names, rownames_fs = 3))

  draw(make_gene_heatmap2(common_df, common_gene_df, "2Gy_vs_0Gy", "common", category_names, rownames_fs = 4))
  draw(make_gene_heatmap2(common_df, common_gene_df, "6Gy_vs_0Gy", "common", category_names, rownames_fs = 3))
  dev.off()
}
