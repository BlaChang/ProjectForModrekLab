#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(purrr)
library(stringr)

create_highlighted_row_anno = function(mat, rownames_fs, red_rownames, blue_rownames){
  # Set stylings for row names and make our selected rows unique
  red_row_idx <- which(rownames(mat) %in% red_rownames)
  blue_row_idx <- which(rownames(mat) %in% blue_rownames)

  fontsizes <- rep(rownames_fs,  nrow(mat))
  fontcolors <- rep('black', nrow(mat))
  fontcolors[red_row_idx] <- 'red'
  fontcolors[blue_row_idx] <- 'blue'
  fontfaces <- rep('plain',nrow(mat))
  fontfaces[red_row_idx] <- 'bold'
  fontfaces[blue_row_idx] <- 'bold'

  # Create text annotation object for displaying row names
  rowAnno <- rowAnnotation(rows = anno_text(rownames(mat), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))
  return(rowAnno)
}
make_gene_heatmap_by_pathway <- function(df,
                                         gene_df,
                                         dose,
                                         name,
                                         category_names = NULL,
                                         filter_df_by_dose = TRUE,
                                         gpar = list(rownames_fs = 14, rowtitle_fs = 14, title_fs = 8, 
                                                     cell_width = unit(10, "mm"), cell_height = unit(4.5, "mm"), cluster_rows = TRUE)){
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
  i = 0
  for (category in category_names){
    filtered_df = df[df$pathway %in% categories[[category]],]

    #The only reason why this line exists is for the validation DF which only has one dose  
    if (filter_df_by_dose){
      filtered_df = filtered_df[filtered_df$DoseComparison == dose, ]
    }

    pathways = unique(filtered_df$common_name)
    for (pathway in pathways){
      gene_list = unique(strsplit(paste(filtered_df[filtered_df$common_name == pathway,"leadingEdge"], collapse = ","), ","))[[1]]
      formatted_gene_list = map_chr(gene_list, str_trim)
      filtered_gene_df = gene_df[gene_df$gene %in% formatted_gene_list & gene_df$DoseComparison == dose & gene_df$p_val < .05,] 
      if (nrow(filtered_gene_df) == 0){
        print(glue("No genes were found in {pathway}"))
        next
      }
      important_genes = list()
      for (cluster in c(0,1,2,3)){
        important_genes = unique(c(important_genes, pick_important_genes(filtered_gene_df)))
      }
      filtered_gene_df = filtered_gene_df[filtered_gene_df$gene %in% important_genes, ]
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
      uniqueRows <- data.frame(read.delim("data/unique/unique_rna_marker_concatenated_output.txt", sep = " "))$gene
      commonRows <- data.frame(read.delim("data/common/common_rna_marker_concatenated_output.txt", sep = " "))$gene

      # Set stylings for row names and make our selected rows unique
      rowAnno <- create_highlighted_row_anno(mat, gpar$rownames_fs, uniqueRows, commonRows) 

      i = i + 1
      ht_list = Heatmap(mat, 
                        name = glue("ht{i}"),
                        cluster_columns = FALSE,
                        show_heatmap_legend = FALSE,
                        width = ncol(mat)* gpar$cell_width, 
                        height = nrow(mat)* gpar$cell_height,
                        row_names_gp = gpar(fontsize = gpar$rownames_fs),
                        row_names_max_width = max_text_width(
                                                             rownames(mat), 
                                                             gp = gpar(fontsize = 12)
                                                             ),
                        show_row_names = FALSE,
                        right_annotation = rowAnno,
                        row_title_gp = gpar(fontsize = gpar$rowtitle_fs),
                        row_title_rot = 0 ,
                        column_title = glue("{dose} Differential Expression for {name} Enriched pathways"),
                        row_title = glue({pathway}),
                        cluster_rows = gpar$cluster_rows,
                        border_gp = gpar(col = "black", lty = 2),
                        col = col_fun)  %v% ht_list 
    }
  }

  ht_list = draw(ht_list)
  decorate_row_title(glue("ht1"), {
                       gp = gpar(col = "red")
  }, slice = 1)
  return(ht_list)
}

#Basically chooses the top two positively regulated genes, and the top two down regulated genes 
pick_important_genes = function(gene_df){
  #positively regulated genes
  positive_gene_df = gene_df[gene_df$avg_log2FC > 0, ]
  positive_gene_df = positive_gene_df[order(positive_gene_df$avg_log2FC, decreasing = TRUE), ]
  positive_genes = head(positive_gene_df, 2)$gene

  negative_gene_df = gene_df[gene_df$avg_log2FC < 0, ]
  negative_gene_df = negative_gene_df[order(negative_gene_df$avg_log2FC, decreasing = FALSE), ]
  negative_genes = head(negative_gene_df, 2)$gene

  return(c(positive_genes, negative_genes))
}

if (sys.nframe() == 0){
  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Innate Immune System", "Adaptive Immune System", "Cytokine Signalling Immune System")
  ### Load in Data
  all_df = data.frame(read.delim("data/all/combinedpathway_with_mainfunction.txt", sep = "\t"))
  unique_df = data.frame(read.delim("data/unique/unique_combinedpathway_with_mainfunction.txt", sep = "\t"))
  common_df = data.frame(read.delim("data/common/common_combinedpathway_with_mainfunction.txt", sep = "\t"))

  all_gene_df = data.frame(read.delim("data/all/rna_marker_concatenated_output.txt", sep = " "))
  unique_gene_df <- data.frame(read.delim("data/unique/unique_rna_marker_concatenated_output.txt", sep = " "))
  common_gene_df <- data.frame(read.delim("data/common/common_rna_marker_concatenated_output.txt", sep = " "))

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  unique_df = unique_df[unique_df$pval <0.05, ]
  common_df = common_df[common_df$pval<0.05, ]

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  unique_df = unique_df[!is.na(unique_df$NES) | !is.na(unique_df$pval), ]
  common_df = common_df[!is.na(common_df$NES) | !is.na(common_df$pval), ]

  pdf("output/genehtbypathway.pdf", height = 11, width = 15)
  draw(make_gene_heatmap_by_pathway(all_df, all_gene_df, "2Gy_vs_0Gy", "all", category_names))
  draw(make_gene_heatmap_by_pathway(all_df, all_gene_df, "6Gy_vs_0Gy", "all", category_names))

  draw(make_gene_heatmap_by_pathway(unique_df, unique_gene_df, "2Gy_vs_0Gy", "unique", category_names))
  draw(make_gene_heatmap_by_pathway(unique_df, unique_gene_df, "6Gy_vs_0Gy", "unique", category_names))

  draw(make_gene_heatmap_by_pathway(common_df, common_gene_df, "2Gy_vs_0Gy", "common", category_names))
  draw(make_gene_heatmap_by_pathway(common_df, common_gene_df, "6Gy_vs_0Gy", "common", category_names))
  print(list_components())
  dev.off()
}
