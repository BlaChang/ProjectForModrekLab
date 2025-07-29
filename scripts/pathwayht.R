#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

create_pathway_heatmap = function(df,
                                   dose,
                                   name,
                                   category_names = NULL,
                                   other = FALSE,
                                   cluster_rows = TRUE,
                                   rownames_fs = 14,
                                   rowtitle_fs = 14) {
  if (class(df) == "list") {
    df = df[[1]]
    dose = dose[[1]]
    category_names = category_names[[1]]
    other = other[[1]]
  }
  df = df[, c("pathway", "Cluster", "NES", "main_function", "DoseComparison", "common_name")]
  casted_df = dcast(df[df$DoseComparison == dose,], common_name ~  Cluster, value.var = "NES", fill = 0) 
  required_columns = c("0", "1", "2", "3")
  for (col in required_columns) {
    if (!col %in% names(casted_df)) {
      casted_df[[col]] = 0
    }
  }
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  mat = as.matrix(casted_df[, c("0", "1", "2", "3")])
  rownames(mat) = casted_df$`rownames(df)`

  ###Basically get main_function in a more easily accessible data format.
  if (is.null(category_names)) {
    category_names = unique(df$main_function)
  }

  categories = list()
  for (category in category_names){
    pathway_list = unique(df[df$main_function %in% category, ]$common_name)
    categories[[category]] = pathway_list
  }


  ht_list = NULL
  for (category in category_names){
    if (category == "None" || category == "Differentiation" || category == "Autophagy" || (category == "Other" & !other)) next
    filtered_casted_df = casted_df[casted_df$common_name %in% categories[[category]],]
    if (dim(filtered_casted_df)[[1]] == 0) next

    rownames = filtered_casted_df$common_name
    mat = data.matrix(filtered_casted_df[, c("0","1","2","3")])
    rownames(mat) = rownames
    colnames(mat) = c("0", "1", "2", "3")
    mat[is.na(mat)] = 0

    if (category == "Extracellular Matrix Formation") category  = "ECM Formation"
    ht_list = Heatmap(mat,
                      cluster_columns = FALSE,
                      cluster_rows = cluster_rows,
                      show_heatmap_legend = FALSE,
                      show_row_dend= FALSE,
                      row_title = glue("{category}"),
                      row_names_gp = gpar(fontsize = rownames_fs),
                      row_names_max_width = max_text_width(
                                                           rownames(mat), 
                                                           gp = gpar(fontsize = rownames_fs)
                                                           ),
                      row_title_gp = gpar(fontsize = rowtitle_fs),
                      column_title = glue("{name} {dose} Enriched Pathways"),
                      row_title_rot = 0 ,
                      border_gp = gpar(col = "black", lty = 2),
                      col = col_fun) %v% ht_list
  }

  return(ht_list)
}


if (sys.nframe() == 0){
  ### Load in Data
  all_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/all/combinedpathway_withmore_mainfunction.txt", sep = "\t"))
  unique_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/unique/unique_combineddepathway_with_mainfunction.txt", sep = "\t"))
  common_df = data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/common/common_combineddepathway_with_mainfunction.txt", sep = "\t"))

  ### Filter By P value less that 0.5
  all_df = all_df[all_df$pval < 0.05, ] %>% filter(!is.na(Cluster))
  unique_df = unique_df[unique_df$pval <0.05, ]
  common_df = common_df[common_df$pval<0.05, ]

  #Remove NA values
  all_df = all_df[!is.na(all_df$NES) | !is.na(all_df$pval), ]
  unique_df = unique_df[!is.na(unique_df$NES) | !is.na(unique_df$pval), ]
  common_df = common_df[!is.na(common_df$NES) | !is.na(common_df$pval), ]

  category_names = c("Apoptosis", "DNA Repair", "Epigenetics", "Extracellular Matrix Formation", "Immune System")
  #category_names = unique(all_df$main_function)


  pdf("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/pathwayheatmap.pdf", width = 17, height = 11)

  draw(create_pathway_heatmap(all_df, dose = "2Gy_vs_0Gy","All", category_names, other = TRUE))
  draw(create_pathway_heatmap(all_df, dose = "6Gy_vs_0Gy", "All",category_names, other = TRUE))

  draw(create_pathway_heatmap(unique_df, dose = "2Gy_vs_0Gy", "Unique", category_names, other = TRUE))
  draw(create_pathway_heatmap(unique_df, dose = "6Gy_vs_0Gy", "Unique", category_names, other = TRUE))

  draw(create_pathway_heatmap(common_df, dose = "2Gy_vs_0Gy", "Common", category_names, other = TRUE))
  draw(create_pathway_heatmap(common_df, dose = "6Gy_vs_0Gy", "Common", category_names, other = TRUE))
  dev.off()

  a = common_df[common_df$DoseComparison == "2Gy_vs_0Gy" & common_df$main_function != "Other", "common_name"]
  b = common_df[common_df$DoseComparison == "6Gy_vs_0Gy" & common_df$main_function != "Other", "common_name"]
  c = intersect(a,b)
  plot_list = list()
  common_dose_common_df = common_df[common_df$common_name %in% c, ]
  plot_list[[1]] = grid.grabExpr(draw(create_pathway_heatmap(common_dose_common_df, dose = "2Gy_vs_0Gy", "Common", category_names, other = TRUE, cluster_rows = FALSE, rownames_fs = 8)))
  plot_list[[2]] = grid.grabExpr(draw(create_pathway_heatmap(common_dose_common_df, dose = "6Gy_vs_0Gy", "Common", category_names, other = TRUE, cluster_rows = FALSE, rownames_fs = 8)))

  # Save to PDF
  pdf("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/commonpathwayheatmap.pdf", width = 17, height = 11)  # Create PDF output file
  grid.arrange(grobs = plot_list, ncol = 2)  # Arrange in 4 columns
  dev.off()
}
