#!/usr/local/bin/Rscript

source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/pathwayht.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/geneht.R", echo = TRUE)
source("/Users/blakechang/Programming/khoi-modrek-lab/figures/scripts/venndiagrams.R", echo = TRUE)

if (sys.nframe() == 0){

  ### Load in Data
  all_df <- data.frame(read.delim("/Users/blakechang/Programming/khoi-modrek-lab/figures/data/all/combinedpathway_withmore_mainfunction.txt", sep = "\t"))
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

  plot_list = list()
  plot_list[[1]] = grid.grabExpr(draw(create_pathway_heatmap(all_df, dose = "2Gy_vs_0Gy","All", category_names, other = TRUE, rownames_fs = 8)))
  plot_list[[2]] = grid.grabExpr(draw(make_gene_heatmap2(all_df, "2Gy_vs_0Gy", "all", category_names, rownames_fs = 4)))

  plot_list2 = list()
  plot_list2[[1]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, dose = "2Gy_vs_0Gy","unique", category_names, other = TRUE, rownames_fs = 8)))
  plot_list2[[2]] = grid.grabExpr(draw(make_gene_heatmap2(unique_df, "2Gy_vs_0Gy", "unique", category_names, rownames_fs = 12)))

  plot_list3 = list()
  plot_list3[[1]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, dose = "2Gy_vs_0Gy","common", category_names, other = TRUE, rownames_fs = 8)))
  plot_list3[[2]] = grid.grabExpr(draw(make_gene_heatmap2(common_df, "2Gy_vs_0Gy", "common", category_names, rownames_fs = 12)))
  # Save to PDF
  
  venn_plot = list()
  venn_plot2 = list()
  venn = venn_diagram(c(0,1,2,3))
  for (i in 0:3){
    venn_plot[[i+1]] = venn_diagram(c(i))
  }
  venn_plot2 = venn_diagram(c(0,1,2,3))
  plot_list4 = list()
  plot_list4[[1]] = grid.arrange(grobs = venn_plot, nrow = 1)
  plot_list4[[2]] = venn_plot2
  plot_list4[[3]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, dose = "2Gy_vs_0Gy","Unique", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list4[[4]] = grid.grabExpr(draw(make_gene_heatmap2(unique_df, "2Gy_vs_0Gy", "Unique", category_names, rownames_fs = 12, rowtitle_fs = 12)))
  plot_list4[[5]] = grid.grabExpr(draw(create_pathway_heatmap(unique_df, dose = "6Gy_vs_0Gy","Unique", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list4[[6]] = grid.grabExpr(draw(make_gene_heatmap2(unique_df, "6Gy_vs_0Gy", "Unique", category_names, rownames_fs = 12)))
  plot_list4[[7]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, dose = "2Gy_vs_0Gy","Common", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list4[[8]] = grid.grabExpr(draw(make_gene_heatmap2(common_df, "2Gy_vs_0Gy", "Common", category_names, rownames_fs = 12, rowtitle_fs = 12)))
  plot_list4[[9]] = grid.grabExpr(draw(create_pathway_heatmap(common_df, dose = "6Gy_vs_0Gy","Common", category_names, other = TRUE, rownames_fs = 5, rowtitle_fs = 6)))
  plot_list4[[10]] = grid.grabExpr(draw(make_gene_heatmap2(common_df, "6Gy_vs_0Gy", "Common", category_names, rownames_fs = 12, rowtitle_fs = 12)))



  pdf("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/pathway-gene-heatmap.pdf", width = 8.5, height = 11)  # Create PDF output file
  grid.arrange(grobs = plot_list, ncol = 2)  # Arrange in 4 columns
  grid.arrange(grobs = plot_list2, ncol = 2)  # Arrange in 4 columns
  grid.arrange(grobs = plot_list3, ncol = 2)  # Arrange in 4 columns
  grid.arrange(grobs = plot_list4, ncol = 2)
  dev.off()
}
