#!/usr/local/bin/Rscript
library("ggvenn")
library("glue")
library("gridExtra")

venn_diagram = function(cluster){
  gene_df = data.frame(read.delim("~/Programming/khoi-modrek-lab/seurat-pipeline/rna_marker_concatenated.txt", sep = " "))
  two = gene_df[gene_df$DoseComparison == "2Gy_vs_0Gy" & gene_df$p_val < 0.05 & gene_df$Cluster %in% cluster, ]$gene
  six = gene_df[gene_df$DoseComparison == "6Gy_vs_0Gy" & gene_df$p_val < 0.05 & gene_df$Cluster %in% cluster, ]$gene
  if (length(cluster) == 4) cluster = "All"
    
  a = list("2Gy DE Genes" = two, "6Gy DE Genes" = six)
  custom_theme = theme(
                       plot.title = element_text(color="black", size=6, face="bold.italic", hjust = 0.5, vjust = 2.0))
  venn = ggvenn(a, c("2Gy DE Genes", "6Gy DE Genes"), set_name_size = 2, stroke_size = 0.5, text_size = 2) + 
      theme_void() + labs(title = glue("2Gy and 6Gy Comparison \n of shared DE genes\n Cluster {cluster} ")) + custom_theme +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))
  return(venn)

}
atac_venn_diagram = function(cluster){
  atac_df = data.frame(read.delim("data/all/atac_marker_concatenated_output.txt", sep = " "))
  two = atac_df[atac_df$DoseComparison == "2Gy_vs_0Gy" & atac_df$p_val < 0.05 & atac_df$Cluster %in% cluster, ]$marker
  six = atac_df[atac_df$DoseComparison == "6Gy_vs_0Gy" & atac_df$p_val < 0.05 & atac_df$Cluster %in% cluster, ]$marker
  if (length(cluster) == 4) cluster = "All"
    
  a = list("2Gy DE Genes" = two, "6Gy DE Genes" = six)
  custom_theme = theme(
                       plot.title = element_text(color="black", size=6, face="bold.italic", hjust = 0.5, vjust = 2.0))
  venn = ggvenn(a, c("2Gy DE Genes", "6Gy DE Genes"), set_name_size = 2, stroke_size = 0.5, text_size = 2) + 
      theme_void() + labs(title = glue("2Gy and 6Gy Comparison of \n Differentially Expressed ATAC Peaks \n Cluster {cluster}")) + custom_theme +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))
  return(venn)

}
if (sys.nframe() == 0){
  venn = venn_diagram(c(0,1,2,3))
  ggsave("output/venndiagrams.pdf", venn, device = "pdf")
  plot_list = list()
  atac_plot_list = list()
  for (i in 0:3){
    plot_list[[i+1]] = venn_diagram(c(i))
    atac_plot_list[[i+1]] = atac_venn_diagram(c(i))
  }
  plot_list[[5]] = venn_diagram(c(0,1,2,3))
  atac_plot_list[[5]] = atac_venn_diagram(c(0,1,2,3))
  pdf("output/venndiagrams.pdf")
  grid.arrange(grobs = plot_list, nrow = 2)
  grid.arrange(grobs = atac_plot_list, nrow = 2)
  dev.off()
}
