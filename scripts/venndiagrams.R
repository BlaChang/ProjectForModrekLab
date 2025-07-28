#!/usr/local/bin/Rscript
library("ggvenn")
library("glue")
library("gridExtra")

venn_diagram = function(cluster){
  gene_df = data.frame(read.delim("~/Programming/khoi-modrek-lab/seurat-pipeline/rna_marker_concatenated.txt", sep = " "))
  two = gene_df[gene_df$DoseComparison == "2Gy_vs_0Gy" & gene_df$p_val < 0.05 & gene_df$Cluster %in% cluster, ]$gene
  six = gene_df[gene_df$DoseComparison == "6Gy_vs_0Gy" & gene_df$p_val < 0.05 & gene_df$Cluster %in% cluster, ]$gene
  if (length(cluster) > 1) cluster = "All"
    
  a = list("2Gy DE Genes" = two, "6Gy DE Genes" = six)
  custom_theme = theme(
                       plot.title = element_text(color="red", size=8, face="bold.italic", hjust = 0.5),
  )
  venn = ggvenn(a, c("2Gy DE Genes", "6Gy DE Genes"), set_name_size = 4) + theme_void()+ labs(title = glue("Venn Diagram of 2Gy vs 6Gy for Cluster {cluster}")) + custom_theme
  return(venn)

}
if (sys.nframe() == 0){
  venn = venn_diagram(c(0,1,2,3))
  ggsave("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/venndiagrams.pdf", venn, device = "pdf")
  plot_list = list()
  for (i in 0:3){
    plot_list[[i+1]] = venn_diagram(c(i))
  }
  plot_list[[5]] = venn_diagram(c(0,1,2,3))
  pdf("/Users/blakechang/Programming/khoi-modrek-lab/figures/output/venndiagrams.pdf")
  grid.arrange(grobs = plot_list, nrow = 2)  + labs(title = "lol")# Arrange in 4 columns
  dev.off()
}
