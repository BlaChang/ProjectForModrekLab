# ProjectForModrekLab

##Structure of this project
.
├── README.md
├── Rplots.pdf
├── data
│   ├── all
│   │   └──. . . . .   # Data/TSV files that are not in this repo, but structured like this  
│   ├── common
│   │   └──. . . . .  
│   ├── unique
│   │   └──. . . . .  
│   └── validation
│       └──. . . . .
├── main.R            #Composes scripts in scripts/ to create Figure 2 and 4
├── scripts           #These are the "Components," which create individual graphs per page 
│   ├── atac.R        #Prints figure 4 but just one figure per page
│   ├── atac_motif_heatmap.R   #Motif Enrichment Heatmap generation
│   ├── geneht.R               #Gene Heatmap Functions and generation 
│   ├── genehtbypathway.R      #Gene Heatmap but divided by category 
│   ├── homergoht.R            #Gene Ontology pathway heatmap for ATAC peaks
│   ├── pathwayht.R            #General pathway heatmap function and generation
│   ├── template.R             #Just Template R file
│   ├── validation.R           #Validation heatmap generation
│   ├── venndiagrams.R         #Both Rnaseq and Atac Venn diagram per cluster generation     
│   └── view_df.R              #Useful Debugging Function(Dependency: Visidata)     
├── output #Notice the names match the script used to generate them.  
│   ├── atac.pdf
│   ├── atac_motif_heatmap.pdf
│   ├── geneht.pdf
│   ├── genehtbypathway.pdf
│   ├── homergoht.pdf
│   ├── main.pdf
│   ├── pathwayht.pdf
│   ├── validation.pdf
│   └── venndiagrams.pdf
├── reactome_addfunction_scripts #Python scripts for querying Reactome API and Msigdb to add a common_name and main_function column
│   └──. . . .                   #Eg. REACTOME_HATS_ACETYLATE_HISTONES  --> common_name: Hats Acetylate Histones, main_function: Epigenetics
      

## Example Usage

For all the scripts to work, your working directory needs to be in repo_folder/, not in scripts/.  

Example Scripts
```
Rscript main.R
open output/main.pdf
```

```
Rscript scripts/atac.R
open output/atac.pdef
```

##Important Notices and Modifying the Code  

###What is the difference between make_gene_heatmap(), make_gene_heatmap2(), make_gene_heatmap3(), and make_gene_heatmap_by_pathway()
- make_gene_heatmap() is just a basic heatmap of Differential Expression of genes (in the leading edge of some pathways)
- make_gene_heatmap2() is a ComplexHeatmap List with these genes separated into main functions 
- make_gene_heatmap_by_pathway() is a ComplexHeatmap List with these genes separated by pathway they contribute to
- make_gene_heatmap3() is like make_gene_heatmap2() but filtered out for the two most up-regulated and two most down-regulated genes per pathway. 








      
      
      
      
      
      
      


##To Do 
- Annotate The code
- Write a Good README
- Send Email to Khoi with the data ZIP file
