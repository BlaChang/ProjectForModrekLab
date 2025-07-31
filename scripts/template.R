#!/usr/local/bin/Rscript
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(glue)
library(gridExtra)
library(patchwork)

source("scripts/pathwaygeneht.R", echo = TRUE)
source("scripts/geneht.R", echo = TRUE)
source("scripts/venndiagrams.R", echo = TRUE)

if (sys.nframe() == 0){

}
