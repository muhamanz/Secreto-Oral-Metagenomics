#Analysis scripts for the manuscript: "Multi-kingdom oral microbiome interactions in early-onset cryptogenic ischemic stroke" by Muhammed Manzoor et al.
#The data in this study are not publicly available due to the sensitive nature of the data derived from human subjects, including personal information. However, researchers who wish to use our data may do so by providing summary statistics and analyses upon request from the SECRETO Oral Study Consortia Data Access Committee (https://ega-archive.org/dacs/EGAC00001003449).

if (!requireNamespace("BiocManager")) {
install.packages("BiocManager") 6	}
#remotes::install_github('microbiome/OMA')
#library(devtools)
#Bioconductor development version
#BiocManager::install("microbiome/mia", version="devel")
packages <- c("ggplot2", "biomformat", "ggthemes", "phyloseq", "vegan", "SpiecEasi", "patchwork", "microbiome", "tidyverse", "reshape2", "survival", "magrittr", "ggnewscale", "propr", "ComplexHeatmap", "maptree", "RColorBrewer", "rms", "viridis", "scales", "data.table")
#Phenotype data is loaded from the included R object
#Load metadata and biom file
biom <- import_biom("cuatroc.biom", parseFunction=parse_taxonomy_greengenes, parallel=TRUE)
metadata <- import_qiime_sample_data("metadata.txt")
#Construct the primary phyloseq object.
phyloseq <- merge_phyloseq(biom, metadata)
#ADD TREE information to the phyloseq object.
random_tree = rtree(ntaxa(phyloseq), rooted=TRUE, tip.label=taxa_names(phyloseq))
phyloseq <- merge_phyloseq(biom, metadata, random_tree)
#remove samples with less than 100000k reads (total)
