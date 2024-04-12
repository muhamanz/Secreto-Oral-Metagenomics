#Analysis scripts for the manuscript: "Multi-kingdom oral microbiome interactions in early-onset cryptogenic ischemic stroke" by Muhammed Manzoor et al.
#The data in this study are not publicly available due to the sensitive nature of the data derived from human subjects, including personal information. However, researchers who wish to use our data may do so by providing summary statistics and analyses upon request from the SECRETO Oral Study Consortia Data Access Committee (https://ega-archive.org/dacs/EGAC00001003449).
#Load packages
#if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
#remotes::install_github('microbiome/OMA')
#library(devtools)
#Bioconductor development version
#BiocManager::install("microbiome/mia", version="devel")
packages <- c("ggplot2", "biomformat", "ggthemes", "phyloseq", "vegan", "SpiecEasi", "patchwork", "microbiome", "tidyverse", " ggsignif", "reshape2", "survival", "magrittr", "ggnewscale", "propr", "ComplexHeatmap", "maptree", "RColorBrewer", "patchwork", "viridis", "miaViz", " file2meco ", "SpiecEasi", " chorddiag ", "scales", "data.table")
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
phyloseq <- subset_samples(phyloseq, sample_sums(phyloseq) > 100000)
# Trim out taxa that have less than 10 reads total:
phyloseq <- subset_taxa(phyloseq, taxa_sums(phyloseq) > 9)
#trim low-read samples
phyloseq <-rarefy_even_depth(phyloseq, sample.size = min(sample_sums(phyloseq)),
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#convert to tse object
tse<- makeTreeSummarizedExperimentFromPhyloseq(phyloseq)
#remove participants who have used antibiotics preceeding 3 months of saliva sampling
tse <- tse[ , colData(tse)$ONEMONTHANTIBIOTICS %in% FALSE]
#set a threshold for prevalance
tse <- subsetByPrevalentTaxa(tse, detection =  0, prevalence = 10/100, as_relative = TRUE) 
# add assay to tse object. To illustrate the use of multiple assays, the relative abundance data can be calculated and stored 
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
assays(tse)
# Add CLR transformation in the data as a new assay
tse <- transformSamples(tse, method = "clr", pseudocount=1)
# Community Diversity: Alpha diversity for observed species and shannon diversity
tse <- mia::estimateRichness(tse, 
                             assay.type = "counts", 
                             index = "observed", 
                             name="observed")
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "shannon", 
                              name = "shannon")
#Alpha diversities can be visualization for figure 1.
# creates data frame from the collected data
df <- as.data.frame(colData(tse))
# Changes old levels with new levels
df$Group  <- factor(df$Group)
ggplot(df, aes(x = Group, y = shannon )) + 
  geom_violin(trim=FALSE)+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1) + theme(text = element_text(size = 5))+ theme_classic() ggplot(df, aes(x = Group, y = observed )) + 
  geom_violin(trim=FALSE)+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1) + theme(text = element_text(size = 5))+ theme_classic()
#Abundance visualization jitter plot
tse <- transformSamples(tse, method = "relabundance")
plotAbundanceDensity(tse, layout = "jitter", assay_name = "relabundance",colour_by = "Group",
                     n = 10, point_size=4,  point_shape=19, point_alpha=0.5) + 
    scale_x_log10(label=scales::percent)
#Subset the feature by kingdom viruses, bacteria, Eukaryota, and archaea
tse_virus <- tse [rowData(tse)$Kingdom %in% c("Viruses"), ]
tse_Bacteria <- tse [rowData(tse)$Kingdom %in% c("Bacteria"), ]
tse_Eukaryota <- tse [rowData(tse)$Kingdom %in% c("Eukaryota"), ]
tse_Archaea <- tse [rowData(tse)$Kingdom %in% c("Archaea"), ]
#subset patients and controls
tse_patinet <- tse [ , colData(tse)$Group %in% c("patient"), ]
tse_control <- tse [ , colData(tse)$Group %in% c("control"), ]
#Agglomerate data at Phylum, Genus and species level
tse_Phylum<- agglomerateByRank(tse, rank="Phylum")
tse_Genus<- agglomerateByRank(tse, rank="Genus")
tse_Species<- agglomerateByRank(tse, rank="Species")
# Comparing communities by beta diversity analysis
# Agglomerate to genus level
tse_genus <- mergeFeaturesByRank(tse,
                                 rank = "Genus")
# Convert to relative abundances
tse_genus <- transformAssay(tse,
                            method = "relabundance",
                            assay.type = "counts")
# Add info on dominant genus per sample
tse_genus <- addPerSampleDominantFeatures(tse_genus,
                                          assay.type = "relabundance",
                                          name = "dominant_taxa")
#perform PCoA with Bray-Curtis dissimilarity.
tse_genus <- runMDS(tse_genus,
                    FUN = vegan::vegdist,
                    name = "PCoA_BC",
                    assay.type = "relabundance")
# Getting the top taxa
top_taxa <- getTopTaxa(tse_genus,top = 10, abund_values = "relabundance")
# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_genus)$dominant_taxa,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
# Storing the previous results as a new column within colData
colData(tse_genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq/sum(most_abundant_freq)*100, 1)
# Retrieving the tse_all_Genus variance
e <- attr(reducedDim(tse_genus, "PCoA"), "eig");
var_explained <- e/sum(e[e>0])*100
# Visualization
plot <-plotReducedDim(tse_genus,"PCoA", colour_by = "most_abundant", shape_by="Group", theme_size = 10, ncomponents = 2) +
  scale_colour_manual(values = c( "violet", "blue", "red","orange","brown","purple", "green", "magenta", "darkgreen", "yellow", "black"),     labels=paste0(names(most_abundant_percent),"(",most_abundant_percent,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       color="")
#Species level plot
tse_Species <- runMDS(tse_Species, FUN = vegan::vegdist,
              name = "PCoA", exprs_values = "relabundance")
# Getting the top taxa
top_taxa <- getTopTaxa(tse_Species,top = 10, abund_values = "relabundance")
# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_Species)$dominant_taxa,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
# Storing the tse_all_Species results as a new column within colData
colData(tse_Species)$most_abundant <- as.character(most_abundant)
# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq/sum(most_abundant_freq)*100, 1)
# Retrieving the tse_Species variance
e <- attr(reducedDim(tse_all_Species, "PCoA"), "eig");
var_explained <- e/sum(e[e>0])*100
# Visualization
plot <-plotReducedDim(tse_Species,"PCoA", colour_by = "most_abundant", shape_by="Group", theme_size = 10, ncomponents = 2) +
  scale_colour_manual(values = c( "violet", "blue", "red","orange","brown","purple", "green", "magenta", "darkgreen", "yellow", "black"),               labels=paste0(names(most_abundant_percent),"(",most_abundant_percent,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       color="")
# Plotting prevalence
altExps(tse) <- splitByRanks(tse)
altExps(tse) <-
   lapply(altExps(tse),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 0/100, sort = FALSE,
                                assay_name = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse,"Phylum"),
                        method="prevalence",
                        top=20L,
                        assay_name="counts")

top_phyla_mean <- getTopTaxa(altExp(tse,"Phylum"),
                             method="mean",
                             top=20L,
                             assay_name="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)
#visualization with plotted with plotRowTree
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
# plotRowTree for archaea
altExps(tse_Archaea) <- splitByRanks(tse_Archaea)
altExps(tse_Archaea) <-
   lapply(altExps(tse_Archaea),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 0/100, sort = FALSE,
                                assay_name = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse_Archaea,"Phylum"),
                        method="prevalence",
                        top=3L,
                        assay_name="counts")
top_phyla_mean <- getTopTaxa(altExp(tse_Archaea,"Phylum"),
                             method="mean",
                             top=3L,
                             assay_name="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse_Archaea)[1:6])
x <- addTaxonomyTree(x)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
#plotRowTree for bacteria
altExps(tse_Bacteria) <- splitByRanks(tse_Bacteria)
altExps(tse_Bacteria) <-
   lapply(altExps(tse_Bacteria),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 0/100, sort = FALSE,
                                assay_name = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse_Bacteria,"Phylum"),
                        method="prevalence",
                        top=35L,
                        assay_name="counts")
top_phyla_mean <- getTopTaxa(altExp(tse_Bacteria,"Phylum"),
                             method="mean",
                             top=35L,
                             assay_name="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
#plotRowTree for virus
altExps(tse_virus) <- splitByRanks(tse_virus)
altExps(tse_virus) <-
   lapply(altExps(tse_virus),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 0/100, sort = FALSE,
                                assay_name = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse_virus,"Phylum"),
                        method="prevalence",
                        top=3L,
                        assay_name="counts")
top_phyla_mean <- getTopTaxa(altExp(tse_virus,"Phylum"),
                             method="mean",
                             top=3L,
                             assay_name="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse_virus)[1:6])
x <- addTaxonomyTree(x)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
#plotRowTree for fungi
altExps(tse_Fungi) <- splitByRanks(tse_Fungi)
altExps(tse_Fungi) <-
   lapply(altExps(tse_Fungi),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 0/100, sort = FALSE,
                                assay_name = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse_Fungi,"Phylum"),
                        method="prevalence",
                        top=2L,
                        assay_name="counts")
top_phyla_mean <- getTopTaxa(altExp(tse_Fungi,"Phylum"),
                             method="mean",
                             top=2L,
                             assay_name="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse_Fungi)[1:6])
x <- addTaxonomyTree(x)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
#Visualizing taxonomic composition for bacteria, viruses, archaea, and fungi
se <- relAbundanceCounts(tse_Bacteria)
se_Genus <- agglomerateByRank(tse_Bacteria, rank ="Genus", onRankOnly=TRUE)
top_taxa <- getTopTaxa(se_Genus,top = 10, assay_name = "relabundance")
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(se)$Genus,
                       function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(se)$Genus <- as.character(genus_renamed)
plotAbundance(se, assay_name="relabundance", rank = "Genus",
           order_rank_by="abund",decreasing = TRUE,
  use_relative = TRUE,
  layout = c("bar", "point"),
  one_facet = TRUE)

se <- relAbundanceCounts(tse_Archaea)
se_Genus <- agglomerateByRank(tse_Archaea, rank ="Genus", onRankOnly=TRUE)
top_taxa <- getTopTaxa(se_Genus,top = 10, assay_name = "relabundance")
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(se)$Genus,
                       function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(se)$Genus <- as.character(genus_renamed)
plotAbundance(se, assay_name="relabundance", rank = "Genus",
           order_rank_by="abund",decreasing = TRUE,
  use_relative = TRUE,
  layout = c("bar", "point"),
  one_facet = TRUE)
se <- relAbundanceCounts(tse_Fungi)
se_Genus <- agglomerateByRank(tse_Fungi, rank ="Genus", onRankOnly=TRUE)
top_taxa <- getTopTaxa(se_Genus,top = 10, assay_name = "relabundance")
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(se)$Genus,
                       function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(se)$Genus <- as.character(genus_renamed)
plotAbundance(se, assay_name="relabundance", rank = "Genus",
           order_rank_by="abund",decreasing = TRUE,
  use_relative = TRUE,
  layout = c("bar", "point"),
  one_facet = TRUE)
se <- relAbundanceCounts(tse_virus)
se_Genus <- agglomerateByRank(tse_virus, rank ="Genus", onRankOnly=TRUE)
top_taxa <- getTopTaxa(se_Genus,top = 10, assay_name = "relabundance")
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(se)$Genus,
                       function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(se)$Genus <- as.character(genus_renamed)

plotAbundance(se, assay_name="relabundance", rank = "Genus",
           order_rank_by="abund",decreasing = TRUE,
  use_relative = TRUE,
  layout = c("bar", "point"),
  one_facet = TRUE)
#Composition heatmap
tse_Phylum <- transformFeatures(tse_Phylum, assay_name = "clr", 
                                       method = "z", name = "clr_z")
# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_Phylum, top = 30)
tse_phylum_subset <- tse_Phylum [top_taxa, ]
# Gets the assay table
mat <- assay(tse_Phylum, "clr_z")
# Creates the heatmap
pheatmap(mat)
# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")
# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)
# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins
# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)
# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 
# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), 
rownames(tse_phylum_subset))), ]
# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters
# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")
# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)
# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins
# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))
# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))
# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)
# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 
# Order data based on 
tse_phylum <- tse_phylum [ , rownames(sample_data)]
# Add sample type data
sample_data$sample_types <- factor(colData(tse_phylum)$Group)
sample_data
# Determines the scaling of colors
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)
pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors)
#Differential abundance analysis using The analysis of composition of microbiomes with bias correction (ANCOM-BC2) (Lin and Peddada 2020)
set.seed(123)
#unadjusted results and species level
ancombc2_out <- ancombc2(tse,
                         assay.type = "counts", tax_level = "Species",
                         fix_formula = "Group",
                         p_adj_method = "fdr",
                         lib_cut = 0,
                         group = "Group", 
                         struc_zero = TRUE, 
                         neg_lb = TRUE,
                         alpha = 0,
                         # multi-group comparison is deactivated automatically
                         global = TRUE)
#save the results
ancombc2_out$res %>%
  dplyr::select(starts_with(c("taxon", "lfc", "q"))) %>%
    head() %>%
  knitr::kable()
# ANCOM-BC2 with confounding
Output2 = ancombc2(data = tse, assay_name = "counts", tax_level = "Species",
                    fix_formula = "age+education+hypertension+smoking+caries+periodontitis+Group", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0, group = "Group", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE)
#save the results
Output2$res %>%
  dplyr::select(starts_with(c("taxon", "lfc", "q"))) %>%
    head() %>%
  knitr::kable()
#subset patients and controls
tse_patinet <- tse [ , colData(tse)$Group %in% c("patient"), ]
tse_control <- tse [ , colData(tse)$Group %in% c("control"), ]
#convert phyloseq into The microtable class
meco_all <- phyloseq2meco(phyloseq)
meco_patinet <- phyloseq2meco(phyloseq_patinet)
meco_control <- phyloseq2meco(phyloseq_control)
#correlation-based network using meco
# The parameter cor_method in trans_network is used to select correlation calculation method.
t1 <- trans_network$new(dataset = meco_control, cor_method = "sparcc",  cal_cor= "SparCC", use_sparcc_method = "SpiecEasi", filter_thres = 0.001) #1036
t2 <- trans_network$new(dataset = meco_patinet, cor_method = "sparcc",  cal_cor= "SparCC", use_sparcc_method = "SpiecEasi", filter_thres = 0.001)
# construct network; require igraph package
t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
t2$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
t2$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t1$res_network
# return t2$res_network
# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")
t2$cal_module(method = "cluster_fast_greedy")
# require rgexf package to be installed
t1$save_network(filepath = "network_control.gexf")
t2$save_network(filepath = "network_patient.gexf")
# calculate network attributes
t1$cal_network_attr()
t1$res_network_attr
t2$cal_network_attr()
t2$res_network_attr
# get node properties
t1$get_node_table(node_roles = TRUE)
t2$get_node_table(node_roles = TRUE)
# return t1$res_node_table # return t2$res_node_table
#Circos plots showed strong (r ≥ 0.8) correlations among multiple bacterial phyla, depicting the linkages between different phyla or within the same phylum
t1$cal_sum_links(taxa_level = "Phylum")
t2$cal_sum_links(taxa_level = "Phylum")
# interactive visualization; require chorddiag package; see https://github.com/mattflor/chorddiag
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
t2$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
t2$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
