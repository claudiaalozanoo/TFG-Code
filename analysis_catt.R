#############################################################################################################################
# Downstream analysis of integrated Seurat object
#############################################################################################################################

# clean up environment
rm(list = ls())
graphics.off()
gc()

# set working dir
setwd("/gpfs/scratch/bsc83/MN4/bsc83/bsc83962")

# Load deps
library(dplyr)
library(DESeq2)
library(Seurat)
library(stringr)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(gridExtra)

# source code
creatio.data.dir <- "/gpfs/projects/bsc83/MN4/bsc83/Data/Creatio/" # creatio.data.dir <- "~/marenostrum/data_creatio/"
creatio.proj.dir <- "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/" # creatio.proj.dir <- "~/marenostrum/projects_creatio/"
source("/gpfs/projects/bsc83/MN4/bsc83/Data/Creatio/code/src/seurat/seurat_src.R")

# integration version
#args <- commandArgs(trailingOnly = TRUE)
#int.ver <- args[1] # 
int.ver="v6"
#clust.res <- args[2] # 
clust.res="1"

# dirs
out.bse.dir <- str_interp("${creatio.proj.dir}integration_human/${int.ver}") # base directory
out.int.dir <- str_interp("${out.bse.dir}/integration") # integration directory
out.clust.dir <- str_interp("${out.bse.dir}/clustering") # surface protein analysis directory
out.surf.dir <- str_interp("${out.bse.dir}/surf_mark_analysis") # clustering directory
out.annot.dir <- str_interp("${out.bse.dir}/markers") # clustering directory
out.cc.dir <- str_interp("${out.bse.dir}/cell_cycle") # cell cycle analysis directory
out.postmit.dir <- str_interp("${out.bse.dir}/postmitotic") # postmitotic analysis directory
scvelo.dir <- str_interp("${out.bse.dir}/scvelo")

# subdirs
out.annot.mark.dir <- str_interp("${out.annot.dir}/wilcox_r${clust.res}/") # annotation marker directory

# create folder
dir.create(out.annot.mark.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.cc.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.postmit.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scvelo.dir, recursive = TRUE, showWarnings = FALSE)

# files
int.file.rds <- str_interp("${out.int.dir}/integration_human.rds") # no cortex or mge
int.file.cc.rds <- str_interp("${out.int.dir}/integration_human_cc_corr.rds")
int.file.new.rds <- str_interp("${out.int.dir}/integration_human_newclust_catt.rds")
mark.clust.file <- str_interp("${out.annot.mark.dir}/markergenes_${clust.res}.rds")
mark.clust.file.new <- str_interp("${out.annot.mark.dir}/markergenes_${clust.res}_new.rds")

# other files
go.file <- "/gpfs/projects/bsc83/MN4/bsc83/Data/Creatio/gsea/msigdb.v7.4.symbols.gmt"
bmrt.tab.file <- str_interp("${creatio.data.dir}/biomart/hsapiens_ensembl_table.tsv")
annot.tab.file <- str_interp("${out.clust.dir}/cluster_annot_${int.ver}_r${clust.res}.txt")

# allocate memory
options(future.globals.maxSize = 8000 * 1024^2)

# set seed for reproducibility
set.seed(1)

###########################################################################################################################
# Read in and assign cluster identities
###########################################################################################################################

# Load Seurat Object
so <- readRDS(int.file.rds)

# select resolution for active identity
Idents(so) <- str_interp("integrated_snn_res.${clust.res}")

# add cell_annot column if available
if (file.exists(annot.tab.file)) {
  # get cluster labels
  clust.labs <- read.table(annot.tab.file, header = TRUE, sep = "\t")

  # Add cell annotation labels
  so@meta.data$cell_annot <- plyr::mapvalues(
    x = so@meta.data[, str_interp("integrated_snn_res.${clust.res}")],
    from = clust.labs$cluster, to = clust.labs$label
  )
} else {
  # assign clustering with desired resolution
  so@meta.data$cell_annot <- so@active.ident

  # refactor
  so@meta.data$cell_annot <- factor(so@meta.data$cell_annot, paste(0:length(unique(so@meta.data$cell_annot))))
}

# fix day levels (missing samples)
pcw.levels <- levels(so$day)[table(so$day) > 0]
so@meta.data$day <- factor(so@meta.data$day, levels = pcw.levels)

# refactor div
so@meta.data$div <- factor(so@meta.data$div, levels = c("8", "12", "16"))
levels(so@meta.data$div) <- c("8 (DIV)", "12 (DIV)", "16 (DIV)")

# merge day and div, refactor, and create palette
so@meta.data$day <- factor(so@meta.data$day, levels = c(levels(so@meta.data$day), levels(so@meta.data$div)))
so@meta.data$day[so@meta.data$dataset == "CHDI"] <- so@meta.data$div[so@meta.data$dataset == "CHDI"]

# color palette for day
col.pal.pcw <- colorRampPalette(c("lightblue", "darkblue"))(length(pcw.levels))
col.pal.div <- colorRampPalette(c("lightgreen", "darkgreen"))(3)
col.pal.day <- c(col.pal.pcw, col.pal.div)

# check umap
DimPlot(so, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.clust.dir}/umap_seurat_clusters_r${clust.res}.pdf"), width = 20, height = 20)

# plot umap with new day colors
DimPlot(so, group.by = "day", label = TRUE, cols = col.pal.day) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = length(col.pal.day), byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.int.dir}/umap_day.png"), width = 20, height = 20)

###########################################################################################################################
# Data visualization Cattaneo
###########################################################################################################################
#subset cattaneo data 
so.cat <- subset(so, dataset =="Cattaneo")

# check umap
DimPlot(so.cat, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
#ggsave(str_interp("${out.clust.dir}/umap_seurat_clusters_r${clust.res}_cat.pdf"), width = 20, height = 20)

# plot umap with new day colors
DimPlot(so.cat, group.by = "day", label = FALSE, cols = col.pal.day) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = length(col.pal.day), byrow = TRUE, override.aes = list(size = 1)))
#ggsave(str_interp("${out.int.dir}/umap_day.png"), width = 20, height = 20)

# plot umap with phase
DimPlot(so.cat, group.by = "Phase", label = FALSE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = length(col.pal.day), byrow = TRUE, override.aes = list(size = 1)))
#ggsave(str_interp("${out.int.dir}/umap_day.png"), width = 20, height = 20)

#barplot distribution of cells in clusters
freq_table <- table(so.cat@meta.data$Phase, so.cat@meta.data$cell_annot)
colors <- c("G1" = "lightblue", "G2M" = "lightgreen", "S" = "lightpink")
barplot(freq_table, 
        main = "Distribution of Cells Along Clusters", 
        xlab = "Clusters", 
        ylab = "Cell number", 
        col = colors,
        ylim = c(0, 27000),
        legend.text = c("S", "G2M", "G1"),
        args.legend = list(x = "topright", bty = "n", fill = colors))
###########################################################################################################################
# Cell cycle per cluster entropy before correction
###########################################################################################################################

clust<- as.character(c(0:31))
entropy_pre_corr <- data.frame(Cluster = character(), Entropy = numeric())


# so.cat.clust0 <- subset(so.cat, cell_annot==clust)
# phase_table <- table(so.cat.clust0@meta.data$Phase)
# phase_probabilities <- as.vector(prop.table(phase_table)) #probabilities vector G1/G2M/S
# log_vector <- log(phase_probabilities, base=3)
# entropy <- sum(-phase_probabilities*log_vector)

for(i in clust){
  so.cat.clust <- subset(so.cat, cell_annot==i)
  phase_table <- table(so.cat.clust@meta.data$Phase)
  phase_probabilities <- as.vector(prop.table(phase_table)) #probabilities vector G1/G2M/S
  log_vector <- log(phase_probabilities, base=3)
  entropy <- sum(-phase_probabilities*log_vector)
  entropy_pre_corr <- rbind(entropy_pre_corr, data.frame(Cluster = i, Entropy = entropy))
}

write.csv(entropy_pre_corr, "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/clustering/initial_entropy.csv", row.names = TRUE)

###########################################################################################################################
# Postmitotic cells
###########################################################################################################################
# 
# # set markers to RNA
# DefaultAssay(so) <- "RNA"
# 
# # 3d umap of cd200
# # plt_3d_umap_gene(so, out.int.dir, "CD200")
# 
# gene.lst <- c("CD200", "TBR1", "PTPRZ1", "DLX6-AS1", "NRXN3")
# 
# for (gene in gene.lst) {
#   # histogram
#   hist(matrix(so@assays$RNA[gene]), breaks = 50)
#   ggsave(str_interp("${out.postmit.dir}/hist_all_${gene}.png"), width = 10, height = 10)
# 
#   # mark cells
#   so$genep <- ifelse(so@assays$RNA[gene] > 0, "Positive", "Negative")
# 
#   # umap
#   DimPlot(so, group.by = "genep")
#   ggsave(str_interp("${out.postmit.dir}/umap_all_${gene}_pos.png"), width = 10, height = 10)
# 
#   DimPlot(subset(so, dataset == "CHDI"), group.by = "genep")
#   ggsave(str_interp("${out.postmit.dir}/umap_chdi_${gene}_pos.png"), width = 10, height = 10)
# 
#   FeaturePlot(so, features = gene)
#   ggsave(str_interp("${out.postmit.dir}/umap_${gene}.png"), width = 10, height = 10)
# 
#   FeaturePlot(subset(so, dataset == "CHDI"), features = gene)
#   ggsave(str_interp("${out.postmit.dir}/umap_chdi_${gene}.png"), width = 10, height = 10)
# }
# 
# # contingency table 1
# tab.cd200.tbr1 <- as.data.frame(table(so@meta.data[so$dataset == "CHDI", c("CD200p", "TBR1p")]))
# write.csv(tab.cd200.tbr1, file = str_interp("${out.postmit.dir}/table_chdi_cd200p_tbr1p.csv"), quote = FALSE, row.names = FALSE)
# 
# # contingency table 2
# tab.cellannot.cd200 <- as.data.frame(table(so@meta.data[so$dataset == "CHDI", c("cell_annot", "CD200p")]))
# write.csv(tab.cellannot.cd200, file = str_interp("${out.postmit.dir}/table_chdi_cellannot_cd200p.csv"), quote = FALSE, row.names = FALSE)
# 
# # contingency table 3
# tab.cellannot.cd200 <- as.data.frame(table(so@meta.data[, c("cell_annot", "CD200p")]))
# write.csv(tab.cellannot.cd200, file = str_interp("${out.postmit.dir}/table_cellannot_cd200p.csv"), quote = FALSE, row.names = FALSE)
# 
# # clear fields
# so$CD200p <- NULL
# so$TBR1p <- NULL

###########################################################################################################################
# Cluster breakdown
###########################################################################################################################

# # get frequency table
# frq.tab <- as.data.frame(table(so@meta.data[, c("dataset","cell_annot")]))
#
# # get percentages for each dataset
# frq.tab <- frq.tab %>% group_by(dataset) %>% mutate(pct=100*Freq/sum(Freq))
#
# # plot dataset barplot
# ggplot(frq.tab, aes(x = cell_annot, group = dataset, fill = dataset, y = Freq)) +
#   geom_col(position = "stack") +
#   theme(legend.position = "top") +
#   xlab("Cluster") +
#   ylab("Cell count")
# ggsave(str_interp("${out.clust.dir}/barplot_dataset_cellcount_r${clust.res}.pdf"), width = 15, height = 8)
#
# # plot dataset barplot
# ggplot(frq.tab, aes(x = cell_annot, group = dataset, fill = dataset, y = pct)) +
#   geom_col(position = "dodge") +
#   theme(legend.position = "top") +
#   xlab("Cluster") +
#   ylab("Dataset percentage (%)")
# ggsave(str_interp("${out.clust.dir}/barplot_dataset_percent_r${clust.res}.pdf"), width = 15, height = 8)
#
# ## split per day
#
# # get frequency table
# frq.tab <- as.data.frame(table(so@meta.data[, c("dataset","cell_annot","day")]))
#
# # get percentages for each dataset
# frq.tab <- frq.tab %>% group_by(dataset) %>% mutate(pct=100*Freq/sum(Freq))
#
# # plot day barplot
# ggplot(frq.tab, aes(x = cell_annot, group = dataset, fill = day, y = Freq)) +
#   geom_col(position = "stack") +
#   theme(legend.position = "top") +
#   xlab("Cluster") +
#   ylab("Cell count") +
#   facet_wrap(~dataset) +
#   scale_fill_manual(values = col.pal.day)
# ggsave(str_interp("${out.clust.dir}/barplot_day_cellcount_r${clust.res}.pdf"), width = 15, height = 8)
#
# # plot percentage
# ggplot(frq.tab, aes(x = cell_annot, group = dataset, fill = day, y = pct)) +
#   geom_col(position = "stack") +
#   theme(legend.position = "top") +
#   xlab("Cluster") +
#   ylab("Dataset percentage (%)") +
#   facet_wrap(~dataset) +
#   scale_fill_manual(values = col.pal.day)
# ggsave(str_interp("${out.clust.dir}/barplot_day_percent_r${clust.res}.pdf"), width = 15, height = 8)

##########################################################################################################################
Cell cycle
##########################################################################################################################

NOTE: we find that first 3 PCA components are related to cell cycle

if (file.exists(int.file.cc.rds)) {
  so.corr <- readRDS(int.file.cc.rds)
} else {
  # recommended approach (see https://satijalab.org/seurat/articles/cell_cycle_vignette)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so.corr <- ScaleData(so, vars.to.regress = "CC.Difference", features = rownames(so))
  so.corr <- RunPCA(so.corr, features = VariableFeatures(so.corr))
  so.corr <- RunUMAP(so.corr, dims = 1:25, verbose = TRUE, label = TRUE, n.components = 3L)
  so.corr <- RunTSNE(so.corr, dims = 1:25, verbose = TRUE, label = TRUE, n.components = 3L)
  saveRDS(so.corr, file = int.file.cc.rds)
}

# seurat list of genes that contribute to cell cycle
cell_cycle_genes <- cc.genes
s_genes <- as.list(cell_cycle_genes$s.genes)
g2m_genes <- as.list(cell_cycle_genes$g2m.genes)
cycle_genes <- c(s_genes, g2m_genes)

## pc comparison

# feature loadings
loadings_pca_df <- as.data.frame(Loadings(so, reduction = "pca"))
#loadings_pca_corr_df <- as.data.frame(Loadings(so.corr, reduction = "pca"))
abs_loadings_pca_df <- abs(loadings_pca_df)
abs_loadings_pca_corr_df <- abs(loadings_pca_corr_df)

contingency table per PC1 non corrected
pc1_df <- abs_loadings_pca_df[order(-abs_loadings_pca_df$PC_1),]
pc1_df <- pc1_df[,1, drop=FALSE]
top_100_df <- head(pc1_df, 100)

gene_in_top_100 <- sapply(rownames(pc1_df), function(gene) gene %in% rownames(top_100_df))
gene_in_list <- rownames(pc1_df) %in% cycle_genes

contingency_table <- table(gene_in_list, gene_in_top_100)

# contingency tables for each PC non corrected SO
contingency_tables <- list()

for (pc in colnames(abs_loadings_pca_df)[1:20]) {
  # sorting
  pc_df <- abs_loadings_pca_df[order(-abs_loadings_pca_df[[pc]]), , drop = FALSE]
  pc_df <- pc_df[, which(colnames(abs_loadings_pca_df) == pc), drop = FALSE]
  top_50_df <- head(pc_df, 50)
  # checks
  gene_in_top_50 <- sapply(rownames(pc_df), function(gene) gene %in% rownames(top_50_df))
  gene_in_list <- rownames(pc_df) %in% cycle_genes

  contingency_tables[[pc]] <- table(gene_in_list, gene_in_top_50)
}
contingency_tables_df <- do.call(rbind, lapply(contingency_tables, as.data.frame))
write.csv(contingency_tables_df, file = "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/cell_cycle/cell_cycle_gene_table.csv", row.names = TRUE)

# contingency tables for each PC corrected SO
contingency_tables_corr <- list()

for (pc in colnames(abs_loadings_pca_corr_df)[1:20]) {
  # sorting
  pc_corr_df <- abs_loadings_pca_corr_df[order(-abs_loadings_pca_corr_df[[pc]]), , drop = FALSE]
  pc_corr_df <- pc_corr_df[, which(colnames(abs_loadings_pca_corr_df) == pc), drop = FALSE]
  top_50_corr_df <- head(pc_corr_df, 50)
  # checks
  gene_in_top_50_corr <- sapply(rownames(pc_corr_df), function(gene) gene %in% rownames(top_50_corr_df))
  gene_in_list_corr <- rownames(pc_corr_df) %in% cycle_genes

  contingency_tables_corr[[pc]] <- table(gene_in_list_corr, gene_in_top_50_corr)
}


# odd ratio corr and non corr
odds_ratio_grouped <- list()
for (i in 1:20) {
  frequencies <- as.data.frame.matrix(contingency_tables[[i]])
  odds_ratio <- frequencies[1, 1] * frequencies[2, 2] / (frequencies[2, 1] * frequencies[1, 2])
  frequencies_corr <- as.data.frame.matrix(contingency_tables_corr[[i]])
  odds_ratio_corr <- frequencies_corr[1, 1] * frequencies_corr[2, 2] / (frequencies_corr[2, 1] * frequencies_corr[1, 2])
  odds_ratio_grouped[[i]] <- c(odds_ratio, odds_ratio_corr)
}


or.df <- data.frame(
  PC = paste(seq(1, 20)),
  Uncorrected = unlist(lapply(odds_ratio_grouped, FUN = function(x) x[1])),
  Corrected = unlist(lapply(odds_ratio_grouped, FUN = function(x) x[2]))
)
or.df <- reshape2::melt(or.df, value.name = "OR")
or.df$PC <- factor(or.df$PC, levels = factor(paste(seq(1, 20))))
ggplot(or.df) +
  geom_col(aes(x = PC, y = OR, fill = variable), position = "dodge")

# pca 1 i 3 sense corregir i corregit
plot1 <- DimPlot(so, reduction = "pca", group.by = "Phase", dims = c(1, 3))
plot2 <- DimPlot(so.corr, reduction = "pca", group.by = "Phase", dims = c(1, 3))
grid.arrange(plot1, plot2, ncol = 2)

# pcs de nomes gens del cicle celular
# cell_cycle_features_loadings <- list()
# for (i in 1:length(s_genes)) {
# cell_cycle_features_loadings[[i]]<-subset(loadings_pca, rownames(loadings_pca) == s_genes[i])
# }

# pca
DimPlot(so, reduction = "pca", group.by = "Phase")
ggsave(str_interp("${out.cc.dir}/pca_cc_uncorrected.pdf"), width = 20, height = 20)
DimPlot(so.corr, reduction = "pca", group.by = "Phase")
ggsave(str_interp("${out.cc.dir}/umap_cc_corrected.pdf"), width = 20, height = 20)

# umap
DimPlot(so, group.by = "Phase", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.cc.dir}/umap_cc_uncorrected.pdf"), width = 20, height = 20)
DimPlot(so.corr, group.by = "Phase", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.cc.dir}/umap_cc_corrected.pdf"), width = 20, height = 20)

# plot datasets
DimPlot(so.corr, group.by = "dataset", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.cc.dir}/umap_dataset_cc_corrected.pdf"), width = 20, height = 20)
DimPlot(so.corr, group.by = "day", label = TRUE, cols = col.pal.day) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.cc.dir}/umap_day_cc_corrected.pdf"), width = 20, height = 20)
DimPlot(so.corr, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
ggsave(str_interp("${out.cc.dir}/umap_cellannot_cc_corrected.pdf"), width = 20, height = 20)

# # plot 3d
# plt_3d_umap_clust(so.corr, "1", outdir = out.cc.dir, pref = "3d_umap_cc_corr")

##########################################################################################################################
GSEA of PCA
##########################################################################################################################

# iterate over the PC
for (pc in unique(colnames(abs_loadings_pca_df))) {
  cat(str_interp("Processing ${pc}..."), sep = "\n")

  # get list for pc
  col_index <- which(colnames(abs_loadings_pca_df) == pc)
  ordered_df <- abs_loadings_pca_df[order(-abs_loadings_pca_df[, col_index]), ]  
  gene_list <- ordered_df[, col_index]
  names(gene_list) <- rownames(ordered_df)
  
  # run GSEA analysis
  gsea.out <- run_gsea(gene_list, go.file, pval = 0.05)

  # prep output
  gsea.tab.out <- as.data.frame(gsea.out$Results)[, c("pathway", "pval", "padj", "NES", "size", "Enrichment")]

  # save output
  out.pref <- str_interp("${out.annot.mark.dir}/gsea_cat/gsea_seurat_${pc}")
  write.table(gsea.tab.out, file = str_interp("${out.pref}_table.txt.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(str_interp("${out.pref}_plot.pdf"), plot = gsea.out$Plot, width = 10, height = 10)
}

###########################################################################################################################
# UMAPs for various PC combinations (studying CC influence)
###########################################################################################################################

## Here we recluster based on PCA components 4-13. This was decided based on the UMAPs that were generated with the python
## notebook where we saw that including more components leads to more noise.

# dimensions to include starting at PC4
so.cat <- subset(so, dataset =="Cattaneo")
max.dim <- c(5,10,15,20)
  
for(i in max.dim){
  
  # recluster
  so.aux <- FindNeighbors(so.cat, reduction = "pca", dims = 1:5, verbose = TRUE)
  so.aux <- FindClusters(so.aux, verbose = TRUE, resolution = 1)
  
  # umap
  so.aux <- RunUMAP(so.aux, n.components=2, dims=1:5, verbose = TRUE, label = TRUE)
  p <- DimPlot(so.aux, reduction = "umap", group.by = "Phase") #same changing to cell_annot for clusters

  # save the plot
  filename <- paste0("/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/clustering/UMAP_initial_phase_dim_",i,"_cc.png")
  ggsave(filename, plot = p, width = 8, height = 6)
}

rm(list = so.aux)
rm(list = so.cat)

##########################################################################################################################
Cluster markers (annotation)
##########################################################################################################################

cat(str_interp("Finding markers"), sep = "\n")

# set markers to RNA
DefaultAssay(so) <- "RNA"

so.cat <- subset(so, dataset == "Cattaneo")

# get cluster markers
clust.lfc.mark <- get_clust_mark(so.cat, mark.clust.file, out.annot.mark.dir, bmrt.tab.file, int.ver, clust.res, p.filt = FALSE)

###########################################################################################################################
# GSEA of clusters
###########################################################################################################################

cat(str_interp("Doing GSEA"), sep = "\n")

# iterate over the clusters using all markers
for (clust in unique(clust.lfc.mark$Cluster)) {
  cat(str_interp("Processing cluster ${clust}..."), sep = "\n")

  # get list for cluster
  gene.tab.clust <- clust.lfc.mark[clust.lfc.mark$Cluster == clust, ]
  gene_list <- gene.tab.clust$avg_log2FC
  names(gene_list) <- gene.tab.clust$Gene

  # run GSEA analysis
  gsea.out <- run_gsea(gene_list, go.file, pval = 0.05)

  # prep output
  gsea.tab.out <- as.data.frame(gsea.out$Results)[, c("pathway", "pval", "padj", "NES", "size", "Enrichment")]

  # save output
  out.pref <- str_interp("${out.annot.mark.dir}/gsea_cat/gsea_seurat_cluster_${clust}")
  write.table(gsea.tab.out, file = str_interp("${out.pref}_table.txt.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(str_interp("${out.pref}_plot.pdf"), plot = gsea.out$Plot, width = 10, height = 10)
}


###########################################################################################################################
# CLUSTERING AN INTEGRATED OBJECT
###########################################################################################################################

## Here we recluster based on PCA components 4-13. This was decided based on the UMAPs that were generated with the python
## notebook where we saw that including more components leads to more noise.

if(file.exists(int.file.new.rds)){
  
  so <- readRDS(int.file.new.rds)
  
} else {
  
  # subset Cattaneo for velocity and cell fate analysis
  so <- subset(so, dataset == 'Cattaneo')
  
  # recluster
  so <- FindNeighbors(so, reduction = "pca", dims = 4:13, verbose = TRUE)
  so <- FindClusters(so, verbose = TRUE, resolution = 1)
  
  # umap
  so <- RunUMAP(so, dims = 4:13, verbose = TRUE, label = TRUE)
  
  # save
  saveRDS(so, file = int.file.new.rds)
  
}

DimPlot(so, reduction = "umap", group.by = "seurat_clusters")

###########################################################################################################################
# Cell cycle per cluster entropy after correction
###########################################################################################################################

clust<- as.character(c(0:20))
entropy_post_corr <- data.frame(Cluster = character(), Entropy = numeric())


# so.cat.clust0 <- subset(so.cat, cell_annot==clust)
# phase_table <- table(so.cat.clust0@meta.data$Phase)
# phase_probabilities <- as.vector(prop.table(phase_table)) #probabilities vector G1/G2M/S
# log_vector <- log(phase_probabilities, base=3)
# entropy <- sum(-phase_probabilities*log_vector)

for(i in clust){
  so.cat.clust <- subset(so, seurat_clusters==i)
  phase_table <- table(so.cat.clust@meta.data$Phase)
  phase_probabilities <- as.vector(prop.table(phase_table)) #probabilities vector G1/G2M/S
  log_vector <- log(phase_probabilities, base=3)
  entropy <- sum(-phase_probabilities*log_vector)
  entropy_post_corr <- rbind(entropy_post_corr, data.frame(Cluster = i, Entropy = entropy))
}

write.csv(entropy_post_corr, "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/clustering/corrected_cc_entropy.csv", row.names = TRUE)
#entropy_post_corr <- read.csv("/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/clustering/corrected_cc_entropy.csv")
#entropy_pre_corr <- read.csv("/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/clustering/initial_entropy.csv")

#barplot distribution of cells in clusters
freq_table <- table(so@meta.data$Phase, so@meta.data$seurat_clusters)
colors <- c("G1" = "lightblue", "G2M" = "lightgreen", "S" = "lightpink")
barplot(freq_table, 
        main = "Distribution of Cells Along Clusters", 
        xlab = "Clusters", 
        ylab = "Cell number", 
        col = colors,
        ylim = c(0, 27000),
        legend.text = c("S", "G2M", "G1"),
        args.legend = list(x = "topright", bty = "n", fill = colors))

#histogram pre post cc entropy
freq_entropy_pre <- table(cut(entropy_pre_corr$Entropy, breaks = seq(0, 1, by = 0.1)))
freq_entropy_post <- table(cut(entropy_post_corr$Entropy, breaks = seq(0, 1, by = 0.1)))
barplot(rbind(as.numeric(freq_entropy_pre), as.numeric(freq_entropy_post)), 
        beside = TRUE,
        names.arg = c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6", "0.6-0.7","0.7-0.8", "0.8-0.9","0.9-1"),
        col = c("skyblue", "lightpink"),
        xlab = "Entropy",
        ylab = "Cluster frequency",
        main = "Entropy distribution")
legend("topright", legend = c("Entropy Pre", "Entropy Post"), fill = c("skyblue", "lightpink"))

#boxplot comparison
boxplot(entropy_pre_corr$Entropy, entropy_post_corr$Entropy, 
        names = c("Pre Cell Cycle Correction", "Post Cell Cycle Correction"), 
        col = c("skyblue", "lightpink"),  
        main = "Boxplot comparison between entropies",   
        xlab = "Clusters", ylab = "Entropy")  

##########################################################################################################################
Cluster markers after removing CC PCA components
##########################################################################################################################

# Here we do marker analysis with new clusters found working with PCs 4-13.

cat(str_interp("Finding markers"), sep = "\n")

# set markers to RNA
DefaultAssay(so) <- "RNA"

so.cat <- subset(so, dataset == "Cattaneo")

# get cluster markers
clust.lfc.mark <- get_clust_mark(so.cat, mark.clust.file.new, out.annot.mark.dir, bmrt.tab.file, int.ver, clust.res, p.filt = FALSE)

###########################################################################################################################
# Removing contaminated clusters
###########################################################################################################################

# remove contaminant populations
clust.keep <- paste(0:20)
clust.keep <- clust.keep[!clust.keep %in% c("2","15","19","20")]
so <- subset(so, seurat_clusters %in% clust.keep)

#recompute umap
so <- RunUMAP(so, dims = 4:13, verbose = TRUE, label = TRUE)
DimPlot(so, reduction = "umap", group.by = "Phase")

###########################################################################################################################
# SAVE FILE FOR SCVELO
###########################################################################################################################

# file for scvelo
fn <- str_interp("${scvelo.dir}/seurat_export_for_scvelo_${int.ver}.tsv")

# store
if (!file.exists(fn)) {
  
  # clustering information
  clust.df <- data.frame(
    cell = rownames(so@meta.data),
    sample = so@meta.data[, c("SampleID")],
    phase = so@meta.data[, c("Phase")],
    cell_annot = so@meta.data[, c("seurat_clusters")]
  )
  
  # add cell embedding info
  pca.df <- as.data.frame(so@reductions$pca@cell.embeddings)
  pca.df <- pca.df[,4:13]
  colnames(pca.df) <- paste("PC_",1:10,sep="")
  pca.df <- tibble::rownames_to_column(pca.df, var = "cell")
  clust.df <- merge(clust.df, pca.df, by = "cell")
  
  # add UMAP information
  umap.df <- as.data.frame(so@reductions$umap@cell.embeddings)
  umap.df <- tibble::rownames_to_column(umap.df, var = "cell")
  clust.df <- merge(clust.df, umap.df, by = "cell")
  
  # store
  write.table(clust.df, file = fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
