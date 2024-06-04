#############################################################################################################################
# Integration of human data
#############################################################################################################################

rm(list = ls())
graphics.off()

# set seed
set.seed(1)

# load deps
library(dplyr)
library(tibble)
library(Seurat)
library(stringr)
library(ggplot2)
library(pheatmap)
library(mixtools)
library(patchwork)
library(data.table)
library(scCustomize)
library(RColorBrewer)

#############################################################################################################################
# Input/output directories/files
#############################################################################################################################

# source code
creatio.data.dir <- "/gpfs/projects/bsc83/MN4/bsc83/Data/Creatio/" # creatio.data.dir <- "~/marenostrum/data_creatio/"
creatio.proj.dir <- "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/" # creatio.proj.dir <- "~/marenostrum/projects_creatio/"
source("/gpfs/projects/bsc83/MN4/bsc83/Data/Creatio/code/src/seurat/seurat_src.R")

# project dirs (assumes there is a cellranger folder)
# scc17_dir <- str_interp("${creatio.proj.dir}/SCCANALS17")
# scc18_dir <- str_interp("${creatio.proj.dir}/SCCANALS18")
chdi_dir <- str_interp("${creatio.proj.dir}/SCC_050810")
catt_dir <- str_interp("${creatio.proj.dir}/CattaneoPaper")
crdf_dir <- str_interp("${creatio.proj.dir}/cardiff_sc_data")
grch38.ref.dir <- str_interp("${creatio.data.dir}/references/refdata-gex-GRCh38-2020-A")

# integration version
args <- commandArgs(trailingOnly = TRUE)
int.ver <- args[1] # int.ver="v6"
clust.res <- args[2] # clust.res="1"

# output dir
out.bse.dir <- str_interp("${creatio.proj.dir}integration_human/${int.ver}/") # base directory
out.qc.dir <- str_interp("${out.bse.dir}/qc") # qc directory
out.int.dir <- str_interp("${out.bse.dir}integration") # integration directory
out.clust.dir <- str_interp("${out.bse.dir}/clustering") # clustering directory
out.mark.dir <- str_interp("${out.bse.dir}/markers/") # marker directory
scvelo.dir <- str_interp("${out.bse.dir}/scvelo")

# markers subdir
grp.plt.dir <- str_interp("${out.mark.dir}/grouped/")
ind.plt.dir <- str_interp("${out.mark.dir}/individual/")

# integrated file
int.file.rds <- str_interp("${out.int.dir}/integration_human.rds")
int.filt.file.rds <- str_interp("${out.int.dir}/integration_human_filt.rds")
int.feats.file <- str_interp("${out.int.dir}/integration_human_features.txt")
mark.clust.file <- str_interp("${out.mark.dir}/markergenes_${clust.res}.rds")

# create folder
dir.create(out.bse.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.qc.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.int.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.clust.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.mark.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scvelo.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(grp.plt.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ind.plt.dir, recursive = TRUE, showWarnings = FALSE)

# other files
bmrt.tab.file <- str_interp("${creatio.data.dir}/biomart/hsapiens_ensembl_table.tsv")

# regressing out 
regout <- c("percent.mt", "percent.rp", "percent.hb") # mito, ribosamal, hemboglobin

# info
print(str_interp("Generating integrated version ${int.ver} in:"))
print(str_interp("${out.bse.dir}"))
print(str_interp("Regressing out ${regout}"))

####################################################################################################################################
# Integration
####################################################################################################################################

if (!file.exists(int.file.rds)) {
  cat("Integrated file does not exist. Generating it...", sep = "\n")

  ## read in transplant data (only 9 human cells are found - so not including)

  # read in transplant data
  # scc17_data <- readin_sn_cellranger(
  #   scc17_dir, grch38.ref.dir, analysis=TRUE, mix.species=TRUE, soupx = TRUE,
  #   anl.spec="GRCh38", sex.col="sex.human", metafile="metadata_transplant"
  # )
  # scc18_data <- readin_sn_cellranger(
  #   scc18_dir, grch38.ref.dir, analysis=TRUE, mix.species=TRUE, soupx = TRUE,
  #   anl.spec="GRCh38", sex.col="sex.human", exclude=c("AY7505","AZ0222","AZ0223","AZ1400","AZ1401")
  # )

  # rename project column
  # scc17_data@meta.data <- rename(scc17_data@meta.data, dataset = project)
  # scc18_data@meta.data <- rename(scc18_data@meta.data, dataset = project)

  ## Cattaneo data

  # read in Cattaneo sc data
  catt_data <- readin_cellranger(catt_dir, doublet.rate = 0.1, soupx = TRUE, subsample = FALSE)

  # add covariate info
  catt_data@meta.data$genotype <- "WT"
  catt_data@meta.data <- catt_data@meta.data %>% rename(day = date)

  # drop unnecessary columns
  columns_to_drop <- c("week", "barcode")
  catt_data@meta.data <- catt_data@meta.data[, !(names(catt_data@meta.data) %in% columns_to_drop)]

  # drop cortical cells
  catt_data <- subset(catt_data, day != "8 week 4 day")

  # drop 7 week samples as well
  catt_data <- subset(catt_data, day != "7 week")

  ## CHDI data

  # read in CHDI sc data (this data is just concatenated)
  load("/gpfs/projects/bsc83/Projects/Creatio/chdi_scrnaseq/HD_DATA_050810.RDATA")
  chdi_data <- CreateSeuratObject(counts = Data_050810, meta.data = Meta_050810, project = "CHDI", min.cells = 3, min.features = 500)

  # rename columns
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(plate_num = Plate_Number)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(SampleID = SampleName)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(cell_line = Cell_line)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(genotype = Cell_type)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(dataset = orig.ident)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(tech_rep = Tech_rep)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(bio_rep = Bio_rep)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(plate = Plate)
  chdi_data@meta.data <- chdi_data@meta.data %>% rename(div = Day)

  # refactor
  chdi_data@meta.data$div <- factor(chdi_data@meta.data$div, levels = c("8", "12", "16"))
  chdi_data@meta.data$genotype <- ifelse(chdi_data@meta.data$genotype == "CTR", "WT", chdi_data@meta.data$genotype)

  # replace SampleID
  chdi_data@meta.data$SampleID <- paste(
    chdi_data@meta.data$cell_line,
    chdi_data@meta.data$genotype,
    chdi_data@meta.data$div,
    chdi_data@meta.data$bio_rep,
    chdi_data@meta.data$tech_rep,
    sep = "_"
  )

  # drop unnecessary columns
  columns_to_drop <- c("flowcell", "lane", "index", "id")
  chdi_data@meta.data <- chdi_data@meta.data[, !(names(chdi_data@meta.data) %in% columns_to_drop)]

  ## Cell cycle

  # get gene list from Tirosh et al 2015 (G2/M phase)
  s.genes <- firstup(cc.genes$s.genes)
  g2m.genes <- firstup(cc.genes$g2m.genes)

  # assign cell cycle score (S.Score, G2M.Score)
  chdi_data <- CellCycleScoring(chdi_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  catt_data <- CellCycleScoring(catt_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  # QC metrics (Cattaneo is done when reading in)
  chdi_data <- comp_qc_metrics(chdi_data)

  # clean environment
  rm(list = c("Data_050810", "Meta_050810"))
  gc()

  ## common genes

  # create list of gene names from both objects (use genes in common between single-cell data)
  total.genes <- list(rownames(catt_data), rownames(chdi_data))

  # Get common gene names
  common.genes <- Reduce(f = intersect, x = total.genes)

  # keep only common genes
  catt_data <- catt_data[common.genes, ]
  chdi_data <- chdi_data[common.genes, ]

  ## qc

  # prefiltering plots
  # qc_plots(catt_dir, catt_data, suffix = "_prefilt")
  # qc_plots(chdi_dir, chdi_data, suffix = "_prefilt")

  # qc subsetting (taken from from F Londoño's code)
  catt_data <- qc_subset(catt_data)
  chdi_data <- qc_subset(chdi_data, nCount_RNA_min = 50000, nCount_RNA_max = 1200000)

  # postfiltering qc
  # qc_plots(catt_dir, catt_data, suffix = "_postfilt")
  # qc_plots(chdi_dir, chdi_data, suffix = "_postfilt")

  # save data as checkpoint
  # save.image(file = str_interp("${bse.dir}/integration_human_midcheckpoint_${int.ver}.0.Rdata"))
  # load(str_interp("${bse.dir}/integration_human_midcheckpoint_${int.ver}.0.Rdata"))

  # normalize data (LogNormalize by default)
  catt_data <- NormalizeData(catt_data, verbose = TRUE)
  chdi_data <- NormalizeData(chdi_data, verbose = TRUE)

  # find variable features
  catt_data <- FindVariableFeatures(catt_data, verbose = TRUE, nfeatures = 2500)
  chdi_data <- FindVariableFeatures(chdi_data, verbose = TRUE, nfeatures = 2500)

  # create list
  scc_data.list <- list(catt = catt_data, chdi = chdi_data)

  # choose reference samples
  # reference_datasets <- which(names(scc_data.list) %in% c("catt"))

  # allocate future memory
  options(future.globals.maxSize = 8000 * 1024^2)

  # find features that are repeatedly ranked as variable within each dataset
  print("Finding integration features (variable genes)...")
  features <- SelectIntegrationFeatures(object.list = scc_data.list, nfeatures = 2500, verbose = TRUE)
  features <- features[grep("^[Mm]t", features, invert = TRUE)] # remove mitochondrial genes
  features <- features[grep("^[Rr]p", features, invert = TRUE)] # remove ribosomal genes

  # store list of genes
  write.table(features, file = int.feats.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # fin integration anchors
  scc.anchors <- FindIntegrationAnchors(
    object.list = scc_data.list,
    # reference = reference_datasets,
    anchor.features = features,
    reduction = "cca",
    dims = 1:25
  )

  # integrate data
  scc_data_int <- IntegrateData(anchorset = scc.anchors, dims = 1:25, verbose = TRUE)

  # scale again together
  scc_data_int <- ScaleData(scc_data_int, verbose = TRUE, vars.to.regress = regout)

  # dimensionality reduction
  scc_data_int <- RunPCA(scc_data_int, verbose = TRUE)
  scc_data_int <- RunUMAP(scc_data_int, dims = 1:25, verbose = TRUE, label = TRUE, n.components = 3L)
  scc_data_int <- RunTSNE(scc_data_int, dims = 1:25, verbose = TRUE, label = TRUE, n.components = 3L)

  # clustering
  scc_data_int <- FindNeighbors(scc_data_int, dims = 1:25, verbose = TRUE)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 0.25)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 0.50)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 0.75)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 1)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 1.25)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 1.5)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 1.75)
  scc_data_int <- FindClusters(scc_data_int, verbose = TRUE, resolution = 2)

  # drop unnecessary columns
  columns_to_drop <- c("orig.ident", "old.ident")
  scc_data_int@meta.data <- scc_data_int@meta.data[, !(names(scc_data_int@meta.data) %in% columns_to_drop)]

  # change name of day values
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "7 week"] <- "7w"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "7 week 4 day"] <- "7w 4d"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "8 week 1 day"] <- "8w 1d"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "8 week 4 day"] <- "8w 4d"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "9 week"] <- "9w"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "9 week 6 day"] <- "9w 6d"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "10 week 6 day"] <- "10w 6d"
  scc_data_int@meta.data$day[scc_data_int@meta.data$day == "11 week 3 day"] <- "11w 3d"

  # refactor day variable
  scc_data_int@meta.data$day <- factor(
    scc_data_int@meta.data$day,
    levels = c("7w", "7w 4d", "8w 1d", "8w 4d", "9w", "9w 6d", "10w 6d", "11w", "11w 3d", "13w")
  )

  # save integrated data
  saveRDS(scc_data_int, file = int.file.rds)
} else {
  cat("Integrated file does exist. Loading it...", sep = "\n")

  # read in integrated file
  scc_data_int <- readRDS(int.file.rds)
}

######################################################################
#CLUSTERING
######################################################################
# select resolution for active identity
Idents(scc_data_int) <- str_interp("integrated_snn_res.${clust.res}")

# assign clustering with desired resolution
scc_data_int@meta.data$cell_annot <- scc_data_int@active.ident
  
# refactor
scc_data_int@meta.data$cell_annot <- factor(scc_data_int@meta.data$cell_annot,paste(0:length(unique(scc_data_int@meta.data$cell_annot))))
  
# fix day levels (missing samples)
pcw.levels <- levels(scc_data_int$day)[table(scc_data_int$day)>0]
scc_data_int@meta.data$day <- factor(scc_data_int@meta.data$day,levels=pcw.levels)

# refactor div
scc_data_int@meta.data$div <- factor(scc_data_int@meta.data$div,levels=c("8","12","16"))
levels(so_chdi@meta.data$div) <- c("8 (DIV)","12 (DIV)","16 (DIV)")

# merge day and div, refactor, and create palette
scc_data_int@meta.data$day <- factor(scc_data_int@meta.data$day,levels=c(levels(scc_data_int@meta.data$day),levels(scc_data_int@meta.data$div)))
scc_data_int@meta.data$day[scc_data_int@meta.data$dataset=="CHDI"] <- scc_data_int@meta.data$div[scc_data_int@meta.data$dataset=="CHDI"]

# color palette for day
col.pal.pcw <- colorRampPalette(c("lightblue", "darkblue"))(length(pcw.levels))
col.pal.div <- colorRampPalette(c("lightgreen", "darkgreen"))(3)
col.pal.day <- c(col.pal.pcw,col.pal.div)

# cell percentage per day in each cluster

so_chdi <- subset(scc_data_int, dataset == "CHDI")

so_chdi_percentage <- so_chdi@meta.data %>%
  group_by(cell_annot, div) %>%
  summarize(cell_num = n())

so_chdi_tot <- so_chdi_percentage %>%
  group_by(cell_annot) %>%
  summarize(cell_total = sum(cell_num))

so_chdi_percentage <- left_join(so_chdi_percentage, so_chdi_tot, by = "cell_annot")

so_chdi_percentage <- so_chdi_percentage %>%
  mutate(percentage = cell_num / cell_total * 100)

write.csv(so_chdi_percentage, "/gpfs/projects/bsc83/Data/Creatio/code/src/seurat/cell_percentage_per_div.csv")


############################################################################################################################
# General plots
############################################################################################################################

# select resolution for active identity
Idents(scc_data_int) <- str_interp("integrated_snn_res.${clust.res}")

# plot each day separately
DimPlot(scc_data_int, split.by = "day", raster = TRUE) + NoLegend()
ggsave(str_interp("${out.int.dir}/umap_split_day.png"), height = 10, width = 25)

# plot each day separately within Cardiff data
DimPlot(subset(scc_data_int, dataset == "Cardiff"), split.by = "day", raster = TRUE) + NoLegend()
ggsave(str_interp("${out.int.dir}/umap_cardiff_split_day.png"), height = 10, width = 15)

# groupby
plt_grpby <- c(
  "seurat_clusters", "dataset", "bio_rep", "tech_rep", "plate_num",
  "genotype", "day", "cell_line", "Phase", "div"
)

# plot for each covariate
for (x in plt_grpby) {
  if (x %in% colnames(scc_data_int@meta.data)) {
    label.flag <- if (x == "seurat_clusters") TRUE else FALSE
    p <- DimPlot(scc_data_int, group.by = x, label = label.flag) + theme(legend.position = "top") +
      guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
    ggsave(str_interp("${out.int.dir}/umap_${x}.png"), plot = p, height = 10, width = 10)
  }
}

# plot for each resolution
for (x in c("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2")) {
  # clustering
  grpby <- str_interp("integrated_snn_res.${x}")
  p <- DimPlot(scc_data_int, group.by = grpby, label = TRUE) + theme(legend.position = "top") +
    guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
  ggsave(str_interp("${out.clust.dir}/umap_seurat_clusters_${x}.png"), plot = p, height = 10, width = 10)

  # 3D UMAPs
  plt_3d_umap_clust(scc_data_int, out.clust.dir, res = x)
}

# get sample metadata
samp.meta <- unique(scc_data_int@meta.data[, c("SampleID", "dataset", "day", "div", "genotype", "cell_line")])
rownames(samp.meta) <- NULL
samp.meta <- samp.meta %>% column_to_rownames("SampleID")

# get clustering breakdown
clst.tab <- table(scc_data_int@meta.data[, c("SampleID", str_interp("integrated_snn_res.${clust.res}"))])
n.clst.tab <- clst.tab / rowSums(clst.tab)
ln.clst.tab <- log(n.clst.tab + 1)
p <- pheatmap(ln.clst.tab, annotation_row = samp.meta, show_rownames = FALSE)
ggsave(str_interp("${out.int.dir}/heatmap_samples_cluster_breakdown_r${clust.res}.png"), plot = p, height = 15, width = 15)

#############################################################################################################################
# Known markers
#############################################################################################################################

# known neuron markers
neurons <- sapply(c(
  "Tubb3", "Map2", "Gadd45g", "Gad2", "Sst", "Npy", "Ppp1r1b", "Bcl11b",
  "Dcx", "Mapt", "Tac1", "Gpr88", "Rarb", "Foxp1", "Foxp2"
), toupper)
early_npcs <- sapply(c("Nes", "Vim", "Btg2", "Gsx2", "Ascl1", "Hes5", "Ttyh1", "Cd99"), toupper)
late_npcs <- sapply(c(
  "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Ebf1", "Isl1", "Sp9", "Sp8",
  "Lhx6", "Nkx2.1", "Helt", "Drd1"
), toupper)
neuroblasts <- sapply(c("Ikzf1", "Ikzf2", "Zfp503", "Zfp521", "Stmn2", "Six3", "Lhx8", "Drd2", "Penk", "Cd200"), toupper)
cortical <- sapply(c("Neurod6", "Neurod2", "Neurog1", "Neurog2", "Tbr1"), toupper)
dev.gen.lst <- list(neurons = neurons, early_npcs = early_npcs, late_npcs = late_npcs, neuroblasts = neuroblasts, cortical = cortical)
plt.feat.umaps(scc_data_int, dev.gen.lst, grp.plt.dir, ind.plt.dir, int.ver)

# other cell types (Olga list)
endothelial <- sapply(c("Pecam1", "Vwf", "Cd31", "Cd144", "Drd1"), toupper)
microglia <- sapply(c("Itgam", "Cx3cr1", "Aif1", "Tmem119", "Cd200r"), toupper)
ependymal <- sapply(c("Foxj1", "Pifo", "Rabl2", "Rsph1"), toupper)
olga.gen.lst <- list(endothelial = endothelial, microglia = microglia, ependymal = ependymal)
plt.feat.umaps(scc_data_int, olga.gen.lst, grp.plt.dir, ind.plt.dir, int.ver)

# markers from Straccia et al
pallium_straccia <- sapply(c("Bcl11b", "Ngn2", "Emx1", "Nr2f1", "Emx2", "Pax6", "Eomes", "Rarg", "Ngn1", "Tbr1", "Slc17a6", "Slc17a7"), toupper)
subpallium_straccia <- sapply(c(
  "Ascl1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Gsx1", "Gsx2", "Nkx2-1", "Rara", "Aldh1a3", "Ebf1", "Etv1",
  "Foxp1", "Foxp2", "Gpr6", "Ikzf1", "Ikzf2", "Isl1", "Meis2", "Sp8", "Znf503", "Calb2", "Gad2", "Lhx1",
  "Npy", "Pvalb", "Sst", "Arpp21", "Calb1", "Drd1", "Drd2", "Gpr88", "Oprm1", "Penk", "Ppp1r1b", "Rarb",
  "Tac1", "Nr2f2", "Gli1", "Lhx6", "Nkx6-2"
), toupper)
early_straccia <- sapply(c("Dach1", "Dbx2", "Gli3", "Msi1", "Nes", "Otx1", "Otx2", "Pou3f1", "Six3", "Sox1", "Sox2", "Tlx1", "Tjp1", "Zbtb16"), toupper)
midbrain_straccia <- sapply(c("Pax2", "Chrm4", "Nr4a2", "Th", "Lmx1b"), toupper)
# spinal_straccia <- sapply(c("Chat", "Foxa2", "Mnx1"),toupper)
pan_neuronal_straccia <- sapply(c("Cdh2", "Dcx", "Map2", "Ncam1", "Rbfox3", "Tubb3"), toupper)
glial_straccia <- sapply(c("Alf1", "Cnp", "Gfap", "Olig1", "Olig2", "S100b"), toupper)
straccia.gen.lst <- list(
  pallium_straccia = pallium_straccia, subpallium_straccia = subpallium_straccia, early_straccia = early_straccia,
  midbrain_straccia = midbrain_straccia, pan_neuronal_straccia = pan_neuronal_straccia,
  glial_straccia = glial_straccia
)
plt.feat.umaps(scc_data_int, straccia.gen.lst, grp.plt.dir, ind.plt.dir, int.ver)

# CHAT interneurons markers (Giralt)
chat_mature <- c("CHAT","VACHT","ACHE","GRIN3A","PRPH","FBN2")
chat_dev <- c("GBX2", "ZIC4", "LHX6", "FGF8","FGF17", "DBX1")
chat.lst <- list(chat_mature=chat_mature, chat_dev=chat_dev)
plt.feat.umaps(scc_data_int, chat.lst, grp.plt.dir, ind.plt.dir, int.ver)

# FeaturePlot(scc_data_int,features = c("FGF8","AIGF","HBGF-8","FGF-8","KAL6","HH6"))
# FeaturePlot(scc_data_int,features = c("DBX1","HLX1"))
