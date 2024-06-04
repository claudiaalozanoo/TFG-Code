#########################################################################################################################################################
# Dependancies
#########################################################################################################################################################
suppressMessages(require(xlsx))
suppressMessages(require(RANN))
suppressMessages(require(Matrix))
suppressMessages(require(dplyr))
suppressMessages(require(fgsea))
suppressMessages(require(SoupX))
suppressMessages(require(qusage))
suppressMessages(require(plotly))
suppressMessages(require(gridExtra))
suppressMessages(require(stringr))
suppressMessages(require(htmlwidgets))
suppressMessages(require(DoubletFinder))
#########################################################################################################################################################
# Read in functions
#########################################################################################################################################################

# reads in cellranger output for single-cell data from project directory
readin_cellranger <- function(project_dir, exclude = c(), subsample = FALSE, doublet.rate = 0.05, soupx = TRUE) {
  
  # testing
  project_name <- basename(project_dir)
  cellranger_dir <- paste(project_dir, "cellranger", sep = "/")
  metadata_path <- paste(project_dir, "metadata.tsv", sep = "/")

  # get samples available in cellranger directory
  samples <- list.dirs(path = cellranger_dir, full.names = FALSE, recursive = FALSE)
  
  # if subsample
  if (subsample) {
    samples <- sample(x = samples, size = 10, replace = FALSE)
  }

  # read in metadata
  sample_metadata <- fread(metadata_path, sep = "\t", colClasses = "character", data.table = FALSE)
  rownames(sample_metadata) <- sample_metadata$barcode

  # remove exclude samples
  keepers <- samples[!samples %in% exclude]

  # remove samples not in metadata
  keepers <- keepers[keepers %in% sample_metadata$barcode]
  miss_meta <- keepers[!keepers %in% sample_metadata$barcode]
  if (length(miss_meta) > 0) {
    print(stringr::str_interp("Excluding ${miss_meta}: missing in metadata"))
  }
  
  # generate list of seurat objects
  seurat_list <- lapply(keepers, function(x) {
    
    if (soupx) {
      
      cat(str_interp("Loading 10X counts for ${x}"),sep="\n")
      
      # load matrix
      sc.data <- load10X(str_interp("${cellranger_dir}/${x}/outs/"))
      
      # load clustering
      if(file.exists(str_interp("${cellranger_dir}/${x}/outs/analysis/clustering/graphclust/clusters.csv"))){
        cat("Found cellranger clustering file",sep="\n")
        cr.clusts.df <- read.csv(str_interp("${cellranger_dir}/${x}/outs/analysis/clustering/graphclust/clusters.csv"),header = TRUE)
      } else if (file.exists(str_interp("${cellranger_dir}/${x}_sec/outs/analysis/clustering/graphclust/clusters.csv"))){
        cat("Found cellranger clustering file",sep="\n")
        cr.clusts.df <- read.csv(str_interp("${cellranger_dir}/${x}_sec/outs/analysis/clustering/graphclust/clusters.csv"),header = TRUE)
      } else {
        cat("Did not find cellranger clustering file for sample {x}",sep="\n")
        return(NULL)
      }
      
      # set clusters
      cr.clusts <- cr.clusts.df$Cluster
      names(cr.clusts) <- cr.clusts.df$Barcode
      sc.data <- setClusters(sc.data,clusters = cr.clusts)
      
      # autoestimate contaminant fraction
      sc.data <- autoEstCont(sc.data, forceAccept =TRUE)
      
      # exclude sample if fraction is larger than 50%
      if (sc.data$fit$rhoEst>0.5){
        return(NULL)
      }
      
      # filter out contaminant cells
      sc.data <- adjustCounts(sc.data)
      
      # create seurat object
      cur_seurat <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = project_name)
      cur_seurat$SampleID <- x
      
    } else {
      
      cat("Skipping SoupX",sep="\n")
      
      # read in
      cur_seurat <- Read10X(data.dir = paste(cellranger_dir, x, "outs/filtered_feature_bc_matrix", sep = "/"))
      cur_seurat <- CreateSeuratObject(counts = cur_seurat, min.cells = 3, min.features = 500, project = project_name)
      cur_seurat$SampleID <- x
      
    }
    
    # get cell metadata
    cell_metadata <- sample_metadata[x, ]
    for (meta in names(cell_metadata)) {
      cur_seurat[[meta]] <- cell_metadata[[meta]]
    }

    # add qc metrics
    cur_seurat <- comp_qc_metrics(cur_seurat)

    # filter out doublets
    if (is.numeric(doublet.rate)) {
      cur_seurat <- filter_out_doublets(cur_seurat, doublet.rate = doublet.rate)
      cur_seurat$doubletfinder_pann <- NULL
      cur_seurat$doubletfinder <- NULL
    }

    # get map estimate of cell sex
    if (("sex" %in% colnames(cell_metadata$meta.data)) && (cell_metadata$sex == "Female_Male")) {
        cur_seurat <- get_sex_map(cur_seurat)
    }
    
    # remove unnecessary columns if they exist
    if("orig.ident" %in% colnames(cur_seurat)) cur_seurat$orig.ident <- NULL
    if("barcode" %in% colnames(cur_seurat)) cur_seurat$barcode <- NULL
    if("id" %in% colnames(cur_seurat)) cur_seurat$id <- NULL

    # add sample prefix to cell ID
    # colnames(cur_seurat) <- paste0(str_interp("${x}_"), colnames(cur_seurat))
    
    # append sample name to
    return(cur_seurat)
  })

  # set names
  seurat_list <- setNames(seurat_list, keepers)

  # remove NULLs
  seurat_list <- seurat_list[lapply(seurat_list, is.null) == FALSE]
  
  # merge objects
  seurat_obj <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

  # store sex breakdown
  if (("sex" %in% colnames(seurat_obj)) && any(sample_metadata$sex == "Female_Male")) {
    
    sex.tab <- as.data.frame.matrix(table(seurat_obj[[c("SampleID", "sex")]]))
    sex.tab$Sample <- rownames(sex.tab)
    sex.tab <- sex.tab[sex.tab$Sample %in% sample_metadata$barcode[sample_metadata$sex == "Female_Male"], ]
    write.table(
      sex.tab[, c("Sample", "Female", "Male")],
      file = stringr::str_interp("${seu_dir}/cell_sex_table.tsv"),
      quote = FALSE, row.names = FALSE
    )
    
    # filter out female cells
    seurat_obj <- subset(seurat_obj, sex == "Male")
  }

  # read in metadata
  rm(seurat_list)
  gc()

  # this needs to be done with all samples: issues in postnatal samples
  # # get gene list from Tirosh et al 2015 (G2/M phase)
  # s.genes <- firstup(cc.genes$s.genes)
  # g2m.genes <- firstup(cc.genes$g2m.genes)
  # 
  # # assign cell cycle score (S.Score, G2M.Score)
  # seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  return(seurat_obj)
}

# reads in cellranger output from single-nuclei data from project directory
readin_sn_cellranger <- function(project_dir, ref, exclude = c(), subsample = FALSE, soupx = TRUE, analysis=FALSE, anl.spec="GRCh38", mix.species=FALSE, metafile=NULL, sex.col="sex") {
  
  # define paths
  project_name <- basename(project_dir)
  cellranger_all_dir <- paste(project_dir, "cellranger_all_reads", sep = "/")
  cellranger_exonic_dir <- paste(project_dir, "cellranger_exonic_reads", sep = "/")
  polya_path <- stringr::str_interp("${ref}/polyA/polyA.tsv.gz")
  
  # read in metadata
  if (is.null(metafile)){
    metadata_path <- str_interp("${project_dir}/metadata.tsv")
  } else {
    metadata_path <- str_interp("${project_dir}/${metafile}.tsv")
  }

  # get samples available in cellranger directory
  samples <- list.dirs(path = cellranger_all_dir, full.names = FALSE, recursive = FALSE)

  # if subsample
  if (subsample) {
    samples <- sample(x = samples, size = 10, replace = FALSE)
  }

  # read in polyA information
  polya.tab <- read.table(polya_path, header = TRUE, sep = "\t")
  polya.tab <- polya.tab[!duplicated(polya.tab), ]
  polya.tab <- polya.tab[, c("gene", "polyA")]
  polya.tab <- aggregate(polyA ~ gene, data = polya.tab, sum)
  rownames(polya.tab) <- polya.tab$gene

  # read in metadata
  sample_metadata <- fread(metadata_path, sep = "\t", colClasses = "character", data.table = FALSE)
  rownames(sample_metadata) <- sample_metadata$barcode

  # remove exclude samples
  keepers <- samples[!samples %in% exclude]

  # remove samples not in metadata
  keepers <- keepers[keepers %in% sample_metadata$barcode]
  miss_meta <- keepers[!keepers %in% sample_metadata$barcode]
  if (length(miss_meta) > 0) {
    cat(stringr::str_interp("Excluding ${miss_meta}: missing in metadata"),sep="\n")
  }

  # generate list of seurat objects
  seurat_list <- lapply(keepers, function(x) {
    
    cat(str_interp("Pre-processing in sample: ${x}"), sep = "\n")

    ## read in data

    if (soupx) {
      
      cat("Running SoupX for (i) all counts and (ii) exonic only", sep = "\n")
      
      # read in and run SoupX      
      seu.all <- readin_sn_soupx_sample(x, cellranger_all_dir, project_name, analysis, anl.spec)
      seu.exon <- readin_sn_soupx_sample(x, cellranger_exonic_dir, project_name, analysis, anl.spec)
      
      # return NULL if either is NULL
      if(is.null(seu.all) || is.null(seu.exon)){
        return(NULL)
      }
      
      
    } else {
      
      # this part is under construction: missing species filter
      
      # exonic and intronic
      seu.all <- Read10X(data.dir = stringr::str_interp("${cellranger_all_dir}/${x}/outs/filtered_feature_bc_matrix"))
      seu.all <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 50, project = project_name)
      seu.all$SampleID <- x
      
      # exonic reads only
      seu.exon <- Read10X(data.dir = stringr::str_interp("${cellranger_exonic_dir}/${x}/outs/filtered_feature_bc_matrix"))
      seu.exon <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 50, project = project_name)
      seu.exon$SampleID <- x
      
    }
    
    # if mixed species keep only relevant genes and adapt names
    if(mix.species){
      
      # get species genes
      spec.genes.all <- grep(anl.spec,rownames(seu.all@assays$RNA@counts))
      spec.genes.exon <- grep(anl.spec,rownames(seu.exon@assays$RNA@counts))
      
      # filter matrices
      seu.all@assays$RNA@counts <- seu.all@assays$RNA@counts[spec.genes.all,]
      seu.all@assays$RNA@data <- seu.all@assays$RNA@data[spec.genes.all,]
      seu.exon@assays$RNA@counts <- seu.exon@assays$RNA@counts[spec.genes.exon,]
      seu.exon@assays$RNA@data <- seu.exon@assays$RNA@data[spec.genes.exon,]
      
      # remove prefix from each gene
      old.rownames <- rownames(seu.all@assays$RNA@counts)
      new.rownames <- sub(str_interp("^${anl.spec}---"), "", old.rownames)
      rownames(seu.all@assays$RNA@counts) <- new.rownames
      rownames(seu.all@assays$RNA@data) <- new.rownames
      old.rownames <- rownames(seu.exon@assays$RNA@counts)
      new.rownames <- sub(str_interp("^${anl.spec}---"), "", old.rownames)
      rownames(seu.exon@assays$RNA@counts) <- new.rownames
      rownames(seu.exon@assays$RNA@data) <- new.rownames
      
    }
    
    # get missing genes and cells in exon only matrix
    miss.genes.exon <- setdiff(rownames(seu.all), rownames(seu.exon))
    miss.cells.exon <- setdiff(colnames(seu.all), colnames(seu.exon))
    
    # missing genes
    miss.genes.exon.mat <- Matrix(0, nrow = length(miss.genes.exon), ncol = length(colnames(seu.exon)), sparse = TRUE)
    rownames(miss.genes.exon.mat) <- miss.genes.exon
    seu.exon@assays$RNA@counts <- rbind2(seu.exon@assays$RNA@counts, miss.genes.exon.mat)

    # missing cells
    miss.cells.exon.mat <- Matrix(0, nrow = nrow(seu.exon@assays$RNA@counts), ncol = length(miss.cells.exon), sparse = TRUE)
    colnames(miss.cells.exon.mat) <- miss.cells.exon
    rownames(miss.cells.exon.mat) <- rownames(seu.exon@assays$RNA@counts)
    seu.exon@assays$RNA@counts <- cbind2(seu.exon@assays$RNA@counts, miss.cells.exon.mat)

    # sort exon only matrix
    seu.exon@assays$RNA@counts <- seu.exon@assays$RNA@counts[rownames(seu.all), ]
    seu.exon@assays$RNA@counts <- seu.exon@assays$RNA@counts[, colnames(seu.all)]
    
    ## pad the exon only matrix

    # get intronic counts
    seu.all@assays$RNA@counts <- seu.all@assays$RNA@counts - seu.exon@assays$RNA@counts

    # normalize intronic counts
    polya.norm <- polya.tab[rownames(seu.all@assays$RNA@counts), c("polyA")]
    rownames.all.mat <- rownames(seu.all@assays$RNA@counts)
    colnames.all.mat <- colnames(seu.all@assays$RNA@counts)
    seu.all@assays$RNA@counts <- Matrix(apply(seu.all@assays$RNA@counts, MARGIN = 2, function(x) x / (1 + polya.norm)), sparse = TRUE)

    # add the exonic reads back
    seu.all@assays$RNA@counts <- seu.all@assays$RNA@counts + seu.exon@assays$RNA@counts

    # check if NAs
    seu.all@assays$RNA@counts[is.na(seu.all@assays$RNA@counts)] <- 0
    
    # make negatives zeros
    seu.all@assays$RNA@counts[seu.all@assays$RNA@counts<0] <- 0
    
    # copy to data slot
    seu.all@assays$RNA@data <- seu.all@assays$RNA@counts
    
    # get cell metadata
    cell_metadata <- sample_metadata[x, ]
    for (meta in names(cell_metadata)) {
      seu.all[[meta]] <- cell_metadata[[meta]]
    }

    # get map estimate of cell sex
    if (sex.col %in% colnames(cell_metadata)){
      if  (cell_metadata[,sex.col] == "Female_Male") {
        cat(str_interp("Getting MAP estimates of cell sex: ${x}"), sep = "\n")
        seu.all <- get_sex_map(seu.all)
      }
    }
    
    # remove unnecessary columns
    seu.all$orig.ident <- NULL
    seu.all$barcode <- NULL
    seu.all$id <- NULL

    return(seu.all)
  })

  # set names
  seurat_list <- setNames(seurat_list, keepers)

  # remove NULLs
  seurat_list <- seurat_list[lapply(seurat_list, is.null) == FALSE]
  
  # merge objects
  if(length(seurat_list)>1){
    seurat_obj <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
  } else if (length(seurat_list)==1){
    seurat_obj <- seurat_list[[1]]
  } else {
    return(NULL)
  }

  # add qc metrics
  seurat_obj <- comp_qc_metrics(seurat_obj)
  
  # this needs to be done with all samples: issues in postnatal samples
  # # get gene list from Tirosh et al 2015 (G2/M phase)
  # s.genes <- firstup(cc.genes$s.genes)
  # g2m.genes <- firstup(cc.genes$g2m.genes)
  # 
  # # assign cell cycle score (S.Score, G2M.Score)
  # seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  # free up memory
  rm(seurat_list)
  gc()
  
  return(seurat_obj)
}

readin_sn_soupx_sample <- function(sample.name, cellranger.dir, project_name, analysis, anl.spec){
  
  # load matrix
  sc <- load10X(str_interp("${cellranger.dir}/${sample.name}/outs/"))
  
  # load clustering
  if (analysis) {
    # load secondary analysis from first run
    cr.clusts.df <- fread(str_interp("${cellranger.dir}/${sample.name}/outs/analysis/clustering/graphclust/clusters.csv"))
  } else {
    # load secondary analysis from second run
    cr.clusts.df <- fread(str_interp("${cellranger.dir}/${sample.name}_sec/outs/analysis/clustering/graphclust/clusters.csv"))
  }
  
  # species classification (if existant)
  spec_class_file <- str_interp("${cellranger.dir}/${sample.name}/outs/analysis/gem_classification.csv")
  if(file.exists(spec_class_file)){
    
    # read in species file
    cell.spec.tab <- fread(spec_class_file)
    
    # plot histogram
    plt.df <- data.frame(Species = cell.spec.tab$call)
    ggplot(plt.df,aes(x=Species)) + geom_bar()
    plt.name <- str_interp("${cellranger.dir}/${sample.name}/outs/analysis/gem_classification_barplot_${sample.name}.png")
    ggsave(plt.name, width=10, height=10)
    
    # plot scatter
    ggplot(cell.spec.tab,aes(x=GRCh38,y=mm10,color=call)) + geom_point() + 
      labs(title="Number of reads aligning each genome and cell classification")
    plt.name <- str_interp("${cellranger.dir}/${sample.name}/outs/analysis/gem_classification_scatter_${sample.name}.png")
    ggsave(plt.name, width=10, height=10)
    
    # remove cells of wrong species
    sc$metaData <- sc$metaData[cell.spec.tab$barcode[cell.spec.tab$call==anl.spec],]
    sc$toc <- sc$toc[, cell.spec.tab$barcode[cell.spec.tab$call==anl.spec], drop = FALSE]
    
    # filter cluster df
    cr.clusts.df <- cr.clusts.df[cr.clusts.df$Barcode %in% rownames(sc$metaData),]
    
  }
  
  # assign clusters 
  cr.clusts <- cr.clusts.df$Cluster
  names(cr.clusts) <- cr.clusts.df$Barcode
  sc <- setClusters(sc,clusters = cr.clusts)
  
  # autoestimate contaminant fraction
  result <- try({
    sc <- autoEstCont(sc, forceAccept =TRUE, soupQuantile = 0.9, priorRho = 0.20, tfidfMin = 0.75)
  })
  
  if (inherits(result, "try-error")) {
    cat(str_interp("An error occurred trying to estimate the contaminant fraction of ${sample.name} using SoupX:"), sep="\n")
    cat(result[1], sep="\n")
    cat(str_interp("Skipping this sample in integration"), sep="\n")
    return(NULL)
  } else {
    cat(str_interp("We were able to estimate the contaminant fraction of ${sample.name} using SoupX: ${sc$fit$rhoEst}"), sep="\n")
  }
  
  # exclude sample if fraction is larger than 50%
  if (sc$fit$rhoEst>0.5){
    cat(str_interp("The contaminant fraction in ${sample.name} is above 50%. Skipping this sample"), sep="\n")
    return(NULL)
  }
  
  # filter out contaminant cells
  sc$toc <- adjustCounts(sc)
  
  # exonic reads only (keep those with at least 100 features)
  seu <- CreateSeuratObject(counts = sc$toc, min.cells = 3, min.features = 50, project = project_name)
  seu$SampleID <- sample.name
  
  return(seu)
  
}

# reads in visium samples in cnag.proj
read.visium <- function(cnag.proj,bsc.proj.dir,qc=FALSE){
  
  # init list
  st.lst <- c()
  var.feats <- c()
  
  # iterate over projects
  for(proj in cnag.proj){
    
    print(str_interp("Working with project ${proj}..."))
    
    # read in metadata
    proj.meta <- fread(file = str_interp("${bsc.proj.dir}/${proj}/metadata.tsv"))
    
    # project directory
    proj.sr.dir <- str_interp("${bsc.proj.dir}/${proj}/spaceranger/")
    
    # read in a sample
    proj.sr.smp <- list.dirs(proj.sr.dir, full.names = FALSE, recursive = FALSE)
    
    # iterate over samples in project  
    for(smp.nme in proj.sr.smp){
      
      print(str_interp("Reading ${smp.nme}..."))
      
      # load
      st.data <- Load10X_Spatial(data.dir = str_interp("${proj.sr.dir}/${smp.nme}/outs/"), slice = smp.nme)
      
      # add metadata
      st.data@meta.data$sex <- "male"
      st.data@meta.data$project <- proj
      st.data@meta.data$sample <- smp.nme
      st.data@meta.data$stage <- proj.meta$stage[proj.meta$barcode==smp.nme]
      st.data@meta.data$litter <- proj.meta$litter[proj.meta$barcode==smp.nme]
      st.data@meta.data$species <- proj.meta$species[proj.meta$barcode==smp.nme]
      st.data@meta.data$modality <- proj.meta$modality[proj.meta$barcode==smp.nme]
      
      # set key for spot segmentation
      st.data@images[[smp.nme]]@key <- "slice1_"
      
      # segment spots
      st.data <- segment.spots(st.data)
      
      # change slice name
      st.data@images[[smp.nme]]@key <- str_interp("${smp.name}_")
      
      # seurat directory
      proj.seu.dir <- str_interp("${bsc.proj.dir}/${proj}/seurat/")
      dir.create(proj.seu.dir, recursive = TRUE, showWarnings = FALSE)
      
      # normalization with NB model
      st.data <- SCTransform(st.data, assay = "Spatial", return.only.var.genes = FALSE, verbose = TRUE)
      
      # add variable features
      var.feats <- c(var.feats,VariableFeatures(st.data))
      
      if (qc){
        
        # also run standard log normalization for comparison
        st.data <- NormalizeData(st.data, verbose = FALSE, assay = "Spatial")
        
        # plot counts
        p1 <- VlnPlot(st.data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
        p2 <- SpatialFeaturePlot(st.data, features = "nCount_Spatial") + theme(legend.position = "right")
        p <- wrap_plots(p1, p2)
        ggsave(filename = str_interp("${proj.seu.dir}/${smp.nme}_ncount.pdf"), plot = p)
        
        # Computes the correlation of the log normalized data and sctransform residuals with the number of UMIs
        st.data <- GroupCorrelation(st.data, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
        st.data <- GroupCorrelation(st.data, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
        
        # plot
        p1 <- GroupCorrelationPlot(st.data, assay = "Spatial", cor = "nCount_Spatial_cor") +
          ggtitle("Log Normalization") + theme(plot.title = element_text(hjust = 0.5))
        p2 <- GroupCorrelationPlot(st.data, assay = "SCT", cor = "nCount_Spatial_cor") +
          ggtitle("SCTransform Normalization") + theme(plot.title = element_text(hjust = 0.5))
        p <- p1 + p2
        ggsave(filename = str_interp("${proj.seu.dir}/${smp.nme}_corr_ncount_norm.pdf"), plot = p)
      }
      
      # add to list
      st.lst <- c(st.lst,st.data)
    }
  }
  
  return(list(seu=st.lst,vf=var.feats))
}

# create seurat dir
create_seurat_dir <- function(project_dir) {
  # seurat directory
  seu_dir <- paste(project_dir, "seurat", sep = "/")

  # check if seurat dir exists
  if (!dir.exists(seu_dir)) {
    dir.create(seu_dir)
  }
}
#########################################################################################################################################################
# Subsetting Seurat object (used in downstream analysis)
#########################################################################################################################################################
so_filter_stg_annot <- function(so, keep.stages, rm.cell_annot){
  
  # subset according to stage
  so <- subset(so, stage %in% keep.stages)
  
  # remove unwanted levels
  so@meta.data$stage <- factor(so@meta.data$stage,levels=keep.stages)
  
  # remove contaminant cells
  so <- subset(so, cell_annot %in% setdiff(levels(so$cell_annot), rm.cell_annot))
  
  # umap to double check
  # DimPlot(so, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
  #   guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
  
  # remove isolated cells
  umap_coords <- as.data.frame(Embeddings(so, reduction = "umap"))
  umap_coords$dist <- compute_nearest_neighbor_distance(umap_coords)
  iso.cells <- rownames(umap_coords)[umap_coords$dist>quantile(umap_coords$dist, 0.99)]
  so <- so[,!(colnames(so) %in% iso.cells)]
  
  # umap to double check
  # DimPlot(so, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
  #   guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
  
  # remove cells in clusters with low abundance
  clust.cnt.tab <- table(so$cell_annot)
  keep.clust <- names(clust.cnt.tab[clust.cnt.tab>50])
  so <- subset(so, cell_annot %in% keep.clust)
  so@meta.data$cell_annot <- factor(so@meta.data$cell_annot,level=unique(so@meta.data$cell_annot))
  
  # umap to double check
  DimPlot(so, group.by = "cell_annot", label = TRUE) + theme(legend.position = "top") +
    guides(color = guide_legend(ncol = 10, byrow = TRUE, override.aes = list(size = 1)))
  
  gc()
  
  return(so)
  
}
#########################################################################################################################################################
# QC
#########################################################################################################################################################

# add qc metrics
comp_qc_metrics <- function(seu.obj){
  
  seu.obj[["percent.mt"]] <- PercentageFeatureSet(seu.obj, pattern = "^[Mm][Tt]") # mitochondrial
  seu.obj[["percent.hb"]] <- PercentageFeatureSet(seu.obj, pattern = "^[Hh][Bb]") # hemoglobin
  seu.obj[["percent.rp"]] <- PercentageFeatureSet(seu.obj, pattern = "^R[Pp]") # ribosomal
  seu.obj[["nCount_RNA"]] <- colSums(x = seu.obj, slot = "counts", na.rm = TRUE)
  seu.obj[["nFeature_RNA"]] <- colSums(x = GetAssayData(object = seu.obj, slot = "counts") > 0, na.rm = TRUE)
  seu.obj[["log10GenesPerUMI"]] <- log10(seu.obj$nFeature_RNA) / log10(seu.obj$nCount_RNA)
  
  return(seu.obj)
}

# doublet detection
filter_out_doublets <- function(seu, doublet.rate = 0.05) {
  
  # preprocess for doublet finder (https://rpubs.com/kenneditodd/doublet_finder_example)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, verbose = F)
  seu <- ScaleData(seu, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = F)
  seu <- RunPCA(seu, verbose = F, npcs = 20)

  # # pK identification (no ground-truth)
  # sweep.list <- paramSweep_v3(seu, PCs = 1:20)
  # sweep.stats <- summarizeSweep(sweep.list)
  # bcmvn <- find.pK(sweep.stats)
  #
  # # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  # bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  # optimal.pk <- bcmvn.max$pK
  # optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

  # find doublets (expect 5% since we recover 6k)
  # see https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets
  nExp <- round(ncol(seu) * doublet.rate)
  seu <- doubletFinder_v3(seu, pK = doublet.rate, nExp = nExp, PCs = 1:10)

  # get column name for doubletfinder label
  colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))] <- "doubletfinder"
  colnames(seu@meta.data)[grepl("pANN_", colnames(seu@meta.data))] <- "doubletfinder_pann"

  # remove doublets
  return(subset(seu, doubletfinder == "Singlet"))
}

# QC plots
qc_plots <- function(project_dir, seurat_obj, suffix = "") {
  # seurat directory
  seu_dir <- paste(project_dir, "seurat", sep = "/")

  # violin plot
  file_name <- paste(seu_dir, "/qcplot_nfeature_ncount_percent", suffix, ".png", sep = "")
  png(filename = file_name, width = 1200, height = 1200)
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2)
  dev.off()

  # pull out metadata
  metadata <- seurat_obj@meta.data

  # Rename columns
  metadata <- metadata %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

  # Visualize the number of cell counts per sample
  metadata %>%
    ggplot(aes(x = SampleID, fill = SampleID)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells")
  ggsave(stringr::str_interp("${seu_dir}/qcplot_cells_per_sample${suffix}.pdf"), height = 10, width = 15)

  # plot UMIs per cell density
  metadata %>%
    ggplot(aes(x = nUMI)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_umi_density${suffix}.pdf"), height = 15, width = 15)

  # Visualize the correlation between genes detected and number of UMIs
  # and determine whether strong presence of cells with low numbers of genes/UMIs
  metadata %>%
    ggplot(aes(x = nUMI, y = nGene, color = percent.mt)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    # geom_vline(xintercept = 500) +
    # geom_hline(yintercept = 250) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_scatter_umi_nGenes_mt${suffix}.pdf"), height = 15, width = 15)

  # Visualize the correlation between genes detected and number of UMIs
  # and determine whether strong presence of cells with low numbers of genes/UMIs
  metadata %>%
    ggplot(aes(x = nUMI, y = nGene, color = percent.hb)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    # geom_vline(xintercept = 500) +
    # geom_hline(yintercept = 250) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_scatter_umi_nGenes_hb${suffix}.pdf"), height = 15, width = 15)

  # Visualize the correlation between genes detected and number of UMIs
  # and determine whether strong presence of cells with low numbers of genes/UMIs
  metadata %>%
    ggplot(aes(x = nUMI, y = nGene, color = percent.rp)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    # geom_vline(xintercept = 500) +
    # geom_hline(yintercept = 250) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_scatter_umi_nGenes_rp${suffix}.pdf"), height = 15, width = 15)

  # Visualize the distribution of mitochondrial gene expression detected per cell
  # This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying
  # cells. We define poor quality samples for mitochondrial counts as cells which surpass the
  # 20% mitochondrial ratio mark, unless of course you are expecting this in your sample.
  metadata %>%
    ggplot(aes(x = percent.mt)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 20) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_mitoratio${suffix}.pdf"), height = 15, width = 15)

  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  metadata %>%
    ggplot(aes(x = log10GenesPerUMI)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    facet_wrap(~SampleID, nrow = 3)
  ggsave(stringr::str_interp("${seu_dir}/qcplot_complexity${suffix}.pdf"), height = 15, width = 15)
}

# subsets the seurat object according to the number of RNA genes, the number of RNA molecules, and the mt and hb percentage
qc_subset <- function(seu.obj, nFeature_RNA_min = 500, nCount_RNA_min=0, nCount_RNA_max = 10000, percent.mt_max = 10, percent.hb_max = 1.5) {
  
  # subset
  seu.obj <- subset(seu.obj, nFeature_RNA >= nFeature_RNA_min)
  seu.obj <- subset(seu.obj, nCount_RNA >= nCount_RNA_min)
  seu.obj <- subset(seu.obj, nCount_RNA <= nCount_RNA_max)
  seu.obj <- subset(seu.obj, percent.mt <= percent.mt_max)
  seu.obj <- subset(seu.obj, percent.hb <= percent.hb_max)

  return(seu.obj)
}

# coverts to lowercase first and then capitalizes first letter
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#########################################################################################################################################################
# 3D UMAPs
#########################################################################################################################################################

plt_3d_umap_clust <- function(so, res, outdir, pref="umap_3d_seurat_clusters"){
  
  # Visualize what headings are called so that you can extract them to form a dataframe
  # Embeddings(object = so, reduction = "umap")
  
  # Prepare a dataframe for cell plotting
  plot.data <- FetchData(object = so, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "cell_annot"))
  
  # Make a column of row name identities (these will be your cell/barcode names)
  plot.data$label <- paste(rownames(plot.data))
  
  # Plot your data, in this example my Seurat object had 21 clusters (0-20)
  fig <- plot_ly(
    data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~cell_annot, 
    colors = c("lightseagreen","gray50","darkgreen","red4","red","turquoise4","black","yellow4","royalblue1",
                              "lightcyan3","peachpuff3","khaki3","gray20","orange2","royalblue4","yellow3","gray80","darkorchid1",
                              "lawngreen","plum2","darkmagenta"),
                              type = "scatter3d", mode = "markers", 
    marker = list(size = 2, width=2), # controls size of points
    text=~label, #This is that extra column we made earlier for which we will use for cell ID
    hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
  
  # xaxis
  axx <- list(
    nticks = 4,
    range = c(-10,10) #select range of xaxis
  )
  
  # yaxis
  axy <- list(
    nticks = 4,
    range = c(-10,10) #select range of yaxis
  )
  
  #zaxis
  axz <- list(
    nticks = 4,
    range = c(-10,10) #select range of zaxis
  )
  
  fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
  fn <- str_interp("${outdir}/${pref}_r${res}.html")
  htmlwidgets::saveWidget(widget = fig_cube, file = fn, selfcontained = TRUE, libdir = NULL)
  
}

plt_3d_umap_gene <- function(so, outdir, gene){
  
  # Visualize what headings are called so that you can extract them to form a dataframe
  # Embeddings(object = so, reduction = "umap")
  
  # create a dataframe
  plot.data <- FetchData(object = so, vars = c("UMAP_1", "UMAP_2", "UMAP_3", str_interp("${gene}")), slot = 'data')
  colnames(plot.data)[4] <- "Expression"
  
  # Add the label column, so that now the column has 'cellname-its expression value'
  plot.data$label <- paste(rownames(plot.data)," - ", plot.data$Expression, sep="")
  
  # Plot your data, in this example my Seurat object had 21 clusters (0-20)
  fig <- plot_ly(
    data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Expression, 
    opacity = .5,colors = c('gray', 'red'), type = "scatter3d", mode = "markers", 
    marker = list(size = 5, width=2), text=~label,hoverinfo="text"
  )
  
  # xaxis
  axx <- list(
    nticks = 4,
    range = c(-10,10) #select range of xaxis
  )
  
  # yaxis
  axy <- list(
    nticks = 4,
    range = c(-10,10) #select range of yaxis
  )
  
  #zaxis
  axz <- list(
    nticks = 4,
    range = c(-10,10) #select range of zaxis
  )
  
  fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
  fn <- str_interp("${outdir}/umap_3d_${gene}.html")
  htmlwidgets::saveWidget(widget = fig_cube, file = fn, selfcontained = TRUE, libdir = NULL)
  
}

#########################################################################################################################################################
# Transcriptomic analysis
#########################################################################################################################################################

# iterates over clusters and identifies marker genes
find_mrkrs <- function(seu_obj) {
  
  # Count the number of cells in each cluster
  cell_counts <- table(Idents(seu_obj))
  
  # find conserved markers for each cluster
  cvd_mrk <- list()
  for (cls in levels(Idents(seu_obj))) {
    if(cell_counts[[cls]]>3){
      cat(str_interp("Finding markers in cluster ${cls}..."))
      cvd_mrk[[cls]] <- FindMarkers(seu_obj, ident.1 = cls, verbose = TRUE)
    } else {
      cat(str_interp("Cluster ${cls} has less than 3 cells. Skipping..."))
    }
  }
  
  # return list of conserved markers
  return(cvd_mrk)
  
}

# iterates over clusters and identifies conserved marker genes grouping by genotype
find_cons_mrkrs <- function(seu_obj) {
  
  # find conserved markers for each cluster
  cvd_mrk <- list()
  for (cls in levels(Idents(seu_obj))) {
    print(cls)
    cvd_mrk[[cls]] <- FindConservedMarkers(seu_obj, ident.1 = cls, grouping.var = "genotype", verbose = TRUE)
  }

  # return list of conserved markers
  return(cvd_mrk)
  
}

# generates marker list for seurat clusters
get_clust_mark <- function(seu.obj, mark.clust.file, mark.dir, bmrt.tab.file, int.ver, res, n = 1000, plotting=FALSE, p.filt=TRUE, format="txt") {
  
  # check if marker gene analysis is done
  if (file.exists(mark.clust.file)) {
    
    cat(stringr::str_interp("Cluster marker file does exist. Reading..."),sep="\n")
    
    # read in integrated file
    cvd_mrk.list <- readRDS(mark.clust.file)
    
  } else {
    
    cat(stringr::str_interp("Cluster marker file does not exist. Generating..."),sep="\n")
    
    # find markers
    cvd_mrk.list <- find_mrkrs(seu.obj)
    
    # save marker gene list
    saveRDS(cvd_mrk.list, file = mark.clust.file)
    
  }
  
  # Extract top markers per cluster
  top.df <- data.frame()
  for (clust.id in names(cvd_mrk.list)) {
    
    cat(str_interp("Cluster: ${clust.id}"),sep="\n")
    
    if (p.filt){
      
      # filter p's and get top n markers in abs value
      top.clust <- cvd_mrk.list[[clust.id]] %>%
        mutate(abs_log2FC = abs(avg_log2FC)) %>%
        filter(p_val_adj <= 0.1) %>%
        top_n(n = n, wt = abs_log2FC)
      
    } else {
      
      # compute absolute value
      top.clust <- cvd_mrk.list[[clust.id]] %>%
        mutate(abs_log2FC = abs(avg_log2FC))
      
    }

    # assign cluster name
    top.clust$Cluster <- clust.id
    
    # create gene column
    top.clust <- tibble::rownames_to_column(top.clust, "Gene")
    
    # sort based on average effect size
    top.clust <- top.clust[order(top.clust$abs_log2FC,decreasing = TRUE),]

    if (plotting){
      
      # plot dir
      plt.dir <- stringr::str_interp("${mark.dir}/plots/")
      dir.create(plt.dir, recursive = TRUE, showWarnings = FALSE)
      
      # function that returns name of plot
      plt.name <- function(x) {
        stringr::str_interp("${plt.dir}/umap_top${n}_seurat_markers_clust${clust.id}_${x}_${int.ver}_res${res}.png")
      }
  
      # function that produces individual featureplots
      lapply(top.clust$Gene[1:10], FUN = function(x){
        p <- FeaturePlot(seu.obj, features = x)
        ggsave(plt.name(x), plot = p, width = 20, height = 20)
      })
      
    }

    # append to table
    top.df <- plyr::rbind.fill(top.df, top.clust)

  }

  # read in biomart table
  bmrt.tab <- read.table(bmrt.tab.file, sep = "\t", header = TRUE, quote = "\"")
  bmrt.tab <- bmrt.tab[bmrt.tab$external_gene_name != "", ]
  bmrt.tab <- bmrt.tab[, c("external_gene_name", "gene_biotype", "description")]
  bmrt.tab <- bmrt.tab[!duplicated(bmrt.tab$external_gene_name), ]
  rownames(bmrt.tab) <- bmrt.tab$external_gene_name
  
  # breakdown of gene class
  top.df$gene_type <- bmrt.tab[top.df$Gene, c("gene_biotype")]
  top.df$gene_description <- bmrt.tab[top.df$Gene, c("description")]
  
  # obtain cluster specific markers
  top.df$Specific <- "No"
  for(c in unique(top.df$Cluster)){
    
    # get cluster pos and neg markers
    c.pos <- top.df$Gene[(top.df$Cluster==c)&(top.df$avg_log2FC>0)]
    c.neg <- top.df$Gene[(top.df$Cluster==c)&(top.df$avg_log2FC<0)]
    
    # get background pos and neg markers
    b.pos <- top.df$Gene[(top.df$Cluster!=c)&(top.df$avg_log2FC>0)]
    b.neg <- top.df$Gene[(top.df$Cluster!=c)&(top.df$avg_log2FC<0)]
    
    # take difference between sets
    spc.pos.mark <- setdiff(c.pos,b.pos)
    spc.neg.mark <- setdiff(c.neg,b.neg)
    
    # get entries in table
    top.df$Specific[(top.df$Cluster==c)&(top.df$Gene %in% spc.pos.mark)] <- "Positive"
    top.df$Specific[(top.df$Cluster==c)&(top.df$Gene %in% spc.neg.mark)] <- "Negative"
    
  }
  
  # store top markers per cluster
  if (p.filt){
    fn <- stringr::str_interp("${mark.dir}/table_top${n}_seurat_markers_clusters_${int.ver}_res${res}")
  } else {
    fn <- stringr::str_interp("${mark.dir}/table_lfc_seurat_markers_clusters_${int.ver}_res${res}")
  }
  
  # write table and excel
  if(format=="txt"){
    fn <- stringr::str_interp("${fn}.txt")
    write.table(top.df, fn, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (format=="xlsx"){
    fn <- stringr::str_interp("${fn}.xlsx")
    write.xlsx(top.df, fn, col.names = TRUE, row.names = FALSE, append = FALSE)
  } else {
    print("Incorrect format for storing table. Choose between: txt or xlsx")
  }

  # return top markers
  return(top.df)
}

# generates conserved marker list for seurat clusters
get_cons_clust_mark <- function(seu.obj, mark.clust.file, mark.dir, bmrt.tab.file, int.ver, res, n = 100, plotting=TRUE) {
  
  # check if marker gene analysis is done
  if (!file.exists(mark.clust.file)) {
    
    cat(stringr::str_interp("Cluster marker file for ${int.ver} with res ${res} does not exist. Generating..."),sep="\n")
    
    # find conserved markers
    cvd_mrk.list <- find_cons_mrkrs(seu.obj)
    
    # save marker gene list
    saveRDS(cvd_mrk.list, file = mark.clust.file)
    
  } else {
    
    cat(stringr::str_interp("Cluster marker file for ${int.ver} with res ${res} does exist. Reading..."),sep="\n")
    
    # read in integrated file
    cvd_mrk.list <- readRDS(mark.clust.file)
    
  }
  
  # Extract top markers per cluster
  top.df <- data.frame()
  for (clust.id in names(cvd_mrk.list)) {
    
    cat(str_interp("Cluster: ${clust.id}"),sep="\n")
    
    hd_fc <- "HD_avg_log2FC" %in% colnames(cvd_mrk.list[[clust.id]])
    wt_fc <- "WT_avg_log2FC" %in% colnames(cvd_mrk.list[[clust.id]])
    
    # get top n markers (average log FC)
    if (hd_fc & wt_fc) {
      top.clust <- cvd_mrk.list[[clust.id]] %>%
        mutate(avg_fc = (HD_avg_log2FC + WT_avg_log2FC) / 2) %>%
        mutate(avg_fc_abs = abs(HD_avg_log2FC + WT_avg_log2FC) / 2) %>%
        top_n(n = n, wt = avg_fc_abs)
    } else if (!hd_fc & wt_fc) {
      top.clust <- cvd_mrk.list[[clust.id]] %>%
        mutate(avg_fc = WT_avg_log2FC) %>%
        mutate(avg_fc_abs = abs(WT_avg_log2FC)) %>%
        top_n(n = n, wt = avg_fc_abs)
    } else if (hd_fc & !wt_fc) {
      top.clust <- cvd_mrk.list[[clust.id]] %>%
        mutate(avg_fc = HD_avg_log2FC) %>%
        mutate(avg_fc_abs = abs(HD_avg_log2FC)) %>%
        top_n(n = n, wt = avg_fc_abs)
    } else {
      next
    }
    
    # assign cluster name
    top.clust$Cluster <- clust.id
    
    # create gene column
    top.clust <- tibble::rownames_to_column(top.clust, "Gene")
    
    # sort based on average effect size
    top.clust <- top.clust[order(top.clust$avg_fc_abs,decreasing = TRUE),]
    
    if (plotting){
      
      # plot dir
      plt.dir <- stringr::str_interp("${mark.dir}/plots/")
      dir.create(plt.dir, recursive = TRUE, showWarnings = FALSE)
      
      # function that returns name of plot
      plt.name <- function(x) {
        stringr::str_interp("${plt.dir}/umap_top${n}_seurat_markers_clust${clust.id}_${x}_${int.ver}_res${res}.png")
      }
      
      # function that produces individual featureplots
      lapply(top.clust$Gene[1:10], FUN = function(x){
        p <- FeaturePlot(seu.obj, features = x)
        ggsave(plt.name(x), plot = p, width = 20, height = 20)
      })
      
    }
    
    # append to table
    top.df <- plyr::rbind.fill(top.df, top.clust)
    
  }
  
  # read in biomart table
  bmrt.tab <- read.table(bmrt.tab.file, sep = "\t", header = TRUE, quote = "\"")
  bmrt.tab <- bmrt.tab[bmrt.tab$external_gene_name != "", ]
  bmrt.tab <- bmrt.tab[, c("external_gene_name", "gene_biotype", "description")]
  bmrt.tab <- bmrt.tab[!duplicated(bmrt.tab$external_gene_name), ]
  rownames(bmrt.tab) <- bmrt.tab$external_gene_name
  
  # breakdown of gene class
  top.df$gene_type <- bmrt.tab[top.df$Gene, c("gene_biotype")]
  top.df$gene_description <- bmrt.tab[top.df$Gene, c("description")]
  
  # store top markers per cluster
  fn <- stringr::str_interp("${mark.dir}/table_top${n}_seurat_markers_clusters_${int.ver}_res${res}.txt")
  write.table(top.df, fn, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # return top markers
  return(top.df)
}

# generates marker list for sc-Embryo, sc-Postnatal, and sn-Postnatal
get_emb_vs_post_mark <- function(seu.obj, mark.sn_vs_sc.file, mark.dir, bmrt.tab.file, int.ver, rev.ver, n = 100) {
  # find conserved markers
  if (!file.exists(mark.sn_vs_sc.file)) {
    print(stringr::str_interp("Embryo vs postnatal marker file for ${int.ver}.${rev.ver} does not exist. Generating..."))

    # define category
    emb.proj <- c("09", "16")
    sn.post.proj <- c("15", "17")
    sc.post.proj <- c("Anderson")
    categ.list <- list(Embryonic = emb.proj, SnPostnatal = sn.post.proj, ScPostnatal = sc.post.proj)

    # get counts per category
    tab.list <- lapply(categ.list, FUN = function(x) {
      tab <- table(seu.obj@meta.data$seurat_clusters[seu.obj@meta.data$project %in% x])
      names(tab[tab > 500])
    })

    # remove empty entries in list
    tab.list <- tab.list[lapply(tab.list,length)>0]
    
    # find conserved markers for each cluster
    cvd_mrk.list <- lapply(tab.list, FUN = function(x) {
      FindConservedMarkers(seu.obj, ident.1 = x, grouping.var = "genotype", verbose = FALSE)
    })

    # store
    saveRDS(cvd_mrk.list, file = mark.sn_vs_sc.file)
  } else {
    print(stringr::str_interp("Embryo vs postnadal marker file for ${int.ver}.${rev.ver} does exist. Reading..."))
    cvd_mrk.list <- readRDS(mark.sn_vs_sc.file)
  }

  # plot dir
  plt.dir <- stringr::str_interp("${mark.dir}/plots/")
  dir.create(plt.dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract top markers per cluster
  top.df <- data.frame()
  for (categ.id in names(cvd_mrk.list)) {
    hd_fc <- "HD_avg_log2FC" %in% colnames(cvd_mrk.list[[categ.id]])
    wt_fc <- "WT_avg_log2FC" %in% colnames(cvd_mrk.list[[categ.id]])

    # get top n markers (average log FC)
    if (hd_fc & wt_fc) {
      top.clust <- cvd_mrk.list[[categ.id]] %>%
        mutate(avg_fc = (HD_avg_log2FC + WT_avg_log2FC) / 2) %>%
        top_n(n = n, wt = avg_fc)
    } else if (!hd_fc & wt_fc) {
      top.clust <- cvd_mrk.list[[categ.id]] %>%
        mutate(avg_fc = WT_avg_log2FC) %>%
        top_n(n = n, wt = avg_fc)
    } else if (hd_fc & !wt_fc) {
      top.clust <- cvd_mrk.list[[categ.id]] %>%
        mutate(avg_fc = HD_avg_log2FC) %>%
        top_n(n = n, wt = avg_fc)
    } else {
      next
    }

    # assign cluster name
    top.clust$Category <- categ.id
    
    # create gene column
    top.clust <- tibble::rownames_to_column(top.clust, "Gene")
    
    # sort based on average effect size
    top.clust <- top.clust[order(top.clust$avg_fc,decreasing = TRUE),]
    
    # function that returns name of plot
    plt.name <- function(x) {
      stringr::str_interp("${plt.dir}/umap_top${n}_seurat_markers_${categ.id}_${x}_${int.ver}.${rev.ver}.png")
    }

    # function that produces individual featureplots
    lapply(top.clust$Gene[1:10], FUN = function(x){
      p <- FeaturePlot(seu.obj, features = x)
      ggsave(plt.name(x), plot = p, width = 20, height = 20)
    })
    
    # append to table
    top.df <- plyr::rbind.fill(top.df, top.clust)
  }

  # read in biomart table
  bmrt.tab <- read.table(bmrt.tab.file, sep = "\t", header = TRUE, quote = "\"")
  bmrt.tab <- bmrt.tab[bmrt.tab$external_gene_name != "", ]
  bmrt.tab <- bmrt.tab[, c("external_gene_name", "gene_biotype", "description")]
  bmrt.tab <- bmrt.tab[!duplicated(bmrt.tab$external_gene_name), ]
  rownames(bmrt.tab) <- bmrt.tab$external_gene_name

  # breakdown of gene class
  top.df$gene_type <- bmrt.tab[top.df$Gene, c("gene_biotype")]
  top.df$gene_description <- bmrt.tab[top.df$Gene, c("description")]

  # store table
  fn <- stringr::str_interp("${mark.dir}/table_top${n}_seurat_markers_groups_${int.ver}.${rev.ver}.txt")
  write.table(top.df, file = fn, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  # return table
  return(top.df)
}

# binomial glm analysis for differential gene expression TFs
glm.bin.diff.exp <- function(seu.obj,gene.name,out.dir,stg.levs=c("E12.5","E14.5","E16.5"),int.ver,res.ver){
  
  # add metadata column using raw counts (counts slot)
  seu.obj@meta.data$gene_pos <- seu.obj@assays$RNA@counts[gene.name,] > 0
  
  # statistical test data
  df.test <- seu.obj@meta.data[c("genotype","stage","cell_annot","gene_pos")]
  df.test$genotype <- factor(df.test$genotype,levels=c("WT","HD"))
  df.test$stage <- factor(df.test$stage,levels=stg.levs)
  
  # glm
  gene.glm <- glm(gene_pos ~ genotype + stage + cell_annot + genotype*stage + 
                    genotype*cell_annot + stage*cell_annot, df.test, family=binomial)
  summ <- summary(gene.glm)
  summ.df <- tibble::rownames_to_column(as.data.frame(summ$coefficients), "Coefficient")
  write.table(
    summ.df,file = str_interp("${out.dir}/glm_${gene.name}_cellannot_${int.ver}_r${res.ver}.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # store OR
  gene.or.glm <- get.or.glm.stage.cellannot(df.test,summ.df)
  write.table(
    gene.or.glm,file = str_interp("${out.dir}/glm_${gene.name}_hd_oddsratio_cellannot_${int.ver}_r${res.ver}.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # plots
  FeaturePlot(seu.obj, feature=gene.name, split.by="genotype")
  ggsave(str_interp("${out.dir}/umap_${gene.name}_exp_${int.ver}.png"), height = 10, width = 10)
  FeaturePlot(seu.obj, feature="gene_pos", split.by="genotype")
  ggsave(str_interp("${out.dir}/umap_${gene.name}_pos_${int.ver}.png"), height = 10, width = 10)
  
  # return OR table
  return(gene.or.glm)
  
} 
#########################################################################################################################################################
# GSEA
#########################################################################################################################################################
# function that runs gene set ontology analysis
run_gsea <- function(gene_list, GO_file, pval=0.05) {
  
  set.seed(54321)
  
  # create GO object
  myGO = fgsea::gmtPathways(GO_file)
  
  # run gsea
  # Run GSEA
  fgRes <- fgsea::fgseaMultilevel(pathways = myGO,stats = gene_list, minSize=15, maxSize=400,scoreType = "std") %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes), pathways = myGO, stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 20), tail(fgRes, n = 20 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"), c("Up-regulated", "Down-regulated"))
  
  g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",title=header)
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}
#########################################################################################################################################################
# Plots
#########################################################################################################################################################

plt.feat.umaps <- function(seu.obj, feat.lst, grp.dir, ind.dir, int.ver) {
  for (x in names(feat.lst)) {
    # grouped
    p <- FeaturePlot(seu.obj, features = feat.lst[[x]]) +
      patchwork::plot_layout(guides = "collect")
    ggsave(stringr::str_interp("${grp.dir}/umap_${x}_${int.ver}.png"), width = 15, height = 15)
    Sys.sleep(0.2)

    # individual
    for (gene in feat.lst[[x]]) {
      if (gene %in% rownames(seu.obj@assays$RNA)) {
        p <- FeaturePlot(seu.obj, features = gene)
        ggsave(stringr::str_interp("${ind.dir}/umap_${x}_${gene}_${int.ver}.png"), width = 15, height = 15)
        Sys.sleep(0.2)
      }
    }
  }
}

#########################################################################################################################################################
# In silico sex deconvolution (only for mouse)
#########################################################################################################################################################

# adds X metadata to seurat object
comp_x <- function(x.seu) {
  
  ## chrY (X1)

  # read in ensembl table
  ensembl.table <- fread("/gpfs/projects/bsc83/Data/Creatio/biomart/mmusculus_ensembl_table.tsv")

  # get chrY genes
  chrY.gene <- ensembl.table$external_gene_name[ensembl.table$chromosome_name == "Y"]

  # overlap in genes
  chrY.list <- chrY.gene[chrY.gene %in% rownames(x.seu@assays$RNA@counts)]

  # remove genes in par region (Kasahara et al 2023)
  par_list <- c("Mid1", "Sts", "Nlgn4", "Asmt")
  chrY.list <- chrY.list[!chrY.list %in% par_list]

  # for each cell calculate percentage of chrY content
  x.seu$x1 <- colSums(x.seu@assays$RNA@counts[rownames(x.seu@assays$RNA@counts) %in% chrY.gene, ])

  ## chrX inactivation (X2)

  # define chrX inactivation genes (only ones that make a difference)
  chrX_inac <- c("Xist", "Tsix")

  # for each cell calculate counts of chrX inactivation genes
  if (sum(rownames(x.seu@assays$RNA@counts) %in% chrX_inac) == 2) {
    # if both genes are in the matrix
    x.seu$x2 <- colSums(x.seu@assays$RNA@counts[rownames(x.seu@assays$RNA@counts) %in% chrX_inac, ])
  } else if (sum(rownames(x.seu@assays$RNA@counts) %in% chrX_inac) == 1) {
    # if only one is
    x.seu$x2 <- x.seu@assays$RNA@counts[rownames(x.seu@assays$RNA@counts) %in% chrX_inac, ]
  } else {
    # if none is
    x.seu$x2 <- 0
  }

  ## rest

  # for each cell calculate percentage of rest of genes
  x.seu$x3 <- colSums(x.seu@assays$RNA@counts) - x.seu$x1 - x.seu$x2

  # return
  return(x.seu)
}

# em algorithm given
emAlg <- function(x, M = 100, k = 2) {
  # Parameters:
  # x: count matrix
  # M: total number of optimization attempts EM algorithm

  # run multiple instances of the EM algorithm
  maxPost <- c()
  maxTheta <- c()
  maxLambda <- c()
  maxloglik <- -100000
  for (i in 1:M) {
    print(stringr::str_interp("i=${i}"))

    # run algorithm
    tryCatch(
      {
        em.out <- mixtools::multmixEM(x, k = k, maxit = 1000, verb = FALSE, epsilon = 1e-8)

        # update parameters if larger loglikelihood
        if (em.out$loglik > maxloglik) {
          maxTheta <- em.out$theta
          maxLambda <- em.out$lambda
          maxloglik <- em.out$loglik
          maxPost <- em.out$posterior
        }
      },
      warning = function(war) {
        # nothing
      },
      error = function(err) {
        # nothing
      }
    )
  }

  return(list(maxPost = maxPost, maxTheta = maxTheta, maxLambda = maxLambda, maxloglik = maxloglik))
}

# computes map estimate
map_est <- function(emAlgOut, thresh = 0.5) {
  
  # take ratio and difference between parameters
  r1 <- emAlgOut$maxTheta[1, 1] / emAlgOut$maxTheta[1, 2]
  r2 <- emAlgOut$maxTheta[2, 1] / emAlgOut$maxTheta[2, 2]
  d1 <- emAlgOut$maxTheta[1, 1] - emAlgOut$maxTheta[1, 2]
  d2 <- emAlgOut$maxTheta[2, 1] - emAlgOut$maxTheta[2, 2]

  # check which are the XY parameters
  if (r1 >= 1 & r2 >= 1) {
    # if no component has XIC counts then assign male
    return(rep("Male", nrow(emAlgOut$maxPost)))
  } else if (d1 > d2) {
    # means Y expression is more probable
    pZ <- emAlgOut$maxPost[, 1]
  } else {
    # means XIC expression is more probable
    pZ <- emAlgOut$maxPost[, 2]
  }

  return(ifelse(pZ >= thresh, "Male", "Female"))
}

# returns df with X metadata
get_x_from_seuobj <- function(x.seu) {
  return(as.matrix(data.frame(x1 = x.seu$x1, x2 = x.seu$x2, x3 = x.seu$x3), ncol = 3))
}

# filter out XY cells
get_sex_map <- function(seu) {
  
  # add X metadata
  seu <- comp_x(seu)

  # create matrix
  Xmat <- matrix(get_x_from_seuobj(seu), ncol = 3)

  # fit model
  emAlgOut <- emAlg(Xmat, M = 500)

  # get posterior for males
  seu[["sex"]] <- map_est(emAlgOut, thresh = 0.9)

  # drop count columns
  seu$x1 <- NULL
  seu$x2 <- NULL
  seu$x3 <- NULL

  # subset based on estimated sex
  return(seu)
  
}

#####################################################################################################################
# Composition analysis
#####################################################################################################################

compAnalysisCellType <- function(sobj, x) {
  
  # define binary
  sobj@meta.data$bin.type <- sobj@meta.data$cell_annot == x

  # get table
  comp.tab <- as.data.frame(sobj@meta.data[c("genotype", "stage", "bin.type")])
  
  # relevel to establish as a reference HD
  comp.tab <- within(comp.tab, genotype <- relevel(factor(genotype), ref = "WT"))
  comp.tab <- within(comp.tab, stage <- relevel(factor(stage), ref = "E12.5"))

  # fit binomial glm (unstable)
  # bin.glm <- glm(formula = bin.type ~ genotype + stage + genotype * stage, data = comp.tab, family = binomial)
  
  # summarize
  tot.count <- as.data.frame(table(comp.tab[,c("genotype","stage")]))
  colnames(tot.count) <- c("genotype", "stage", "total")
  succ.count <- as.data.frame(table(comp.tab))
  succ.count <- succ.count[succ.count$bin.type=="TRUE",c("genotype", "stage","Freq")]
  colnames(succ.count) <- c("genotype", "stage", "success")
  glm.tab <- merge(succ.count,tot.count)
  
  # add pseudocounts for stability
  glm.tab$success <- glm.tab$success+1
  glm.tab$total <- glm.tab$total+1
  
  # fit binomial regression
  bin.glm <- glm(formula = cbind(success, total - success) ~ genotype + stage + genotype * stage, data = glm.tab, family = binomial)

  # extract coefficients
  summ <- summary(bin.glm)
  summ.df <- tibble::rownames_to_column(as.data.frame(summ$coefficients), "Coefficient")

  # add cell type column
  summ.df$CellType <- x

  # return dataframe
  return(summ.df)
}

# calculate odds ratio in glm ~ genotype * stage
get_odds_ratio_mat_glm <- function(model, stages) {
  
  # Initialize an empty table to store the results
  or.mat <- matrix(NA, nrow = 1, ncol = length(stages))
  colnames(or.mat) <- stages
  
  # get coefficient for genotype beta
  coef_gt <- model[model$Coefficient=="genotypeHD", "Estimate"]
  p_value_gt <- model[model$Coefficient=="genotypeHD", "Pr(>|z|)"]
  
  # Calculate odds ratios for each combination of genotype and stage levels
  for (j in 1:length(stages)) {
    
    if (stages[j] == "E12.5") {
      
      # For E12.5 stage, use the coefficient for HD genotype only
      or.mat[1, j] <- exp((p_value_gt <= 0.05) * coef_gt)
      
    } else {
      
      coeff.name <- paste0("genotypeHD:stage", stages[j])
      
      if (coeff.name %in% model$Coefficient){
        
        # Get the coefficient and p-value for the interaction term
        coef_gt_stage <- model[model$Coefficient==coeff.name, "Estimate"]
        p_value_gt_stage <- model[model$Coefficient==coeff.name,"Pr(>|z|)"]
        
        # Add odds ratio to the table if p-value is significant
        or.mat[1, j] <- exp((p_value_gt <= 0.05) * coef_gt +
                              (p_value_gt_stage<=0.05) * coef_gt_stage)
        
      }
      
    }
  }
  
  # make into dataframe
  or.df <- as.data.frame(or.mat)
  
  return(or.df)
}

# calculate probabilities for wt
get_probs_day_glm <- function(model, stages) {
  
  # Initialize an empty table to store the results
  prob.mat <- matrix(NA, nrow = 1, ncol = length(stages))
  colnames(prob.mat) <- stages
  
  # get coefficient for genotype beta
  coef_0 <- model[model$Coefficient=="(Intercept)", "Estimate"]
  pval_0 <- model[model$Coefficient=="(Intercept)", "Pr(>|z|)"]
  coef_1 <- model[model$Coefficient=="stageE14.5", "Estimate"]
  pval_1 <- model[model$Coefficient=="stageE14.5", "Pr(>|z|)"]
  coef_2 <- model[model$Coefficient=="stageE16.5", "Estimate"]
  pval_2 <- model[model$Coefficient=="stageE16.5", "Pr(>|z|)"]
  coef_3 <- model[model$Coefficient=="stage1 month", "Estimate"]
  pval_3 <- model[model$Coefficient=="stage1 month", "Pr(>|z|)"]
  coef_4 <- model[model$Coefficient=="stage5 months", "Estimate"]
  pval_4 <- model[model$Coefficient=="stage5 months", "Pr(>|z|)"]
  
  # set to zero non significant coefficients
  coef_0 <- (pval_0<=0.05) * coef_0
  coef_1 <- (pval_1<=0.05) * coef_1
  coef_2 <- (pval_2<=0.05) * coef_2
  coef_3 <- (pval_3<=0.05) * coef_3
  coef_4 <- (pval_4<=0.05) * coef_4
  
  # Calculate probabilities (inverse logit)
  prob.mat[1,1] <- exp(coef_0)/(exp(coef_0)+1)
  prob.mat[1,2] <- exp(coef_0+coef_1)/(exp(coef_0+coef_1)+1)
  prob.mat[1,3] <- exp(coef_0+coef_2)/(exp(coef_0+coef_2)+1)
  prob.mat[1,4] <- exp(coef_0+coef_3)/(exp(coef_0+coef_3)+1)
  prob.mat[1,5] <- exp(coef_0+coef_4)/(exp(coef_0+coef_4)+1)
  
  # return dataframe
  return(as.data.frame(prob.mat))
  
}

# get p-values glm ~ genotype * stage
get_pvals_mat_glm <- function(model, stages) {
  
  # Initialize an empty table to store the results
  p.mat <- matrix(NA, nrow = 1, ncol = length(stages))
  colnames(p.mat) <- stages
  
  # get coefficient for genotype beta
  p_value_gt <- model[model$Coefficient=="genotypeHD", "Pr(>|z|)"]
  
  # Calculate odds ratios for each combination of genotype and stage levels
  for (j in 1:length(stages)) {
    
    if (stages[j] == "E12.5") {
      
      # For E12.5 stage, use the coefficient for HD genotype only
      p.mat[1, j] <- (p_value_gt <= 0.05)
      
    } else {
      
      coeff.name <- paste0("genotypeHD:stage", stages[j])
      
      if (coeff.name %in% model$Coefficient){
        
        # Get the coefficient and p-value for the interaction term
        p_value_gt_stage <- model[model$Coefficient==coeff.name,"Pr(>|z|)"]
        
        # Add odds ratio to the table if p-value is significant
        p.mat[1, j] <- (p_value_gt_stage<=0.05)
      }
      
    }
  }
  
  # make into dataframe
  p.df <- as.data.frame(p.mat)
  
  return(p.df)
}

# odds ratio per stage and cell_annot
get.or.glm.stage.cellannot <- function(df.test,glm.df){
  
  # Initialize a data frame with column names from the named list
  or.df <- data.frame(matrix(nrow=length(levels(df.test$cell_annot)), 
                             ncol = length(levels(df.test$stage))))
  
  colnames(or.df) <- levels(df.test$stage)
  rownames(or.df) <- levels(df.test$cell_annot)
  
  for(stg in levels(df.test$stage)){
    
    # init OR vector for stage
    stg.or <- list()
    
    for(annot in levels(df.test$cell_annot)){
      
      # check if first level
      frst.stg.lev <- stg==levels(df.test$stage)[1]
      frst.annot.lev <- annot==levels(df.test$cell_annot)[1]
      
      gt.coeff <- glm.df$Coefficient==str_interp("genotypeHD")
      gt.stg.coeff <- glm.df$Coefficient==str_interp("genotypeHD:stage${stg}")
      gt.annot.coeff <- glm.df$Coefficient==str_interp("genotypeHD:cell_annot${annot}")
      
      # get OR
      if (frst.stg.lev & frst.annot.lev) {
        or.df[[annot,stg]] <- exp(
          glm.df$Estimate[gt.coeff]*(glm.df[gt.coeff,"Pr(>|z|)"]<0.05)
        )
      } else if(frst.stg.lev & !frst.annot.lev){
        or.df[[annot,stg]] <- exp(
          glm.df$Estimate[gt.coeff]*(glm.df[gt.coeff,"Pr(>|z|)"]<0.05)+
          glm.df$Estimate[gt.annot.coeff]*(glm.df[gt.annot.coeff,"Pr(>|z|)"]<0.05)
        )
      } else if (!frst.stg.lev & frst.annot.lev){
        or.df[[annot,stg]] <- exp(
          glm.df$Estimate[gt.coeff]*(glm.df[gt.coeff,"Pr(>|z|)"]<0.05)+
          glm.df$Estimate[gt.stg.coeff]*(glm.df[gt.stg.coeff,"Pr(>|z|)"]<0.05)
        )
      } else {
        or.df[[annot,stg]] <- exp(
          glm.df$Estimate[gt.coeff]*(glm.df[gt.coeff,"Pr(>|z|)"]<0.05)+
          glm.df$Estimate[gt.stg.coeff]*(glm.df[gt.stg.coeff,"Pr(>|z|)"]<0.05) +
          glm.df$Estimate[gt.annot.coeff]*(glm.df[gt.annot.coeff,"Pr(>|z|)"]<0.05)
        )
      }
      
    }
    
  }
  
  or.df <- tibble::rownames_to_column(or.df, "cell_annot")
  
  return(or.df)
  
}

#####################################################################################################################
# Trajectory analysis
#####################################################################################################################
#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# compute nearest neighbor distance
compute_nearest_neighbor_distance <- function(coordinates) {
  
  # Compute the nearest neighbor distance for each point
  distances <- RANN::nn2(coordinates, k = 10)$nn.dists
  nearest_neighbor_distance <- rowMeans(distances)
  
  return(nearest_neighbor_distance)
}

#########################################################################################################################################################
# Spatial transcriptomics
#########################################################################################################################################################

segment.spots <- function(st.data){
  
  # get sample name
  smp.name <- unique(st.data$sample)
    
  # segment spots
  if(smp.name %in% c("BA2142","BA2143","BA2144","BA2145","BA2146")){
    st.data@meta.data$genotype <- "WT"
    st.data.hd <- subset(st.data, slice1_imagerow > 275, invert = TRUE)
    st.data@meta.data[colnames(st.data.hd),"genotype"] <- "HD"
  } else if(smp.name %in% c("BA2147","BA2148")){
    st.data@meta.data$genotype <- "WT"
    st.data.hd <- subset(st.data, slice1_imagerow > 310, invert = TRUE)
    st.data@meta.data[colnames(st.data.hd),"genotype"] <- "HD"
  } else if(smp.name %in% c("BA2149")){
    st.data@meta.data$genotype <- "WT"
    st.data.hd <- subset(st.data, slice1_imagerow > 225, invert = TRUE)
    st.data@meta.data[colnames(st.data.hd),"genotype"] <- "HD"
  } else {
    print("Sample not recognized in segmentation dictionary.") 
  }
  
  return(st.data)
  
}
