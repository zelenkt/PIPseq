#!/usr/bin/env Rscript

#can be run directly from command line by: Rscript 1-Single_Cell_RNAseq_Analysis.R

####---- User Input ----####
## Set local Github repository as working directory
setwd("/Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN")


## Specify organism
organism <- "Mouse"

# Specify the input and output folder paths
input_path <- "2-Single_Cell_RNAseq_Pipeline/Input/pipseeker_v3.3"
output_path<- "2-Single_Cell_RNAseq_Pipeline/Output"


# Specify which cancer should be processed and adjust the metadata provided

Project_Names <- c("ko_ovarian_mouse_tumor_cd8_ot1_rep1","ko_ovarian_mouse_tumor_cd8_ot1_rep2","wt_ovarian_mouse_tumor_cd8_ot1_rep1","wt_ovarian_mouse_tumor_cd8_ot1_rep2")
condition <- c("KO","KO","WT","WT")
disease <- c("ovarT","ovarT","ovarT","ovarT")
sample <- c("OTK1","OTK2","OTW1","OTW2")
# conditions for differential analysis - it will be calculated as condition1/condition2 so KO should be first
condition1 <- "KOovarT"
condition2 <- "WTovarT"
# Setup extension for labeling of integrated data
extension <- "ovarT_integrated_noRefUsed"


# Project_Names <- c("ko_melanoma_mouse_tumor_cd8_pmel_rep1","ko_melanoma_mouse_tumor_cd8_pmel_rep2","wt_melanoma_mouse_tumor_cd8_pmel_rep1","wt_melanoma_mouse_tumor_cd8_pmel_rep2")
# condition <- c("KO","KO","WT","WT")
# disease <- c("melanT","melanT","melanT","melanT")
# sample <- c("MTK1","MTK2","MTW1","MTW2")
# condition1 <- "KOmelanT"
# condition2 <- "WTmelanT"
# # Setup extension for labeling
# extension <- "melanT_integrated_noRefUsed"
# 




# list of genes to be included for plots
gene.list <- list(genes=rev(c("Tcf7","Bach2","Jun","Il7r","Sell","Ccr5","Ccr7","Cxcr3","Dpp4","Id3","Bcl2","Eomes","Tox","Tox2","Bhlhe40","Nr4a1","Nr4a2","Egr2","Pdcd1","Tigit","Lag3","Havcr2","Ctla4","Prdm1","Gzmb","Gzmc","Gzmf","Prf1","Ncr1","Ly6c2","Klrb1c","Klrd1","Klre1","Klrg1","Klrk1","Fcer1g","Fcgr3","Tyrobp")))

# genes used to infer signatures
g0 <- list(	stemness=c("Il7r","Tcf7", "Lef1", "Bach2", "Bcl6","Slamf6","Ccr7"), effector=c("Id2","Prdm1","Fasl", "Gzmb", "Prf1"), exhaustion =c("Havcr2", "Lag3", "Pdcd1", "Tigit", "Ctla4", "Nr4a1", "Nr4a2","Nr4a3","Egr2"),innate=c("Ncr1", "Ly6c2", "Klrb1c", "Klrd1", "Klre1", "Klrg1","Klrk1","Fcer1g","Tyrobp")
)
# , "Dpp4"
# "Tbx21",
gene.signatures.list <- list(g0=g0)

# basic annotation for WT clusters
annot_wt = rev(c("Il7r","Tcf7","Lef1","Bach2","Bcl6","Slamf6","Ccr5","Dpp4","Bcl2","Tbx21","Id2","Gzma","Gzmb","Prf1","Ctla4","Pdcd1","Havcr2","Nr4a1","Nr4a2","Ccl3","Egr2"))

# this is for feature and density plots
# feature_genes_list <- c("Tcf7","Bach2","Havcr2","Ctla4","Pdcd1","Jun","Tyrobp","Fcer1g","Ncr1","Hif1a") 
feature_genes_list <- gene.list[[1]]


# fcer1g related set of genes
fcer1g.signature =c("Cd247","Lat", "Ide","Xcl1","Ccr2","Cxcr3","Napsa","Adamts9","Serpinb6a","B2m","Ncr1","Klrb1c","Il2rb","Cd7","Tyrobp","Fcer1g","Otx2","Zbtb12","Mlx","Nr5a1","Otx1","Sohlh2","Fli1","Etv2","Ets2","Zbtb2","Tbx21","Tcf7")





# Setup working resolution for downstream analyses
working_res <- "integrated_snn_res.0.7"


# Specify whether the integrative analysis will be run and which condition should be used as a reference (if any)
integration <- TRUE

# Select a reference for eference-based integration. Note that we do not recommended it for standard-size datasets.
reference <- NULL
# reference <- "WT" 

# Specify what variables should be regressed out
regress.out.list <- list(c("nFeature_RNA","percent.mt"),c("nFeature_RNA","percent.mt","S.Score","G2M.Score"))
regress.out <- regress.out.list[[1]]




### Specify additional parameters
# number of variable features for integration
varFeatures <- 3000

# lower and upper filtering thresholds for nFeatures and upper threshold for mitochondrial content
MINnFeature <- 500
MAXnFeature <- 20000
mtPerc <- 20


# Specify the number of randomly sampled cells to subset for interactive shiny app (skipped with NULL parameter)
num_cells <- 1000
# num_cells <- NULL

# Define the file names for count file and/or for raw matrix data from CellRanger / PIPseeker
Count_file_name <- "Count_File.txt.gz"
Matrix_name <- "matrix.mtx.gz"
Barcodes_name <- "barcodes.tsv.gz"
Features_name <- "features.tsv.gz"

# Define the file name for clinical/meta data (OPTIONAL)
meta_name <- "Meta_File_Example_GSE116256_AML921A-D0.txt.gz"


seed <- 42


# Setup range of resultions for cluster generation
resolution <- c(0.1,0.2,0.3,0.4,0.5,0.7,0.85,1,1.2,1.35,1.5,1.75,2,2.5)

####---- Load Packages ----####

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix", "fields","SeuratDisk","tibble","SeuratData","SingleR","SingleCellExperiment", "ggplot2", "cowplot", "ggrepel", "UCell", "ggpubr","stringr","SeuratWrappers","reticulate","tidyr","scCustomize","Nebulosa","egg","reshape","reshape2")
invisible(lapply(packages, library, character.only = TRUE))



#### Define functions

# use these patches to save and load objects https://github.com/whtns/chevreul/issues/42
#' Patch of SeuratDisk::SaveH5Seurat function
#'
#' The "Assay5" class attribute "RNA" needs to be converted to a standard "Assay"
#' class for compatibility with SeuratDisk. It requires to make a temporary copy
#' so the size of the object grows bigger.
#'
#' @param object the Seurat object
#' @param filename the file path where to save the Seurat object
#' @param verbose SaveH5Seurat verbosity
#' @param overwrite whether to overwrite an existing file
#' 
#' @return NULL
SaveH5SeuratObject <- function(
    object,
    filename,
    verbose = TRUE,
    overwrite = TRUE
) {
  
  # add copy of "RNA" 
  object[["RNA-tmp"]] <- CreateAssayObject(counts = object[["RNA"]]$counts)
  # remove original
  object[["RNA"]] <- NULL
  # export
  SaveH5Seurat(object, filename = filename, overwrite, verbose)
  
  return(NULL)
}

#' Patch of SeuratDisk::LoadH5Seurat function
#'
#' The "Assay" class attribute "RNA" needs to be converted to the new "Assay5"
#' class. It requires to make a temporary copy so the size of the object grows
#' bigger.
#'
#' @param filename the file path where to save the Seurat object
#' @param verbose LoadH5Seurat verbosity
#' 
#' @return NULL
LoadH5SeuratObject <- function(
    filename,
    verbose = TRUE
) {
  
  # load h5 data
  object <- LoadH5Seurat(filename)
  # create "Assay5" class from old "Assay" class
  keys <- Key(object)
  slotID <- names(keys)[startsWith(keys, "rna")]
  object[["RNA"]] <- CreateAssay5Object(counts = object[[slotID]]$counts)
  # delete "Assay" class
  object[[slotID]] <- NULL
  # reorder assays list
  object@assays <- object@assays[sort(names(object@assays))]
  
  return(object)
}






# Define a function for quality control (QC) using Seurat
qc.seurat <- function(seurat, species, MINnFeature, MAXnFeature, mtPerc) {
  mt.pattern <- case_when(
    species == "Human" ~ "^MT-",
    species == "Mouse" ~ "^mt-",
    TRUE ~ "^MT-"
  )
  ribo.pattern <- case_when(
    species == "Human" ~ "^RP[LS]",
    species == "Mouse" ~ "^Rp[ls]",
    TRUE ~ "^RP[LS]"
  )
  
  # Calculate percentage of mitochondrial and ribosomal genes
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern, assay = "RNA")
  seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern, assay = "RNA")
  
  # Filter cells based on QC criteria
  seurat[, seurat[["percent.mt"]] <= mtPerc & seurat[["nFeature_RNA"]] >= MINnFeature & seurat[["nFeature_RNA"]] <= MAXnFeature]
}

# Define a function to annotate and save data
annotate.save.data <- function(seurat_obj, organism, Project_Name, meta, output_folder) {
  
  if (organism=="Mouse") { 
    
    # Load mouse reference databases from celldex
    Immgen.ref <- celldex::ImmGenData()
    mouseRNAseq.ref <- celldex::MouseRNAseqData()
    # Convert Seurat object to SingleCellExperiment
    sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
    
    # Auto-annotate cell types using reference databases from celldex
    ImmGen.main <- SingleR(test = sce, assay.type.test = 1, ref = Immgen.ref, labels = Immgen.ref$label.main)
    ImmGen.fine <- SingleR(test = sce, assay.type.test = 1, ref = Immgen.ref, labels = Immgen.ref$label.fine)
    mouseRNAseq.main <- SingleR(test = sce, assay.type.test = 1, ref = mouseRNAseq.ref, labels = mouseRNAseq.ref$label.main)
    mouseRNAseq.fine <- SingleR(test = sce, assay.type.test = 1, ref = mouseRNAseq.ref, labels = mouseRNAseq.ref$label.fine)
    
    # Add the celldex annotations to the metadata of the Seurat object
    seurat_obj@meta.data$ImmGen.main <- ImmGen.main$pruned.labels
    seurat_obj@meta.data$ImmGen.fine <- ImmGen.fine$pruned.labels
    seurat_obj@meta.data$mouseRNAseq.main <- mouseRNAseq.main$pruned.labels
    seurat_obj@meta.data$mouseRNAseq.fine <- mouseRNAseq.fine$pruned.labels
    
  } else if (organism=="Human") {
    
    # Load reference databases from celldex
    hpca.ref <- celldex::HumanPrimaryCellAtlasData()
    dice.ref <- celldex::DatabaseImmuneCellExpressionData()
    blueprint.ref <- celldex::BlueprintEncodeData()
    monaco.ref <- celldex::MonacoImmuneData()
    northern.ref <- celldex::NovershternHematopoieticData()
    
    # Convert Seurat object to SingleCellExperiment
    sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
    
    # Auto-annotate cell types using reference databases from celldex
    hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
    hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
    dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
    dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
    blue.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
    blue.fine <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.fine)
    monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
    monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
    northern.main <- SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.main)
    northern.fine <- SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.fine)
    
    # Add the celldex annotations to the metadata of the Seurat object
    seurat_obj@meta.data$hpca.main <- hpca.main$pruned.labels
    seurat_obj@meta.data$hpca.fine <- hpca.fine$pruned.labels
    seurat_obj@meta.data$dice.main <- dice.main$pruned.labels
    seurat_obj@meta.data$dice.fine <- dice.fine$pruned.labels
    seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
    seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
    seurat_obj@meta.data$northern.main <- northern.main$pruned.labels
    seurat_obj@meta.data$northern.fine <- northern.fine$pruned.labels
    seurat_obj@meta.data$blue.main <- blue.main$pruned.labels
    seurat_obj@meta.data$blue.fine <- blue.fine$pruned.labels
    
    
  } else {
    print("Please check what organism you used")
  }
  
  # Add UMAP coordinates to the metadata
  UMAP <- as.data.frame(Embeddings(object = seurat_obj[["umap"]]))
  seurat_obj@meta.data$UMAP_1 <- UMAP$UMAP_1
  seurat_obj@meta.data$UMAP_2 <- UMAP$UMAP_2
  
  # Add tSNE coordinates to the metadata
  tsne <- as.data.frame(Embeddings(object = seurat_obj[["tsne"]]))
  seurat_obj@meta.data$tSNE_1 <- tsne$tSNE_1
  seurat_obj@meta.data$tSNE_2 <- tsne$tSNE_2
  
  
  meta_out <- seurat_obj@meta.data
  meta_out$Cell <- rownames(meta_out)
  meta_out <- meta_out %>% relocate(Cell)
  if (!is.null(meta)) {
    meta_out <- merge(meta,meta_out, by = "Cell", all.y = T)
    meta_obj <- meta_out
    rownames(meta_obj) <- meta_obj$Cell
    meta_obj <- meta_obj[,-1]
    seurat_obj@meta.data <- meta_obj
  }
  
  # Write metadata to a file
  write_delim(meta_out, file = file.path(output_folder, paste0(Project_Name,"_metafile_with_annotation.txt")),
              delim = "\t")
  
  DefaultAssay(seurat_obj) <-  'RNA' # temporarily making 'RNA' active assay
  seurat_obj <- FindVariableFeatures(seurat_obj)
  DefaultAssay(seurat_obj) <-  'SCT' # returning 'SCT' as the default assay
  
  
  # Write seurat H5 file 
  SaveSeuratRds(object = seurat_obj,file = file.path(output_folder, paste0(Project_Name,".Rds")))
  
  SaveH5SeuratObject(seurat_obj, filename = file.path(output_folder, paste0(Project_Name,"_h5friendly")),
               overwrite = TRUE, verbose = TRUE)
  
  
  
  
  # write out count data
  for (j in names(seurat_obj@assays)) {
      if (j =="RNA") {
      print("Exporting from RNA assay")
      # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
      # scaled_counts <- as.data.frame(seurat_obj@assays[[j]]@scale.data)
      raw_counts <- as.data.frame(seurat_obj@assays[[j]]@layers$counts)
      normalized_data <- as.data.frame(seurat_obj@assays[[j]]@layers$data)
      } else if (j =="SCT") {
      print("Exporting from SCT assay")
      raw_counts <- as.data.frame(seurat_obj@assays[[j]]@counts)
      normalized_data <- as.data.frame(seurat_obj@assays[[j]]@data)
      } else {print("Make sure a correct assay is generated")}
  
    
    
    # Add row names as a column for Alyssa's Umap app
    # scaled_counts <- tibble::rownames_to_column(scaled_counts, var = "Genes")
    raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
    normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
    
    # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
    # write.table(scaled_counts, file = paste0(output_folder, "/", Project_Name,"_",j,"_ScaledCounts.txt"),
    #             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(raw_counts, file = paste0(output_folder, "/", Project_Name,"_",j,"_raw_counts.txt"),
                sep = "\t", row.names = FALSE, col.names = TRUE)
    write.table(normalized_data, file = paste0(output_folder, "/", Project_Name,"_",j,"_normalized_counts.txt"),
                sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

integrate.save.data <- function(output_path, output_path.list, extension, Project_Names, varFeatures, integration, reference, condition, resolution, seed) {

  
  # Retrive save seurat objects
  for (i in 1:length(Project_Names)){
    Project_Name <- Project_Names[i]
    if (!file.exists(paste0(output_path.list[i]))) {
      print(paste0("Seurat object for ",Project_Name," does not exist, repeat the first part of the analysis..."))
    } else {
      print(paste0(Project_Name," seurat object exists and it will be loaded..."))
      seurat_obj <- LoadSeuratRds(paste0(output_path.list[i]))
      assign(paste0("seurat_obj",i),seurat_obj)
    }
  }
  
  # create a list of seurat objects
  seurat.object.list <- list()
  for (i in 1:length(Project_Names)){
    seurat.object.list <- append(seurat.object.list,eval(parse(text=paste0("seurat_obj",i))))
  }
  
  
  # integrate data
  if (integration == TRUE) {
    print(paste0("Integrative analysis will be performed..."))
    
    # perform integrated analysis
    features <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = varFeatures)
    seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = features)
    
    if (!(is.null(reference))) {
      print(paste0(reference," will be used as a reference for integration..."))
      
      anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = features, reference = grep(reference,condition), reduction = "rpca")
      
    } else {
      
      print(paste0("No reference has been set, all data will be used for integration..."))
      
      anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = features)
    }
    
    integrated.seurat.object <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)
    
    integrated.seurat.object <- RunPCA(integrated.seurat.object, verbose = FALSE, seed.use = seed) %>%
      RunUMAP(verbose = FALSE, dims = 1:50,seed.use = seed) %>% 
      RunTSNE(reduction = "pca", dims = 1:50, seed.use = seed) %>% 
      FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>% 
      FindClusters(resolution = resolution, verbose = FALSE) 
    
    
    DefaultAssay(integrated.seurat.object) <- "RNA"
    integrated.seurat.object <- NormalizeData(integrated.seurat.object, verbose = FALSE)
    integrated.seurat.object <- FindVariableFeatures(integrated.seurat.object, selection.method = "vst", nfeatures = 5000)
    
    integrated.seurat.object@meta.data$diseaseCondition <- paste0(integrated.seurat.object@meta.data$disease,integrated.seurat.object@meta.data$condition)
    
    DefaultAssay(integrated.seurat.object) <- "integrated"
    

    # Add UMAP coordinates to the metadata
    UMAP <- as.data.frame(Embeddings(object = integrated.seurat.object[["umap"]]))
    integrated.seurat.object@meta.data$UMAP_1 <- UMAP$UMAP_1
    integrated.seurat.object@meta.data$UMAP_2 <- UMAP$UMAP_2
      
    # Add tSNE coordinates to the metadata
    tsne <- as.data.frame(Embeddings(object = integrated.seurat.object[["tsne"]]))
    integrated.seurat.object@meta.data$tSNE_1 <- tsne$tSNE_1
    integrated.seurat.object@meta.data$tSNE_2 <- tsne$tSNE_2
      

    SaveH5SeuratObject(integrated.seurat.object, filename = file.path(integrated.path,paste0("01_",extension,"_seurat_obj")), overwrite = TRUE, verbose = TRUE)
    SaveSeuratRds(integrated.seurat.object, file = file.path(integrated.path,paste0("01_",extension,"_seurat_obj.Rds")))
    
    # write out count data
    joinedLayers <- JoinLayers(integrated.seurat.object, assay = "RNA")
    raw_counts <- as.data.frame(joinedLayers@assays$RNA@layers$counts)
    normalized_data <- as.data.frame(joinedLayers@assays$RNA@layers$data)
    metadata <- as.data.frame(joinedLayers@meta.data)
    
    
    # raw_counts <- as.data.frame(integrated.seurat.object@assays$RNA@counts)
    # normalized_data <- as.data.frame(integrated.seurat.object@assays$RNA@data)
    # metadata <- as.data.frame(integrated.seurat.object@meta.data)
    
    # Add row names as a column for the Umap app
    raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
    normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
    metadata <- tibble::rownames_to_column(metadata, var = "Cell")
    
    # Write raw_counts, normalized_data, and metadata to separate files
    write.table(raw_counts, file = paste0(integrated.path,"/",extension,"_raw_counts.txt"),sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(normalized_data, file = paste0(integrated.path,"/",extension,"_normalized_counts.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(metadata, file = paste0(integrated.path,"/",extension,"_metafile.txt"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    return(integrated.seurat.object)
    
  } else {
    print("Integration analysis will be skipped....")
  }
}


####---- Run Script ----####
set.seed(seed)


# Process the raw counts for all the samples and perform QC filtering
for (i in 1:length(Project_Names)){
  Project_Name <- Project_Names[i]
  
  output_folder <-  paste0(output_path,"/",Project_Name,"/Single_Cell_RNAseq_Output")
  if (!file.exists(output_folder)) {
    dir.create(output_folder,recursive = T)
  }
  
  
  meta_file <- paste0(input_path,"/",Project_Name,"/",meta_name)
  Count_file <-  paste0(input_path,"/",Project_Name,"/",Count_file_name)
  Matrix_File <- paste0(input_path,"/",Project_Name,"/",Matrix_name)
  Barcode_File <- paste0(input_path,"/",Project_Name,"/",Barcodes_name)
  Feature_File <- paste0(input_path,"/",Project_Name,"/",Features_name)
  
  if (file.exists(meta_file)) {
    meta <- as.data.frame(read_delim(meta_file, delim = '\t', col_names = T))
    meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
    colnames(meta)[1] <- "Cell"
  } else {meta <- NULL}
  
  
  
  # Read the count file
  if (file.exists(Count_file)) {
    counts <- read.csv2(Count_file, sep = "", header = TRUE, row.names = 1)
    colnames(counts) <- gsub("[[:punct:]]",".",colnames(counts))
  } else if (file.exists(Matrix_File)) {
    counts <- ReadMtx(mtx = Matrix_File,  #mtx file
                      features = Feature_File, #feature file
                      cells = Barcode_File) #barcode file
  } else {
    print("Please define a correct input file")
  }
  
  
  # Create a Seurat object and perform QC
  seurat_obj <- CreateSeuratObject(counts = counts, min.cells=3)
  
  seurat_obj <- qc.seurat(seurat_obj, organism, 100, 30000, 100)
  
  # save QC graph before filtering 
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA","percent.rp", "percent.mt"), ncol = 4)
  ggsave2(paste0(output_path,"/",Project_Name,"/QC_prefilt_",Project_Name,".pdf"), width = 7, height = 5)
  
  seurat_obj <- qc.seurat(seurat_obj, organism, MINnFeature, MAXnFeature, mtPerc)
  
  # save QC graph post filtering
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA","percent.rp", "percent.mt"), ncol = 4)
  ggsave2(paste0(output_path,"/",Project_Name,"/QC_postfilt_",Project_Name,".pdf"), width = 7, height = 5)
  
  # Analyze data by Seurat if h5seurat file does not exist...
  if (!file.exists(paste0(output_folder,"/",Project_Name,"_h5friendly.h5seurat"))) {
    print(paste0("Processing scRNAseq raw counts for ",Project_Name," by Seurat..."))
    
    
    # add sample related metadata
    seurat_obj@meta.data$condition <- paste0(condition[i])
    seurat_obj@meta.data$disease <- paste0(disease[i])
    seurat_obj@meta.data$sample <- paste0(sample[i])
    seurat_obj@meta.data$conditionDisease <- paste0(condition[i],disease[i])
    seurat_obj@meta.data$diseaseCondition <- paste0(disease[i],condition[i])
    # append sample name prefix to barcodes, automatically adds _ separator
    seurat_obj <- RenameCells(seurat_obj, add.cell.id= paste0(sample[i]))
    
    DefaultAssay(seurat_obj) <- "RNA"
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- CellCycleScoring(seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, verbose = FALSE, seed.use = seed)
    
    
    seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE, seed.use = seed, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress=regress.out)

    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, seed.use = seed) %>% 
      RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE, seed.use = seed) %>% 
      RunTSNE(reduction = "pca", dims = 1:50, seed.use = seed) %>% 
      FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>% 
      FindClusters(resolution = resolution, verbose = FALSE) 
    
    # Find doublets using DoubletFinder
    suppressMessages(require(DoubletFinder))
    nExp <- round(ncol(seurat_obj) * 0.05)  # expect 5% doublets
    seurat_obj <- doubletFinder(seurat_obj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = T)
    
    # rename doublet column
    df <- names(seurat_obj@meta.data)
    names(seurat_obj@meta.data)[names(seurat_obj@meta.data) == paste0(df[grep("^DF.classifications_",names(seurat_obj@meta.data))])] <- "doublets"
    
    # annotate and save data for individual datasets
    seurat_obj <- annotate.save.data(seurat_obj, organism, Project_Name, meta, output_folder)
    
    # # save the seurat object to the environment if needed
    #  assign(paste0("seurat_obj",i),seurat_obj)

  } else {
    print(paste0(Project_Name," raw counts were already processed, skipping..."))
  }
}


# integrate and save the data
output_path.list <- c()
for (i in 1:length(Project_Names)){
  path <- paste0(output_path,"/",Project_Names[i],"/Single_Cell_RNAseq_Output/",Project_Names[i],".Rds")
  output_path.list <- append(output_path.list,path)
}

integrated.path <- paste0(output_path,"/",extension)
if (!file.exists(integrated.path)) {
  dir.create(integrated.path,recursive = T)
}

integrated.seurat.object <- integrate.save.data(output_path, output_path.list, extension, Project_Names, varFeatures, integration, reference, condition, resolution, seed)






# ######### generate more detailed graphs for expression and signatures
# integrated.seurat.object <- LoadH5Seurat(paste0(integrated.path,"/01_",extension,"_seurat_obj.h5seurat"))
# integrated.seurat.object <- LoadSeuratRds(paste0(integrated.path,"/01_",extension,"_seurat_obj.Rds"))

# subset only small number of cells to be used for shiny app
if (is.null(num_cells)) {
      print(paste0("The number of cells to be subset was not defined, only the full dataset/datasets will be saved..."))
            
  } else {
      
    print(paste0(num_cells," cells will be randomly selected and saved to be used for interactive shiny application..."))

    # Randomly sample cells
    random_cells <- sample(1:ncol(integrated.seurat.object), size = num_cells, replace = FALSE)
    subset_seurat_obj <- integrated.seurat.object[,random_cells]

    # # write out count data
    # for (j in 1:length(names(subset_seurat_obj@assays))) {
    #   # write out count data
    #   joinedLayers <- JoinLayers(subset_seurat_obj, assay =  names(subset_seurat_obj@assays)[j])
    #   
      joinedLayers <- JoinLayers(subset_seurat_obj, assay = "RNA")
      
      raw_counts <- as.data.frame(joinedLayers@assays$RNA@layers$counts)
      normalized_data <- as.data.frame(joinedLayers@assays$RNA@layers$data)
      metadata <- as.data.frame(joinedLayers@meta.data)

      # Add row names as a column for the Umap app
      raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
      normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
      metadata <- tibble::rownames_to_column(metadata, var = "Cell")
      # # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
      # raw_counts <- as.data.frame(subset_seurat_obj@assays[[j]]@counts)
      # normalized_data <- as.data.frame(subset_seurat_obj@assays[[j]]@data)
      # metadata <- as.data.frame(subset_seurat_obj@meta.data)
      # 
      # # Add row names as a column for Alyssa Umap app
      # raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
      # normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
      # metadata <- tibble::rownames_to_column(metadata, var = "Cell")
      
      # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
      write.table(raw_counts, file = paste0(integrated.path, "/", extension,"_",num_cells,"_RNA_raw_counts.txt"),
                  sep = "\t", row.names = FALSE, col.names = TRUE)
      write.table(normalized_data, file = paste0(integrated.path, "/", extension,"_",num_cells,"_RNA_normalized_counts.txt"),
                  sep = "\t", row.names = FALSE, col.names = TRUE)
      write.table(metadata, file = paste0(integrated.path,"/",extension,"_",num_cells,"_RNA_metafile.txt"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
}




figures.path <- paste0(integrated.path,"/figures_raw")
if (!file.exists(figures.path)) {
  dir.create(figures.path,recursive = T)
}

  # joinLayers is needed when using Seurat v5
  integrated.seurat.object <- JoinLayers(integrated.seurat.object, assay = "RNA")
  working_assay <- "RNA"
  DefaultAssay(integrated.seurat.object) <- "RNA"
  Idents(integrated.seurat.object) <- working_res
  markers <- FindAllMarkers(integrated.seurat.object, assay = working_assay, random.seed = seed) # only.pos = TRUE
  
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene))  %>% #filters out all ribosomal and mitochondrial genes
    slice_head(n = 10) %>%
    ungroup() -> top10
  top10 = top10[!duplicated(top10$gene),]
  
  
  write.table(markers, file = paste0(figures.path,"/",extension,"_allMarkers_combined_cond_RNA_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  

  DotPlot(integrated.seurat.object, features = top10$gene,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip() +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_dotplot_allMarkers_condDis_noRiboMt_assay",working_assay,"_",working_res,".svg"), width = 10, height = 20)
  
  for (g in 1:length(gene.list)){
  genes <- gene.list[g]

  DotPlot(integrated.seurat.object, features = genes,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip()  +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_dotplot_",names(gene.list[g]),"_condDis_assay",working_assay,"_",working_res,".svg"), width = 10, height = 9)
  }

  
  # fix the umap coordinates for feature and signature plots https://github.com/satijalab/seurat/issues/2507 https://satijalab.org/seurat/articles/essential_commands.html
  DefaultAssay(integrated.seurat.object) <- "integrated"
  new.embeddings = Embeddings(integrated.seurat.object, reduction = "umap")
  new_reduction <- CreateDimReducObject(embeddings = new.embeddings, key = "umap")
  integrated.seurat.object[["umap"]] <- new_reduction

  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = working_res) + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_clusters_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_cond_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = "Phase") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_phase_",working_res,".svg"), width = 7, height = 5)

  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = "sample") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_sample_",working_res,".svg"), width = 7, height = 5)

  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = working_res, split.by = "diseaseCondition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-diseaseCondition_",working_res,".svg"), width = 12, height = 5)
  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, split.by = "disease",group.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-disease_groupCond_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T, group.by = "conditionDisease") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_conditionDisease_",working_res,".svg"), width = 7, height = 5)

  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T,group.by = working_res, split.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-condition_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(integrated.seurat.object, reduction= "umap", label = T, repel = T,group.by = working_res, split.by = "disease") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-disease_",working_res,".svg"), width = 7, height = 5)

  

  DefaultAssay(integrated.seurat.object) <- "RNA"
 
  # differential analysis between two conditions
  integrated.seurat.object$conditionDisease.clust <- paste(integrated.seurat.object$conditionDisease, integrated.seurat.object[[working_res]][,1],sep = "_")
  
  markers.path <- paste0(figures.path,"/markers")
  if (!file.exists(markers.path)) {
    dir.create(markers.path,recursive = T)
  }
  
  for (i in 0:(length(table(integrated.seurat.object@meta.data[[working_res]]))-1))  {
    tryCatch({
      DefaultAssay(integrated.seurat.object) <- working_assay
      Idents(integrated.seurat.object) <- "conditionDisease.clust"
      # finds markers within each cluster differential between the two conditions / replicates
      
      markers <- FindMarkers(integrated.seurat.object, assay = working_assay, ident.1 = paste0(condition1,"_",i), ident.2 = paste0(condition2,"_",i), verbose = T, random.seed = seed)
      
      markers$gene <- rownames(markers)
      
      write.table(markers, file = paste0(markers.path,"/cl_",i,"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      for (g in 1:length(gene.list)){
      genes <- as.data.frame(gene.list[g])
      markers$label <- markers$gene %in% genes[,1]
      
      markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
      
      ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
        data = subset(markers, label == TRUE),
        aes(label = gene),
        size = 7,
        box.padding = unit(3, "lines"),
        point.padding = unit(0.25, "lines"),
        max.overlaps = Inf) +  scale_color_manual(values = c("#7ac5cdff","black","#f08080ff" )) + ggtitle(paste0(disease," CD8 T cells, cluster ",i," ",condition1,"/",condition2))
      
      ggsave(paste0(markers.path,"/volcano_cl",i,"_",names(gene.list[g]),"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  


############## Calculate and plot signature scores
signatures.path <- paste0(figures.path,"/signatures")
if (!file.exists(signatures.path)) {
  dir.create(signatures.path,recursive = T)
}

DefaultAssay(integrated.seurat.object) <- "RNA"
Idents(integrated.seurat.object) <- working_res

for (i in 1:(length(gene.signatures.list))) {
  signatures <- AddModuleScore_UCell(integrated.seurat.object, features = gene.signatures.list[[i]])
  signature.names <- paste0(names(gene.signatures.list[[i]]), "_UCell")
  
  # generate violin plots
  VlnPlot(signatures, features = signature.names, group.by = working_res, split.by = "diseaseCondition") # order: melanKO, melanWT, ovarKO, ovarWT
  ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 20, height = 5)
  
  # generate feature plots
  FeaturePlot(signatures, features = signature.names,split.by = "diseaseCondition", ncol = 3, order = T,pt.size=0.9)
  ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_featplot_split-disCond_",working_res,".svg"), width = 17, height = 25)
  
}








#######################################
# after investigating the data, perform cluster filtering and merging of clusters

reclustered.path <- paste0(integrated.path,"/reclustered")
if (!file.exists(reclustered.path)) {
  dir.create(reclustered.path,recursive = T)
}

# SaveSeuratRds(integrated.seurat.object, file = file.path(integrated.path, paste0("03_",extension,"_seurat_obj.Rds")))

# seurat_obj_melan <- LoadSeuratRds(file.path(integrated.path, paste0("03_",extension,"_seurat_obj.Rds")))
# 
# seurat_obj_ovar <- LoadSeuratRds(paste0("2-Single_Cell_RNAseq_Pipeline/Output/ovarT_integrated_noRefUsed/03_ovarT_integrated_noRefUsed_seurat_obj.Rds"))

     
reg.status.list <- c("notRegCC","regCC")
reg.status <- reg.status.list[1]

seurat_obj <- integrated.seurat.object

seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj$tSNE_1 <- seurat_obj@reductions$tsne@cell.embeddings[,1]
seurat_obj$tSNE_2 <- seurat_obj@reductions$tsne@cell.embeddings[,2]

working_assay <- "RNA"



# plot raw DimPlots
DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.3",split.by = "condition", reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_splitWTKO_nofilt_res03_",reg.status,"_",extension,".svg"), width = 9, height = 5.5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.4",split.by = "condition", reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_splitWTKO_nofilt_res04_",reg.status,"_",extension,".svg"), width = 9, height = 5.5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.5",split.by = "condition", reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_splitWTKO_nofilt_res05_",reg.status,"_",extension,".svg"), width = 9, height = 5.5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.7",split.by = "condition",reduction = "tsne",  label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_splitWTKO_nofilt_res07_",reg.status,"_",extension,".svg"), width = 9, height = 5.5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.1.5",split.by = "condition",reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_splitWTKO_nofilt_res15_",reg.status,"_",extension,".svg"), width = 9, height = 5.5)

DimPlot(object = seurat_obj, group.by = "condition",reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/01_dimplot_WTKOboth_nofilt_",reg.status,"_",extension,".svg"), width = 6, height = 6)

seurat_obj$FinalClusterNameCondition <- paste(seurat_obj$integrated_snn_res.1.5,seurat_obj$condition,sep = "_")
Idents(seurat_obj) <- "FinalClusterNameCondition"


# see the clusters with low number of cells. Remove clusters with less than 10 cells.
table(seurat_obj$FinalClusterNameCondition)
df <- as.data.frame(table(seurat_obj$FinalClusterNameCondition))
df <- subset(df,Freq<=10)

if (reg.status=="regCC") { 
  seurat_obj = subset(x = seurat_obj, idents = as.vector(df$Var1[-1]), invert = T) 
  
} else if (reg.status=="notRegCC") {
  seurat_obj = subset(x = seurat_obj, idents = as.vector(df$Var1), invert = T) # for noCCregress us as.vector(df$Var1)
  
} else {
  print("Check if correct regression status is selected")
}


# remove the contaminated clusters in melanoma dataset
if (disease[1]=="ovarT") {
  print("No additional filtering needed for ovarian tumor samples")
  
} else if (disease[1]=="melanT") {
  # this indicates that cluster 10+13 at 1.5 resolution are likely contaminated with melanoma cells so these two clusters will be removed
  DotPlot(seurat_obj, features = c("Mitf","S100b","Mcam","Mlana","Trac","Cd3e","Cd8a","Cd8b1"),group.by = "FinalClusterNameCondition",assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  ggsave2(paste0(reclustered.path,"/01_bubble_bothWTKO_melanGenes_filt_res15_",reg.status,"_",extension,".svg"), width = 10, height = 4)
  
  
  print("Clusters contaminated with melanoma-specific genes will be removed")
  Idents(seurat_obj) <- "integrated_snn_res.1.5"
  seurat_obj = subset(x = seurat_obj, idents = "10", invert = T)
  seurat_obj = subset(x = seurat_obj, idents = "13", invert = T)
  
  DotPlot(seurat_obj, features = c("Mitf","S100b","Mcam","Mlana","Trac","Cd3e","Cd8a","Cd8b1"),group.by = "FinalClusterNameCondition",assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  ggsave2(paste0(reclustered.path,"/01b_bubble_bothWTKO_melanGenes_filt_res15_",reg.status,"_",extension,".svg"), width = 10, height = 4)
  
  # plot filtered DimPlots, note that tSNE reflects much better the clustering than UMAP so we stick to that
  DimPlot(object = seurat_obj, group.by = "integrated_snn_res.1.5",reduction = "tsne",split.by = "condition", label = T)
  ggsave2(paste0(reclustered.path,"/01b_dimplot_splitWTKO_filt_res15_",reg.status,"_",extension,".svg"), width = 8, height = 5)
  
  DimPlot(object = seurat_obj, group.by = "integrated_snn_res.1.5", reduction = "tsne",label = T)
  ggsave2(paste0(reclustered.path,"/01b_dimplot_noSplit_filt_res15_",reg.status,"_",extension,".svg"), width = 5, height = 5)

} else {
  print("Check if correct disease is selected")
}


# keep a backup for testing
seurat_obj_back <- seurat_obj # main backup
# seurat_obj <- seurat_obj_back


# rerun clustering, umaps and tsne plot - I did that for testing but I prefer the shape of the original tSNE plot for visualization purposes so I add back the original tSNE coordinates
DefaultAssay(seurat_obj) <-  'integrated'
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE,npcs = 50, seed.use = seed) %>%
  RunUMAP(verbose = FALSE, dims = 1:50,seed.use = seed, min.dist = 0.7) %>% 
  RunTSNE(reduction = "pca", dims = 1:50, seed.use = seed) %>% FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%  FindClusters(resolution = resolution, verbose = FALSE) 
# # check the primary sources of heterogeneity and number of dimensions to be used..even at 50 we still see some important genes somewhat separating different groups of cells
# DimHeatmap(seurat_obj, dims = 50:40, cells = 1000, balanced = TRUE)


########################## testing
# check the original tSNE coordinates
DimPlot(object = seurat_obj_back, group.by = "integrated_snn_res.0.5",split.by = "condition",reduction = "tsne", label = T)

# check the newly created tSNE coordinates
seurat_obj2 <- seurat_obj
seurat_obj2$tSNE_1 <- seurat_obj2@reductions$tsne@cell.embeddings[,1]
seurat_obj2$tSNE_2 <- seurat_obj2@reductions$tsne@cell.embeddings[,2]
DimPlot(object = seurat_obj2, group.by = "integrated_snn_res.0.5",split.by = "condition",reduction = "tsne", label = T)

rm(seurat_obj2)
rm(seurat_obj_back)
########################## end testing

# keep the original tSNE coordinates for ovarian model but use the updated coordinates for melanoma

if (disease[1]=="ovarT") {
  # # for ovarian
  # this replaces the new tSNE coorindates by the original ones
  seurat_obj@reductions$tsne@cell.embeddings[,1] <- seurat_obj$tSNE_1
  seurat_obj@reductions$tsne@cell.embeddings[,2]  <- seurat_obj$tSNE_2
  
} else if (disease[1]=="melanT") {
  seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
  seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
  seurat_obj$tSNE_1 <- seurat_obj@reductions$tsne@cell.embeddings[,1]
  seurat_obj$tSNE_2 <- seurat_obj@reductions$tsne@cell.embeddings[,2]

} else {
  print("Check if correct disease is selected")
}

# re-scale data
DefaultAssay(seurat_obj) <-  'RNA'
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress=regress.out)

Idents(seurat_obj) <- "condition"
seurat_obj_wt = subset(x = seurat_obj, idents = "KO", invert = T)

seurat_obj_wt$FinalClusterNameCondition <- paste(seurat_obj_wt$integrated_snn_res.0.4,seurat_obj_wt$condition,sep = "_")
Idents(seurat_obj_wt) <- "FinalClusterNameCondition"

seurat_obj$FinalClusterNameCondition <- paste(seurat_obj$integrated_snn_res.0.4,seurat_obj$condition,sep = "_")
Idents(seurat_obj) <- "FinalClusterNameCondition"

DotPlot(seurat_obj_wt, features = gene.list,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_WTonly_genelist_filt_res04_",reg.status,"_",extension,".svg"), width = 5, height = 6)

DotPlot(seurat_obj_wt, features = annot_wt,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_WTonly_WTannot_filt_res04_",reg.status,"_",extension,".svg"), width = 5, height = 4.5)

DotPlot(seurat_obj, features = gene.list,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_bothWTKO_genelist_filt_res04_",reg.status,"_",extension,".svg"), width = 9.3, height = 6.5)


seurat_obj$FinalClusterNameCondition <- paste(seurat_obj$integrated_snn_res.1.5,seurat_obj$condition,sep = "_")
Idents(seurat_obj) <- "FinalClusterNameCondition"
seurat_obj_wt$FinalClusterNameCondition <- paste(seurat_obj_wt$integrated_snn_res.1.5,seurat_obj_wt$condition,sep = "_")
Idents(seurat_obj_wt) <- "FinalClusterNameCondition"


DotPlot(seurat_obj_wt, features = gene.list,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_WTonly_genelist_filt_res15_",reg.status,"_",extension,".svg"), width =10, height = 6)

DotPlot(seurat_obj_wt, features = annot_wt,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_WTonly_WTannot_filt_res15_",reg.status,"_",extension,".svg"), width = 10, height = 4.5)

DotPlot(seurat_obj, features = gene.list,group.by = "FinalClusterNameCondition", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
ggsave2(paste0(reclustered.path,"/02_bubble_bothWTKO_genelist_filt_res15_",reg.status,"_",extension,".svg"), width = 10, height = 6)





# plot filtered DimPlots, note that tSNE reflects much better the clustering than UMAP so we stick to that
DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.4",split.by = "condition",reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_splitWTKO_filt_res04_",reg.status,"_",extension,".svg"), width = 8, height = 5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.0.5",reduction = "tsne",split.by = "condition", label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_splitWTKO_filt_res05_",reg.status,"_",extension,".svg"), width = 8, height = 5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.1.5",reduction = "tsne",split.by = "condition", label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_splitWTKO_filt_res15_",reg.status,"_",extension,".svg"), width = 8, height = 5)

DimPlot(object = seurat_obj, group.by = "integrated_snn_res.1.5", reduction = "tsne",label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_noSplit_filt_res15_",reg.status,"_",extension,".svg"), width = 5, height = 5)

DimPlot(object = seurat_obj, group.by = "Phase",reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_Phase_filt_",reg.status,"_",extension,".svg"), width = 5, height = 5)

DimPlot(object = seurat_obj, group.by = "sample",reduction = "tsne", label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_Sample_filt_",reg.status,"_",extension,".svg"), width = 5, height = 5)

DimPlot(object = seurat_obj, group.by = "sample",reduction = "tsne",split.by = "condition",  label = T)
ggsave2(paste0(reclustered.path,"/02_dimplot_SampleSplit_filt_",reg.status,"_",extension,".svg"), width = 8, height = 5)

DimPlot(object = seurat_obj, group.by = "condition",reduction = "tsne",label = F)
ggsave2(paste0(reclustered.path,"/02_dimplot_WTKOboth_filt_",reg.status,"_",extension,".svg"), width = 5, height = 5)


# separately prepared markers for Phase2 - grouped cell cycle phases
seurat_obj$Phase2 <- paste(ifelse(seurat_obj$Phase == "G1","G1","S/G2/M"))

seurat_obj$Phase2_cond <- paste(seurat_obj$Phase2, seurat_obj$condition, sep = "_")
DimPlot_scCustom(seurat_object = seurat_obj, group.by = "Phase2",reduction = "tsne", label = T)+ ggtitle("Cell cycle phase")
ggsave2(paste0(reclustered.path,"/04_dimplot_Phase2_noSplit_filt_resMix_",extension,".svg"), width = 5.1, height = 5)

Idents(seurat_obj) <- "Phase2"
markers <- FindAllMarkers(seurat_obj, assay = working_assay, random.seed = seed, only.pos = TRUE) # only.pos = TRUE
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 10) %>%
  ungroup() -> top

top = top[!duplicated(top$gene),]
seurat_obj$Phase2_cond <- factor(seurat_obj$Phase2_cond , levels=c("G1_WT","G1_KO","S/G2/M_WT","S/G2/M_KO"))
DotPlot(seurat_obj, features = top$gene,group.by = "Phase2_cond", assay=working_assay, scale.min = 0,scale.max =70 ,scale.by = 'size') + coord_flip() +RotatedAxis()+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
ggsave(paste0(reclustered.path,"/02_dot_topMarkersPhase2_KOWTsepar.svg"), width = 4.3, height = 6)

seurat_obj$Phase2_FinClNamCond <- paste(seurat_obj$Phase2, seurat_obj$FinalClusterNameCondition, sep = "_")
table(seurat_obj$Phase2_FinClNamCond)




### FIRST ANNOTATE CLUSTERS AND THEN RENAME THEM

seurat_obj_reclust <- seurat_obj
# seurat_obj_reclust <- seurat_obj_melan
# seurat_obj_reclust <- seurat_obj_ovar
# seurat_obj_ovar_1_25 <- seurat_obj_reclust
# seurat_obj_melan_1_25 <- seurat_obj_reclust

if (reg.status=="regCC") { 
  print("We don't have a correct definition for clusters with CC genes regressed out")
  
} else if (reg.status=="notRegCC") {
  
  # rename clusters depending on the cancer type
  if (disease[1]=="ovarT") {
    # # for ovarian
    print("Ovarian clusters will be renamed")
    
    seurat_obj_reclust$reclustering <- paste(ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 4, "Tpex1", ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 6, "Tpex2",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 13, "Tpex2",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 5, "Tinex2",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 8, "Tinex2",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 11, "Tinex2",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 0, "Tinex1",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 9, "Tinex1",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 2, "Teffex",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 14, "Teffex",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 1, "Ttex",ifelse(seurat_obj_reclust$integrated_snn_res.1.5 == 10, "Ttex","Tinex2")))))))))))))
    

    seurat_obj_reclust$FinClusterNameCondition2 <- paste(seurat_obj_reclust$reclustering,seurat_obj_reclust$condition,sep = "_")
    Idents(seurat_obj_reclust) <- "FinClusterNameCondition2"
    
    ########## EXTRACT BARCODES and prepare for velocity analysis
    # extract filtered barcodes from the seurat file - this is needed to prepare files for scVelo - analysis of RNA velocity - see the python notebook for more details
    df <- seurat_obj_reclust@meta.data %>%
      rownames_to_column("barcodes") %>%
      select(barcodes)
    
    barcode <-data.frame(df, str_split_fixed(df$barcodes, '_', 2))
    colnames(barcode) <- c("full_barcode","id","barcode")
    
    # table(barcode$id)
    
    for (r in 1:length(sample)){ 
      filt.barcodes <- subset(barcode, id==sample[r])[,3]
      write.table(filt.barcodes, file = paste0(reclustered.path,"/barcodes.filt.CLEAN_BARCODE.",sample[r],"_10-15-24.txt"), sep = "\t", row.names = FALSE, col.names = F,quote=F)
      filt.barcodes <- subset(barcode, id==sample[r])[,1]
      write.table(filt.barcodes, file = paste0(reclustered.path,"/barcodes.filt.FULL_BARCODE.",sample[r],"_10-15-24.txt"), sep = "\t", row.names = FALSE, col.names = F,quote=F)
      
    }
    
    # just this is enough to save file as loom file
    SaveLoom(seurat_obj_reclust, filename = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj_LOOM")), overwrite = TRUE, verbose = TRUE)
    
    
    Idents(seurat_obj_reclust) <- "reclustering"
    seurat_obj_reclust$reclustering <- factor(seurat_obj_reclust$reclustering, levels=c("Tpex1","Tpex2","Tinex1","Tinex2","Teffex","Ttex"))
    
    Idents(seurat_obj_reclust) <- "condition"
    seurat_obj_reclust$condition <- factor(seurat_obj_reclust$condition, levels=c("WT","KO"))
    
    
    seurat_obj_reclust$FinClusterNameCondition2 <- factor(seurat_obj_reclust$FinClusterNameCondition2, levels=c("Tpex1_WT","Tpex1_KO","Tpex2_WT","Tpex2_KO","Tinex1_WT","Tinex1_KO","Tinex2_WT","Tinex2_KO","Teffex_WT","Teffex_KO","Ttex_WT","Ttex_KO"))
    Idents(seurat_obj_reclust) <- "reclustering"
    SaveSeuratRds(seurat_obj_reclust, file = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))
    
    seurat_obj_ovar <- seurat_obj_reclust
    
    
    
    
  } else if (disease[1]=="melanT") {
    print("Melanoma clusters will be renamed")
   
    seurat_obj_reclust$reclustering <- paste(ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 0, "Tpex1", ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 3, "Tinex1",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 1, "Tinex2",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 4, "Tinex2",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 2, "Tpex2",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 5, "Ttex",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 6, "Tinex1",ifelse(seurat_obj_reclust$integrated_snn_res.0.5 == 7, "Tpex1","na")))))))))
    
    seurat_obj_reclust$FinClusterNameCondition2 <- paste(seurat_obj_reclust$reclustering,seurat_obj_reclust$condition,sep = "_")
    Idents(seurat_obj_reclust) <- "FinClusterNameCondition2"

    
    ########## EXTRACT BARCODES and prepare for velocity analysis
    # extract filtered barcodes from the seurat file - this is needed to prepare files for scVelo - analysis of RNA velocity - see the python notebook for more details
    df <- seurat_obj_reclust@meta.data %>%
      rownames_to_column("barcodes") %>%
      select(barcodes)
    
    barcode <-data.frame(df, str_split_fixed(df$barcodes, '_', 2))
    colnames(barcode) <- c("full_barcode","id","barcode")
    
    # table(barcode$id)
    
    for (r in 1:length(sample)){ 
      filt.barcodes <- subset(barcode, id==sample[r])[,3]
      write.table(filt.barcodes, file = paste0(reclustered.path,"/barcodes.filt.CLEAN_BARCODE.",sample[r],"_10-15-24.txt"), sep = "\t", row.names = FALSE, col.names = F,quote=F)
      filt.barcodes <- subset(barcode, id==sample[r])[,1]
      write.table(filt.barcodes, file = paste0(reclustered.path,"/barcodes.filt.FULL_BARCODE.",sample[r],"_10-15-24.txt"), sep = "\t", row.names = FALSE, col.names = F,quote=F)
      
    }
    
    # just this is enough to save file as loom file
    SaveLoom(seurat_obj_reclust, filename = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj_LOOM")), overwrite = TRUE, verbose = TRUE)
    
    
    seurat_obj_reclust$FinClusterNameCondition2 <- factor(seurat_obj_reclust$FinClusterNameCondition2, levels=c("Tpex1_WT","Tpex1_KO","Tpex2_WT","Tpex2_KO","Tinex1_WT","Tinex1_KO","Tinex2_WT","Tinex2_KO","Ttex_WT","Ttex_KO"))
    Idents(seurat_obj_reclust) <- "reclustering"
    
    SaveSeuratRds(seurat_obj_reclust, file = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))
    
    
    seurat_obj_melan <- seurat_obj_reclust
    
                                                
  } else {
    print("Check if correct disease is selected")
  }
  

} else {
  print("Check if correct regression status is selected")
}





# SaveSeuratRds(seurat_obj_reclust, file = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))
# 
# SaveSeuratRds(seurat_obj_reclust, file = file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))
# 

# seurat_obj_melan <- LoadSeuratRds(file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))
# 
# seurat_obj_ovar <- LoadSeuratRds(file.path(integrated.path, paste0("03_reclustered_",extension,"_seurat_obj.Rds")))




if (disease[1]=="ovarT") {
  # # for ovarian
  print("Ovarian data will be plotted")
  
  seurat_obj_reclust <- seurat_obj_ovar
  
  
  ### for OVAR
  # # DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "reclustering",split.by = "condition",reduction = "tsne", label = T, colors_use = c("#00BA38","#00BFC4","#619CFF","#F8766D","#B79F00","#F564E3"))
  # for ovar - updated effex
  DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "reclustering",split.by = "condition",reduction = "tsne", label = T, colors_use = c("#00BA38","#F8766D","#00BFC4","#B79F00","#F564E3","#619CFF"))
  ggsave2(paste0(reclustered.path,"/09_ovar_dimplot_reclust_splitWTKO_filt_resMix_",reg.status,"_",extension,".svg"), width = 8.5, height = 4.5)
  
  
  DotPlot(seurat_obj_reclust, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/09f_ovar_bubble_reclust_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 6, height = 10)
  
  DotPlot(seurat_obj_reclust, features = fcer1g.signature, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/09_ovar_bubble_reclust_fcer1g_WT-KO_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 6, height = 4.8)
  
  
  Stacked_VlnPlot(seurat_object = seurat_obj_reclust, features = c("Cd3e","Cd8a","Tcf7","Bach2","Il7r","Pdcd1", "Ctla4", "Havcr2","Lag3"), x_lab_rotate = TRUE, colors_use = c("white","#ff8686"), group.by = "FinClusterNameCondition2",split.by = "condition")
  ggsave2(paste0(reclustered.path,"/09b_ovar_vlnplot_basicTcellAnnot_split-Cond.svg"), width = 6.5, height =7)
  
  Stacked_VlnPlot(seurat_object = seurat_obj_reclust, features = feature_genes_list, x_lab_rotate = TRUE, colors_use = c("white","#ff8686"), group.by = "FinClusterNameCondition2",split.by = "condition")
  ggsave2(paste0(reclustered.path,"/09g_ovar_vlnplot_featGeneList_split-Cond.svg"), width = 6.5, height =8)
  
  Stacked_VlnPlot(seurat_object = seurat_obj_reclust, features = rev(gene.list[[1]]), x_lab_rotate = TRUE, colors_use = c("white","#ff8686"), group.by = "FinClusterNameCondition2",split.by = "condition")
  ggsave2(paste0(reclustered.path,"/09g_ovar_vlnplot_geneList_split-Cond.svg"), width = 6.2, height =25)
  
  
  
  
  
  ############## Calculate and plot signature scores
  signatures.path <- paste0(reclustered.path,"/signatures")
  if (!file.exists(signatures.path)) {
    dir.create(signatures.path,recursive = T)
  }
  
  working_res <- "FinClusterNameCondition2"
  DefaultAssay(seurat_obj_reclust) <- "RNA"
  Idents(seurat_obj_reclust) <- working_res
  
  for (i in 1:(length(gene.signatures.list))) {
    signatures <- AddModuleScore_UCell(seurat_obj_reclust, features = gene.signatures.list[[i]])
    signature.names <- paste0(names(gene.signatures.list[[i]]), "_UCell")
    
    # generate violin plots
    VlnPlot(signatures, features = signature.names, group.by = working_res, split.by = "diseaseCondition",cols = c("#ff8686","white"),ncol=4)
    ggsave2(paste0(signatures.path,"/signaturesE_ovar_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 26, height = 4.5)
    
    # Stacked_VlnPlot(seurat_object = signatures, features = signature.names, x_lab_rotate = TRUE, colors_use = c("#ff8686","white"), group.by = working_res, split.by = "diseaseCondition")
    # ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 20, height = 5)
    
    # generate feature plots
    # FeaturePlot(signatures, features = signature.names,reduction = "tsne",split.by = "diseaseCondition", ncol = 3, order = T,pt.size=0.9)
    # ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_featplot_split-disCond_",working_res,".svg"), width = 7, height = 9)
    
    FeaturePlot_scCustom(seurat_object = signatures, features = signature.names,reduction = "tsne",split.by = "diseaseCondition")
    ggsave2(paste0(signatures.path,"/signaturesE_ovar_",names(gene.signatures.list[i]),"_featplot_split-disCond_",working_res,".svg"), width = 7, height = 9)
    
  }
  

  
  seurat_obj_reclust$Phase2 <- paste(ifelse(seurat_obj_reclust$Phase == "G1","G1","S/G2/M"))
  seurat_obj_reclust$Phase2_cond <- paste(seurat_obj_reclust$Phase2, seurat_obj_reclust$conditionDisease, sep = "_")
  seurat_obj_reclust$Phase2_cond_cl <- paste(seurat_obj_reclust$Phase2_cond, seurat_obj_reclust$reclustering, sep = "_")
  seurat_obj_reclust$Phase2_cl <- paste(seurat_obj_reclust$Phase2, seurat_obj_reclust$reclustering, sep = "_")
  
  table(seurat_obj_reclust$Phase2)
  # ovarian
  # G1 S/G2/M 
  # 1473   1603 

  table(seurat_obj_reclust$Phase2_cond)
  # ovarian
  # G1_KO     G1_WT S/G2/M_KO S/G2/M_WT 
  # 747       726       613       990 

  table(seurat_obj_reclust$Phase2_cl)
  # ovarian
  # G1_Teffex     G1_Tinex1     G1_Tinex2      G1_Tpex1      G1_Tpex2       G1_Ttex S/G2/M_Teffex S/G2/M_Tinex1 S/G2/M_Tinex2 
  # 53           442           165           269           143           401           285           106           952 
  # S/G2/M_Tpex1  S/G2/M_Tpex2   S/G2/M_Ttex 
  # 19           169            72 
  

  
  # CELL CYCLE PHASE regardless of genotype - combined together
  # ovarian
  data <- data.frame(
    label = c("G1", "S/G2/M"),
    Tpex1 = c(269,19 ),
    Tpex2 = c(143, 169),
    Tinex1 = c(442,106),
    Tinex2 = c(165, 952),
    Teffex = c(53, 285),
    Ttex = c(401, 72)
  )
  data_long <- pivot_longer(data, cols = -label)
  data_long$name <- factor(data_long$name, levels=c("Tpex1","Tpex2","Tinex1","Tinex2","Teffex","Ttex"))
  
  
  p1 <- ggplot(data_long, aes(fill=label, y=value, x=name, width=.7)) + geom_bar(position="fill", stat="identity")+ theme_bw(base_size = 12) + ylab("Cell phase proportion [%]") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=12), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=12), axis.text.x = element_text(color="black",size=12), axis.text.y = element_text(color="black",size=12)) + theme(legend.text = element_text(color="black",size = 12), legend.title = element_text(color="black",size = 12)) + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = c("navy","orange"))
  ggsave(paste0(reclustered.path,"/09_barplot_reclustered_proportions_CellPhase_WTKOtogether_resMix.svg"), p1, width = 5, height = 3)
  
  
  # create % proportion for each cluster FOR CLUSTERS
  Idents(seurat_obj_reclust) <- "condition"
  seurat_obj_reclust_wt = subset(x = seurat_obj_reclust, idents = "KO", invert = T)
  
  Idents(seurat_obj_reclust_wt) = "reclustering"
  wt_cluster_sizes <- table(Idents(seurat_obj_reclust_wt))
  wt_clusters <- names(wt_cluster_sizes)[wt_cluster_sizes > 2]
  wt_cluster_sizes / sum(wt_cluster_sizes)
  
  seurat_obj_reclust_ko = subset(x = seurat_obj_reclust, idents = "WT", invert = T)
  
  Idents(seurat_obj_reclust_ko) = "reclustering"
  ko_cluster_sizes <- table(Idents(seurat_obj_reclust_ko))
  ko_clusters <- names(ko_cluster_sizes)[ko_cluster_sizes > 2]
  ko_cluster_sizes / sum(ko_cluster_sizes)

  # # ovarian
  data <- data.frame(
    label = c("WT", "KO"),
    Tpex1 = c(0.09440559,0.09264706 ),
    Tinex1 = c(0.18123543,0.17426471 ),
    Ttex1 = c(0.12762238, 0.18676471),
    Tpex2 = c(0.09731935, 0.10661765),
    Tinex2 = c(0.39277389, 0.32573529),
    Ttex2 = c(0.10664336, 0.11397059)
  )
  
  data_long <- pivot_longer(data, cols = -label)
  data_long$name <- factor(data_long$name, levels=c("Tpex1","Tpex2","Tinex1","Tinex2","Teffex","Ttex"))
  
  
  
  data_long$label <- factor(data_long$label, levels=c("WT","KO"))
  p1 <- ggplot(data_long, aes(x=name, y = value*100, fill = label)) + geom_bar(position="dodge", stat="identity", width=0.5) + theme_bw(base_size = 12) + ylab("Cluster proportion [%]") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=12), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=12), axis.text.x = element_text(color="black",size=12), axis.text.y = element_text(color="black",size=12)) + theme(legend.text = element_text(color="black",size = 12), legend.title = element_text(color="black",size = 12)) + scale_y_continuous(expand = c(0, 0))+ scale_fill_manual(values = c("#00BFC4","#F8766D"))  # "#7ac5cdff","#ff8686ff"
  
  ggsave(paste0(reclustered.path,"/04_barplot_reclustered_proportions_WT-KO_StemExhaust_Filt_resMix.svg"), p1, width = 5, height = 3)
  
  
  p2<- DotPlot(seurat_obj_reclust, features = rev(gene.list),group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  
  
  graph <- ggarrange(p1, p2, ncol=1, nrow=2, heights=c(1, 6))
  ggsave(paste0(reclustered.path,"/bubble_reclust_withBAR_WT-KO_StemExhaust_Filt_resMix.svg"), graph, width = 7, height = 14)
  
  DotPlot(seurat_obj_reclust, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_WT-KO_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 6, height = 8.8)
  
  
  DotPlot(seurat_obj_reclust_wt, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_WTonly_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 5, height = 8.8)
  
  DotPlot(seurat_obj_ovar, features = c("Mapk1","Mapk8","Mapk14","Mapk6","Mapk3","Mapk9","Mapk7","Cd300c","Cd300lb","Cd300c2","Clec4e","Clec6a","Cd300ld","Tyrobp","Hcst","Syk","Fcer1g", "Akt1", "Sting1", "Ccl2", "Ccr2", "Ly6c2", "Pik3ca","Nfatc1","Nfatc2","Nfatc3","Nfatc4","Nfatc5", "Ptpn6","Lyn","Fyn","Src","Bcl11a","Bcl11b"), group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave(paste0(reclustered.path,"/05_dotplot_ovarian_SIGNALING.svg"), width = 6.5, height = 7.5)
  
  
  
  sc <- c("notRegressCC_KOWT"=seurat_obj_reclust,"notRegressCC_WTonly"=seurat_obj_reclust_wt)
  cl <- c("reclustering")
  
  markers.path <- paste0(reclustered.path,"/markers")
  if (!file.exists(markers.path)) {
    dir.create(markers.path,recursive = T)
  }
  
  for (s in 1:length(sc)) {
    integrated.seurat.object <- sc[[s]]
    # joinLayers is needed when using Seurat v5
    integrated.seurat.object <- JoinLayers(integrated.seurat.object, assay = "RNA")
    working_assay <- "RNA"
    DefaultAssay(integrated.seurat.object) <- working_assay
    
    for (c in 1:1) {
      Idents(integrated.seurat.object) <- cl[c]
      markers <- FindAllMarkers(integrated.seurat.object, assay = working_assay, random.seed = seed, only.pos = TRUE) # only.pos = TRUE
      
      # m <- markers
      # markers <- m
      
      markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 10) %>%
        ungroup() -> top
      
      top = top[!duplicated(top$gene),]
      
      markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>%
        ungroup() -> markers2
      
      #filters out all ribosomal and mitochondrial genes
      
      write.table(top, file = paste0(markers.path,"/topMarkers_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      write.table(markers, file = paste0(markers.path,"/allMarkers_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      write.table(markers2, file = paste0(markers.path,"/allMarkersFiltRBMT_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      DotPlot(integrated.seurat.object, features = top$gene,group.by = cl[c], assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip() +RotatedAxis()
      ggsave(paste0(markers.path,"/dot_topMarkers_",names(sc)[[s]],"_",cl[c],".svg"), width = 7, height = 14)
      
      # genes <- as.data.frame(gene.list[1])
      top$label <- paste0(top$cluster,top$gene)
      markers$label2 <- paste0(markers$cluster,markers$gene)
      markers$label <- markers$label2 %in% top[[8]]
      
      markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
      
      ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val))) + geom_point(aes(color = significant))+ geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
        data = subset(markers, label == TRUE),
        aes(label = gene),
        size = 4,
        box.padding = unit(3, "lines"),
        point.padding = unit(0.25, "lines"),
        max.overlaps = Inf) +  scale_color_manual(values = c("black","grey" )) + ggtitle(paste0(disease," CD8 T cells, ",names(sc)[[s]],"_",cl[c]," WT only"))+facet_wrap(~cluster)
      ggsave(paste0(markers.path,"/volcano_topMarkers_",names(sc)[[s]],"_",cl[c],".svg"), width = 15, height = 10)
      
    }
  }
  
  
  # differential analysis between two conditions
  working_res <- "reclustering"
  seurat_obj_reclust$conditionDisease.clust <- paste(seurat_obj_reclust$conditionDisease, seurat_obj_reclust[[working_res]][,1],sep = "_")
  
  for (k in 1:(length(table(seurat_obj_reclust@meta.data[[working_res]]))))  {
    tryCatch({
      DefaultAssay(seurat_obj_reclust) <- working_assay
      Idents(seurat_obj_reclust) <- "conditionDisease.clust"
      # finds markers within each cluster differential between the two conditions / replicates
      i <- names(table(seurat_obj_reclust[[working_res]])[k])
      
      markers <- FindMarkers(seurat_obj_reclust, assay = working_assay, ident.1 = paste0(condition1,"_",i), ident.2 = paste0(condition2,"_",i), verbose = T, random.seed = seed)
      
      markers$gene <- rownames(markers)
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>%
        ungroup() -> markers
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 20) %>%
        ungroup() -> top
      
      write.table(markers, file = paste0(markers.path,"/cl_",i,"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      top = top[!duplicated(top$gene),]
      
      for (g in 1:length(gene.list)){
        genes <- as.data.frame(gene.list[g])
        markers$label <- markers$gene %in% genes[,1]
        
        
        markers$label <- markers$gene %in% top[[6]]
        
        
        
        markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
        
        
        ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
          data = subset(markers, label == TRUE),
          aes(label = gene),
          size = 5,
          box.padding = unit(3, "lines"),
          point.padding = unit(0.25, "lines"),
          max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","black","#7ac5cdff" )) + ggtitle(paste0(disease," CD8 T cells, cluster ",i," ",condition1,"/",condition2))
        
        ggsave(paste0(markers.path,"/volcano_cl",i,"_",names(gene.list[g]),"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
      }
      
      
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  
  
  
  # this was used to identify commonly deregulated genes in ovarian and melanoma and also for integration with epigenetics data...this is simple comparison between KO and WT regardless of any cluster
  for (k in 1:(length(table(seurat_obj_reclust$disease))))  {
    tryCatch({
      DefaultAssay(seurat_obj_reclust) <- working_assay
      Idents(seurat_obj_reclust) <- "conditionDisease"
      # finds markers within each cluster differential between the two conditions / replicates
      i <- names(table(seurat_obj_reclust$disease))[[k]]
      
      markers <- FindMarkers(seurat_obj_reclust, assay = working_assay, ident.1 = names(table(seurat_obj_reclust$conditionDisease))[[1]], ident.2 = names(table(seurat_obj_reclust$conditionDisease))[[2]], verbose = T, random.seed = seed,logfc.threshold = 0)
      
      markers$gene <- rownames(markers)
      
      markers %>% dplyr::filter(p_val < 0.05) %>% dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% 
        ungroup() -> markers_filtered
      
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 20) %>%
        ungroup() -> top
      
      write.table(markers_filtered, file = paste0(markers.path,"/allKOvsWT_",i,"_assay",working_assay,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      # the problem is that some gene symbols were presented twice (different ensembl ID) in the original features dataset
      # seurat deals with that and makes the gene symbols unique by adding .1 which can be verified for example on genes Aldoa, Ddit3, Nnt, Nrg1 by looking at them in exported all genes from the seurat object
      # https://github.com/satijalab/seurat/issues/978
      # t <- as.data.frame(rownames(seurat_obj@assays$RNA@counts))
      
      ######### Alternatively CONSIDER PLAYING WITH ReadMtx function of seurat to use the unique cell names (ensembl id can be set by changing value 2 to 1 in feature.column = 1...or consider removing the duplicates or labeling them somehow
      # https://github.com/satijalab/seurat/issues/1867
      
      # add ensembl ID extracted from features.tsv file (it contains the same genes for wt and ko samples so it can be done only once)
      Project_Name <- Project_Names[1]
      Feature_File <- paste0(input_path,"/",Project_Name,"/",Features_name)
      Feature <- read.csv2(Feature_File, sep = "", header = F)
      Feature <- Feature[,1:2]
      colnames(Feature) <- c("geneID","geneSymbol_with_dupl")
      Feature$gene <- make.unique(Feature$geneSymbol_with_dupl)
      
      markers_annotated <- merge(Feature, markers,by="gene")
      
      write.table(markers_annotated, file = paste0(markers.path,"/01-31-2025_allKOvsWT_ANNOTATED_",i,"_assay",working_assay,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      ### THESE WERE FURTHER USED FOR GSEA BY A SEPARATE SCRIPT in /Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN/2-Single_Cell_RNAseq_Pipeline/Output/Rpathways 
      
      
      top = top[!duplicated(top$gene),]
      
      for (g in 1:length(gene.list)){
        genes <- as.data.frame(gene.list[g])
        markers_filtered$label <- markers_filtered$gene %in% genes[,1]
        markers_filtered$label <- markers_filtered$gene %in% top[[6]]
        
        markers_filtered$significant <- ifelse(markers_filtered$label == TRUE, "Highlight", ifelse(markers_filtered$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
        
        ggplot(markers_filtered %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
          data = subset(markers_filtered, label == TRUE),
          aes(label = gene),
          size = 5,
          box.padding = unit(3, "lines"),
          point.padding = unit(0.25, "lines"),
          max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","black","#7ac5cdff" )) + ggtitle(paste0(i," CD8 TILs, KO vs WT"))
        
        ggsave(paste0(markers.path,"/volcano_allKOvsWT_",i,"_",names(gene.list[g]),"_assay",working_assay,".svg"), width = 8, height = 4)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  
  
  
    # clustered dot plot
  svg(paste0(reclustered.path,"/06_clustered_dotplot_ovarian.svg"), width = 6.5, height = 7.5, family ="Arial")
  Clustered_DotPlot(seurat_obj_reclust,features = gene.list[[1]], group.by = "FinClusterNameCondition2", assay=working_assay, cluster_ident = F)
  dev.off()
  
  # Clustered_DotPlot(seurat_obj_reclust,features = fcer1g.signature, group.by = "FinClusterNameCondition2", assay=working_assay, cluster_ident = F)
  
  

  # plot density and feature plots for individual genes
  # https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html
  
  
  for (i in 1:(length(feature_genes_list))) {
    # i <- 1
    feature <- feature_genes_list[i]
    
    p1= Plot_Density_Custom(seurat_object = seurat_obj_reclust_wt, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in WT"))
    p2= Plot_Density_Custom(seurat_object = seurat_obj_reclust_ko, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in KO"))
    p1+p2
    ggsave2(paste0(signatures.path,"/featplot_",feature,"_split-disCond_",working_res,".svg"), width = 7, height = 3)
    
    p1= FeaturePlot_scCustom(seurat_object = seurat_obj_reclust_wt, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in WT"))
    p2= FeaturePlot_scCustom(seurat_object = seurat_obj_reclust_ko, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in KO"))
    plot_grid(p1, p2)
    ggsave2(paste0(signatures.path,"/featplotNOTDENSITY_",feature,"_split-disCond_",working_res,".svg"), width = 7, height = 3)
    
  }
  
  
  DotPlot(seurat_obj_reclust_wt, features = annot_wt, group.by = "integrated_snn_res.0.5", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/05_bubble_WTonlyAnnot_filt_res05_",reg.status,"_",extension,".svg"), width = 4.5, height = 5.5)
  
  DotPlot(seurat_obj_reclust_wt, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/05_bubble_reclust_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4.5, height = 5.5)
  
  
  Idents(seurat_obj_reclust_wt) <- "reclustering"
  obj = subset(x = seurat_obj_reclust_wt, idents = "Tpex2", invert = T)
  obj = subset(x = obj, idents = "Tinex2", invert = T)
  obj = subset(x = obj, idents = "Teffex", invert = T)

  DotPlot(obj, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_G1_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4, height = 6)
  
  obj = subset(x = seurat_obj_reclust_wt, idents = "Tpex1", invert = T)
  obj = subset(x = obj, idents = "Tinex1", invert = T)
  obj = subset(x = obj, idents = "Ttex", invert = T)
  DotPlot(obj, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_SG2M_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4, height = 6)
  
  
  
  
  
  # # this creates heatmap including all cells
  # DoHeatmap(seurat_obj_reclust, features = annot_wt, assay=working_assay,group.by = "reclustering")
  
  # prepare for heatmap visualization using average expression values per ident (cluster or genotype)
  # As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
  # obj.averages <- AverageExpression(seurat_obj_reclust, return.seurat = TRUE)
  
  obj.aggregated <- AggregateExpression(
    seurat_obj_reclust,
    assays = working_assay,
    return.seurat = TRUE, group.by = "FinClusterNameCondition2"
  )
  
  DoHeatmap(obj.aggregated, features = annot_wt, label = T ,draw.lines = F, size =2)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_annotWT_clustCond.svg"), width = 4.5, height = 3)
  
  DoHeatmap(obj.aggregated, features = rev(gene.list[[1]]), label = T ,draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_geneList_clustCond.svg"), width = 5, height = 5.5)
  
  
  obj.aggregated <- AggregateExpression(
    seurat_obj_reclust,
    assays = working_assay,
    return.seurat = TRUE, group.by = "reclustering"
  )
  DoHeatmap(obj.aggregated, features = annot_wt, label = T ,draw.lines = F)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_annotWT_reclust.svg"), width = 4.5, height = 3)
  
  
  DoHeatmap(obj.aggregated, features = rev(gene.list[[1]]), label = T ,draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_geneList_reclust.svg"), width = 5, height = 5.5)
  
  rm(obj.aggregated)
  
  
  
  
  
  
  
  
} else if (disease[1]=="melanT") {
  print("Melanoma data will be plotted")
  
  seurat_obj_reclust <- seurat_obj_melan
  
  # for melanoma USE
  DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "reclustering",split.by = "condition",reduction = "tsne", label = T, colors_use = c("#00BFC4","#B79F00","#00BA38","#F8766D","#619CFF"))
  ggsave2(paste0(reclustered.path,"/04_dimplot_reclust_splitWTKO_filt_resMix_",reg.status,"_",extension,".svg"), width = 8.5, height = 4.5)
  
  
  # DimPlot(seurat_obj_reclust, label = T, repel = T,reduction = "tsne", group.by = "Phase") + ggtitle("Cell cycle phase")
  # ggsave2(paste0(reclustered.path,"/04_dimplot_reclust_CCphase_filt_resMix_",reg.status,"_",extension,".svg"), width = 4.5, height = 4.5)
  
  DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "reclustering",reduction = "tsne", label = T, colors_use = c("#00BFC4","#B79F00","#00BA38","#F8766D","#619CFF"))
  ggsave2(paste0(reclustered.path,"/04_dimplot_noSplit_filt_resMix_",reg.status,"_",extension,".svg"), width = 5, height = 5)
  
  DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "condition",reduction = "tsne", label = T, colors_use = c("#00BFC4","#F8766D"))
  ggsave2(paste0(reclustered.path,"/04_dimplot_WTKOboth_filt_",reg.status,"_",extension,".svg"), width = 5, height = 5)
  
  
  DimPlot_scCustom(seurat_object = seurat_obj_reclust, group.by = "sample",reduction = "tsne", label = T, colors_use = c("#00BFC4","#619CFF","#F8766D","#B79F00"))
  ggsave2(paste0(reclustered.path,"/04_dimplot_Samples_filt_resMix_",reg.status,"_",extension,".svg"), width = 5, height = 5)
  
  
  
  Stacked_VlnPlot(seurat_object = seurat_obj_reclust, features = feature_genes_list, x_lab_rotate = TRUE, colors_use = c("white","#ff8686"), group.by = "FinClusterNameCondition2",split.by = "condition")
  ggsave2(paste0(reclustered.path,"/09g_melan_vlnplot_featGeneList_split-Cond.svg"), width = 6, height =8)
  
  Stacked_VlnPlot(seurat_object = seurat_obj_reclust, features = rev(gene.list[[1]]), x_lab_rotate = TRUE, colors_use = c("white","#ff8686"), group.by = "FinClusterNameCondition2",split.by = "condition")
  ggsave2(paste0(reclustered.path,"/09g_melan_vlnplot_geneList_split-Cond.svg"), width = 6, height =25)
  
  
  
  DotPlot(seurat_obj_reclust, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/09f_melan_bubble_reclust_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 5.5, height = 10)
  
  
  ############## Calculate and plot signature scores
  signatures.path <- paste0(reclustered.path,"/signatures")
  if (!file.exists(signatures.path)) {
    dir.create(signatures.path,recursive = T)
  }
  
  working_res <- "FinClusterNameCondition2"
  DefaultAssay(seurat_obj_reclust) <- "RNA"
  Idents(seurat_obj_reclust) <- working_res
  
  for (i in 1:(length(gene.signatures.list))) {
    # i=1
    signatures <- AddModuleScore_UCell(seurat_obj_reclust, features = gene.signatures.list[[i]])
    signature.names <- paste0(names(gene.signatures.list[[i]]), "_UCell")
    
    # generate violin plots
    VlnPlot(signatures, features = signature.names, group.by = working_res, split.by = "diseaseCondition",cols = c("#ff8686","white"),ncol=4)
    ggsave2(paste0(signatures.path,"/signaturesE_melan_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 26, height = 4.5)
    
    # Stacked_VlnPlot(seurat_object = signatures, features = signature.names, x_lab_rotate = TRUE, colors_use = c("#ff8686","white"), group.by = working_res, split.by = "diseaseCondition")
    # ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 20, height = 5)
    
    # generate feature plots
    # FeaturePlot(signatures, features = signature.names,reduction = "tsne",split.by = "diseaseCondition", ncol = 3, order = T,pt.size=0.9)
    # ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_featplot_split-disCond_",working_res,".svg"), width = 7, height = 9)
    
    FeaturePlot_scCustom(seurat_object = signatures, features = signature.names,reduction = "tsne",split.by = "diseaseCondition")
    ggsave2(paste0(signatures.path,"/signaturesE_melan_",names(gene.signatures.list[i]),"_featplot_split-disCond_",working_res,".svg"), width = 7, height = 9)
    
  }
  
  
  seurat_obj_reclust$Phase2 <- paste(ifelse(seurat_obj_reclust$Phase == "G1","G1","S/G2/M"))
  seurat_obj_reclust$Phase2_cond <- paste(seurat_obj_reclust$Phase2, seurat_obj_reclust$conditionDisease, sep = "_")
  seurat_obj_reclust$Phase2_cond_cl <- paste(seurat_obj_reclust$Phase2_cond, seurat_obj_reclust$reclustering, sep = "_")
  seurat_obj_reclust$Phase2_cl <- paste(seurat_obj_reclust$Phase2, seurat_obj_reclust$reclustering, sep = "_")
  
  table(seurat_obj_reclust$Phase2)
  # # melanoma 
  # G1 S/G2/M 
  # 1788   1598 
  
  table(seurat_obj_reclust$Phase2_cond)
  # G1_KOmelanT     G1_WTmelanT S/G2/M_KOmelanT S/G2/M_WTmelanT 
  # 359            1429             281            1317
  
  table(seurat_obj_reclust$Phase2_cl)

  # # melanoma 
  # G1_Tinex1     G1_Tinex2      G1_Tpex1      G1_Tpex2       G1_Ttex S/G2/M_Tinex1 S/G2/M_Tinex2  S/G2/M_Tpex1 
  # 564           228           495           392           109           142           745           429 
  # S/G2/M_Tpex2   S/G2/M_Ttex 
  # 220            62 
  
  
  
  # CELL CYCLE PHASE regardless of genotype - combined together
  # melanoma
  data <- data.frame(
    label = c("G1", "S/G2/M"),
    Tpex1 = c(495,429),
    Tpex2 = c(392,220),
    Tinex1 = c(564,142),
    Tinex2 = c(228,745),
    Ttex = c(109,62)
  )
  data_long <- pivot_longer(data, cols = -label)
  data_long$name <- factor(data_long$name, levels=c("Tpex1","Tpex2","Tinex1","Tinex2","Ttex"))
  
  p1 <- ggplot(data_long, aes(fill=label, y=value, x=name, width=.7)) + geom_bar(position="fill", stat="identity")+ theme_bw(base_size = 12) + ylab("Cell phase proportion [%]") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=12), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=12), axis.text.x = element_text(color="black",size=12), axis.text.y = element_text(color="black",size=12)) + theme(legend.text = element_text(color="black",size = 12), legend.title = element_text(color="black",size = 12)) + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = c("navy","orange"))
  ggsave(paste0(reclustered.path,"/09_barplot_reclustered_proportions_CellPhase_WTKOtogether_resMix.svg"), p1, width = 5, height = 3)
  
  
  
  # create % proportion for each cluster FOR CLUSTERS
  Idents(seurat_obj_reclust) <- "condition"
  seurat_obj_reclust_wt = subset(x = seurat_obj_reclust, idents = "KO", invert = T)
  
  Idents(seurat_obj_reclust_wt) = "reclustering"
  wt_cluster_sizes <- table(Idents(seurat_obj_reclust_wt))
  wt_clusters <- names(wt_cluster_sizes)[wt_cluster_sizes > 2]
  wt_cluster_sizes / sum(wt_cluster_sizes)
  
  seurat_obj_reclust_ko = subset(x = seurat_obj_reclust, idents = "WT", invert = T)
  
  Idents(seurat_obj_reclust_ko) = "reclustering"
  ko_cluster_sizes <- table(Idents(seurat_obj_reclust_ko))
  ko_clusters <- names(ko_cluster_sizes)[ko_cluster_sizes > 2]
  ko_cluster_sizes / sum(ko_cluster_sizes)
  
  # # melanoma
  data <- data.frame(
    label = c("WT", "KO"),
    Tpex1 = c(0.27312454,0.2718750),
    Tinex1 = c(0.20648216, 0.2171875),
    Ttex = c(0.04260743, 0.0843750),
    Tpex2 = c(0.18536052,0.1609375 ),
    Tinex2 = c(0.29242535,0.2656250 )
  )
  data_long <- pivot_longer(data, cols = -label)
  data_long$name <- factor(data_long$name, levels=c("Tpex1","Tpex2","Tinex1","Tinex2","Ttex"))
  

  
  data_long$label <- factor(data_long$label, levels=c("WT","KO"))
  p1 <- ggplot(data_long, aes(x=name, y = value*100, fill = label)) + geom_bar(position="dodge", stat="identity", width=0.5) + theme_bw(base_size = 12) + ylab("Cluster proportion [%]") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=12), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=12), axis.text.x = element_text(color="black",size=12), axis.text.y = element_text(color="black",size=12)) + theme(legend.text = element_text(color="black",size = 12), legend.title = element_text(color="black",size = 12)) + scale_y_continuous(expand = c(0, 0))+ scale_fill_manual(values = c("#00BFC4","#F8766D"))  # "#7ac5cdff","#ff8686ff"
  
  ggsave(paste0(reclustered.path,"/04_barplot_reclustered_proportions_WT-KO_StemExhaust_Filt_resMix.svg"), p1, width = 5, height = 3)
  
  
  p2<- DotPlot(seurat_obj_reclust, features = rev(gene.list),group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  
  
  graph <- ggarrange(p1, p2, ncol=1, nrow=2, heights=c(1, 6))
  ggsave(paste0(reclustered.path,"/bubble_reclust_withBAR_WT-KO_StemExhaust_Filt_resMix.svg"), graph, width = 7, height = 14)
  
  DotPlot(seurat_obj_reclust, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_WT-KO_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 6, height = 8.8)
  
  
  DotPlot(seurat_obj_reclust_wt, features = gene.list, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_WTonly_geneList_filt_resMix_",reg.status,"_",extension,".svg"), width = 5, height = 8.8)
  

  DotPlot(seurat_obj_reclust, features = c("Mapk1","Mapk8","Mapk14","Mapk6","Mapk3","Mapk9","Mapk7","Cd300c","Cd300lb","Cd300c2","Clec4e","Clec6a","Cd300ld","Tyrobp","Hcst","Syk","Fcer1g", "Akt1", "Sting1", "Ccl2", "Ccr2", "Ly6c2", "Pik3ca","Nfatc1","Nfatc2","Nfatc3","Nfatc4","Nfatc5", "Ptpn6","Lyn","Fyn","Src","Bcl11a","Bcl11b", "Fasl","Tnfsf10"), group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave(paste0(reclustered.path,"/05_dotplot_melanoma_SIGNALING.svg"), width = 6, height = 7.5)
  
  
  
  
  sc <- c("notRegressCC_KOWT"=seurat_obj_reclust,"notRegressCC_WTonly"=seurat_obj_reclust_wt)
  cl <- c("reclustering")
  
  markers.path <- paste0(reclustered.path,"/markers")
  if (!file.exists(markers.path)) {
    dir.create(markers.path,recursive = T)
  }
  
  for (s in 1:length(sc)) {
    integrated.seurat.object <- sc[[s]]
    # joinLayers is needed when using Seurat v5
    integrated.seurat.object <- JoinLayers(integrated.seurat.object, assay = "RNA")
    working_assay <- "RNA"
    DefaultAssay(integrated.seurat.object) <- working_assay
    
    for (c in 1:1) {
      Idents(integrated.seurat.object) <- cl[c]
      markers <- FindAllMarkers(integrated.seurat.object, assay = working_assay, random.seed = seed, only.pos = TRUE) # only.pos = TRUE
      
      # m <- markers
      # markers <- m
      
      markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 10) %>%
        ungroup() -> top
      
      top = top[!duplicated(top$gene),]
      
      markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>%
        ungroup() -> markers2
      
      #filters out all ribosomal and mitochondrial genes
      
      write.table(top, file = paste0(markers.path,"/topMarkers_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      write.table(markers, file = paste0(markers.path,"/allMarkers_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      write.table(markers2, file = paste0(markers.path,"/allMarkersFiltRBMT_WT_RNA_",names(sc)[[s]],"_",cl[c],".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      DotPlot(integrated.seurat.object, features = top$gene,group.by = cl[c], assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip() +RotatedAxis()
      ggsave(paste0(markers.path,"/dot_topMarkers_",names(sc)[[s]],"_",cl[c],".svg"), width = 7, height = 14)
      
      # genes <- as.data.frame(gene.list[1])
      top$label <- paste0(top$cluster,top$gene)
      markers$label2 <- paste0(markers$cluster,markers$gene)
      markers$label <- markers$label2 %in% top[[8]]
      
      markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
      
      ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val))) + geom_point(aes(color = significant))+ geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
        data = subset(markers, label == TRUE),
        aes(label = gene),
        size = 4,
        box.padding = unit(3, "lines"),
        point.padding = unit(0.25, "lines"),
        max.overlaps = Inf) +  scale_color_manual(values = c("black","grey" )) + ggtitle(paste0(disease," CD8 T cells, ",names(sc)[[s]],"_",cl[c]," WT only"))+facet_wrap(~cluster)
      ggsave(paste0(markers.path,"/volcano_topMarkers_",names(sc)[[s]],"_",cl[c],".svg"), width = 15, height = 10)
      
    }
  }
  
  
  # differential analysis between two conditions per cluster
  working_res <- "reclustering"
  seurat_obj_reclust$conditionDisease.clust <- paste(seurat_obj_reclust$conditionDisease, seurat_obj_reclust[[working_res]][,1],sep = "_")
  
  for (k in 1:(length(table(seurat_obj_reclust@meta.data[[working_res]]))))  {
    tryCatch({
      DefaultAssay(seurat_obj_reclust) <- working_assay
      Idents(seurat_obj_reclust) <- "conditionDisease.clust"
      # finds markers within each cluster differential between the two conditions / replicates
      i <- names(table(seurat_obj_reclust[[working_res]])[k])
      
      markers <- FindMarkers(seurat_obj_reclust, assay = working_assay, ident.1 = paste0(condition1,"_",i), ident.2 = paste0(condition2,"_",i), verbose = T, random.seed = seed)
      
      markers$gene <- rownames(markers)
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>%
        ungroup() -> markers
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 20) %>%
        ungroup() -> top
      
      write.table(markers, file = paste0(markers.path,"/cl_",i,"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      top = top[!duplicated(top$gene),]
      
      for (g in 1:length(gene.list)){
        genes <- as.data.frame(gene.list[g])
        markers$label <- markers$gene %in% genes[,1]
        
        
        markers$label <- markers$gene %in% top[[6]]
        
        
        
        markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
        
        
        ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
          data = subset(markers, label == TRUE),
          aes(label = gene),
          size = 5,
          box.padding = unit(3, "lines"),
          point.padding = unit(0.25, "lines"),
          max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","black","#7ac5cdff" )) + ggtitle(paste0(disease," CD8 T cells, cluster ",i," ",condition1,"/",condition2))
        
        ggsave(paste0(markers.path,"/volcano_cl",i,"_",names(gene.list[g]),"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  
  # this was used to identify commonly deregulated genes in ovarian and melanoma and also for integration with epigenetics data...this is simple comparison between KO and WT regardless of any cluster
  for (k in 1:(length(table(seurat_obj_reclust$disease))))  {
    tryCatch({
      DefaultAssay(seurat_obj_reclust) <- working_assay
      Idents(seurat_obj_reclust) <- "conditionDisease"
      # finds markers within each cluster differential between the two conditions / replicates
      i <- names(table(seurat_obj_reclust$disease))[[k]]
      
      markers <- FindMarkers(seurat_obj_reclust, assay = working_assay, ident.1 = names(table(seurat_obj_reclust$conditionDisease))[[1]], ident.2 = names(table(seurat_obj_reclust$conditionDisease))[[2]], verbose = T, random.seed = seed,logfc.threshold = 0)
      
      markers$gene <- rownames(markers)
      
      markers %>% dplyr::filter(p_val < 0.05) %>% dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% 
        ungroup() -> markers_filtered
      
      
      markers %>%
        dplyr::filter(p_val < 0.05) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene)) %>% slice_head(n = 20) %>%
        ungroup() -> top
      
      write.table(markers_filtered, file = paste0(markers.path,"/allKOvsWT_",i,"_assay",working_assay,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      # the problem is that some gene symbols were presented twice (different ensembl ID) in the original features dataset
      # seurat deals with that and makes the gene symbols unique by adding .1 which can be verified for example on genes Aldoa, Ddit3, Nnt, Nrg1 by looking at them in exported all genes from the seurat object
      # https://github.com/satijalab/seurat/issues/978
      # t <- as.data.frame(rownames(seurat_obj@assays$RNA@counts))
      
      ######### Alternatively CONSIDER PLAYING WITH ReadMtx function of seurat to use the unique cell names (ensembl id can be set by changing value 2 to 1 in feature.column = 1...or consider removing the duplicates or labeling them somehow
      # https://github.com/satijalab/seurat/issues/1867
      
      # add ensembl ID extracted from features.tsv file (it contains the same genes for wt and ko samples so it can be done only once)
      Project_Name <- Project_Names[1]
      Feature_File <- paste0(input_path,"/",Project_Name,"/",Features_name)
      Feature <- read.csv2(Feature_File, sep = "", header = F)
      Feature <- Feature[,1:2]
      colnames(Feature) <- c("geneID","geneSymbol_with_dupl")
      Feature$gene <- make.unique(Feature$geneSymbol_with_dupl)
      
      markers_annotated <- merge(Feature, markers,by="gene")
      
      write.table(markers_annotated, file = paste0(markers.path,"/01-31-2025_allKOvsWT_ANNOTATED_",i,"_assay",working_assay,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      ### THESE WERE FURTHER USED FOR GSEA BY A SEPARATE SCRIPT in /Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN/2-Single_Cell_RNAseq_Pipeline/Output/Rpathways 
      
      
      top = top[!duplicated(top$gene),]
      
      for (g in 1:length(gene.list)){
        genes <- as.data.frame(gene.list[g])
        markers_filtered$label <- markers_filtered$gene %in% genes[,1]
        markers_filtered$label <- markers_filtered$gene %in% top[[6]]
        
        markers_filtered$significant <- ifelse(markers_filtered$label == TRUE, "Highlight", ifelse(markers_filtered$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
        
        ggplot(markers_filtered %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
          data = subset(markers_filtered, label == TRUE),
          aes(label = gene),
          size = 5,
          box.padding = unit(3, "lines"),
          point.padding = unit(0.25, "lines"),
          max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","black","#7ac5cdff" )) + ggtitle(paste0(i," CD8 TILs, KO vs WT"))
        
        ggsave(paste0(markers.path,"/volcano_allKOvsWT_",i,"_",names(gene.list[g]),"_assay",working_assay,".svg"), width = 8, height = 4)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  
  
  
  
  
  
  # clustered dot plot
  svg(paste0(reclustered.path,"/06_clustered_dotplot_ovarian.svg"), width = 6.5, height = 7.5, family ="Arial")
  Clustered_DotPlot(seurat_obj_reclust,features = gene.list[[1]], group.by = "FinClusterNameCondition2", assay=working_assay, cluster_ident = F)
  dev.off()
  
  # Clustered_DotPlot(seurat_obj_reclust,features = fcer1g.signature, group.by = "FinClusterNameCondition2", assay=working_assay, cluster_ident = F)
  
  
  # plot density and feature plots for individual genes
  # https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html

  for (i in 1:(length(feature_genes_list))) {
    # i <- 1
    feature <- feature_genes_list[i]
    
    p1= Plot_Density_Custom(seurat_object = seurat_obj_reclust_wt, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in WT"))
    p2= Plot_Density_Custom(seurat_object = seurat_obj_reclust_ko, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in KO"))
    p1+p2
    ggsave2(paste0(signatures.path,"/featplot_",feature,"_split-disCond_",working_res,".svg"), width = 7, height = 3)
    
    p1= FeaturePlot_scCustom(seurat_object = seurat_obj_reclust_wt, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in WT"))
    p2= FeaturePlot_scCustom(seurat_object = seurat_obj_reclust_ko, features = feature,reduction = "tsne")+ggtitle(paste0(feature," in KO"))
    plot_grid(p1, p2)
    ggsave2(paste0(signatures.path,"/featplotNOTDENSITY_",feature,"_split-disCond_",working_res,".svg"), width = 7, height = 3)
    
  }
  
  
  
  DotPlot(seurat_obj_reclust_wt, features = annot_wt, group.by = "integrated_snn_res.0.5", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/05_bubble_WTonlyAnnot_filt_res05_",reg.status,"_",extension,".svg"), width = 4.5, height = 5.5)
  
  DotPlot(seurat_obj_reclust_wt, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/05_bubble_reclust_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4.5, height = 5.5)
  
  
  Idents(seurat_obj_reclust_wt) <- "reclustering"
  obj = subset(x = seurat_obj_reclust_wt, idents = "Tpex2", invert = T)
  obj = subset(x = obj, idents = "Tinex2", invert = T)
  DotPlot(obj, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_G1_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4, height = 6)
  
  obj = subset(x = seurat_obj_reclust_wt, idents = "Tpex1", invert = T)
  obj = subset(x = obj, idents = "Tinex1", invert = T)
  obj = subset(x = obj, idents = "Ttex", invert = T)
  DotPlot(obj, features = annot_wt, group.by = "FinClusterNameCondition2", assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")
  ggsave2(paste0(reclustered.path,"/04_bubble_reclust_SG2M_WTonlyAnnot_filt_resMix_",reg.status,"_",extension,".svg"), width = 4, height = 6)
  
  
  
  
  
  # # this creates heatmap including all cells
  # DoHeatmap(seurat_obj_reclust, features = annot_wt, assay=working_assay,group.by = "reclustering")
  
  # prepare for heatmap visualization using average expression values per ident (cluster or genotype)
  # As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
  # obj.averages <- AverageExpression(seurat_obj_reclust, return.seurat = TRUE)
  
  obj.aggregated <- AggregateExpression(
    seurat_obj_reclust,
    assays = working_assay,
    return.seurat = TRUE, group.by = "FinClusterNameCondition2"
  )
  
  DoHeatmap(obj.aggregated, features = annot_wt, label = T ,draw.lines = F, size =2)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_annotWT_clustCond.svg"), width = 4.5, height = 3)
  
  DoHeatmap(obj.aggregated, features = rev(gene.list[[1]]), label = T ,draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_geneList_clustCond.svg"), width = 5, height = 5.5)
  
  
  obj.aggregated <- AggregateExpression(
    seurat_obj_reclust,
    assays = working_assay,
    return.seurat = TRUE, group.by = "reclustering"
  )
  DoHeatmap(obj.aggregated, features = annot_wt, label = T ,draw.lines = F)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_annotWT_reclust.svg"), width = 4.5, height = 3)
  
  
  DoHeatmap(obj.aggregated, features = rev(gene.list[[1]]), label = T ,draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
  ggsave(paste0(reclustered.path,"/06_aggregatedExpr_heatmap_geneList_reclust.svg"), width = 5, height = 5.5)
  
  rm(obj.aggregated)
  
  
  
  
  
} else {
  print("Check if correct disease is selected")
}










## NEXT WE COMPARE EXPRESSION BETWEEN MELANOMA AND OVARIAN

    ### load and plot markers combined for melanoma and ovarian models, separately for each cluster
    markers_aggregated <- as.data.frame(read_delim(paste0(markers.path,"/markers_melanT-ovarT_KOvsWT_clusters.txt"), delim = '\t', col_names = T))
    markers_aggregated$rank <- sign(as.numeric(markers_aggregated$log2FoldChange)) * (-log10(as.numeric(markers_aggregated$pvalue)))
    library(reshape)
    library(reshape2)
    markers_aggregated_casted = dcast(markers_aggregated, geneSymbol+cluster+description~disease, mean,  value.var = c("rank")) 
    markers_aggregated_casted <- na.omit(markers_aggregated_casted) #remove NAs
    
    markers_aggregated_casted <- markers_aggregated_casted[order(-markers_aggregated_casted$melanT),]
    
    markers_aggregated_casted$geneCl <- paste0(markers_aggregated_casted$geneSymbol,"_",markers_aggregated_casted$cluster)
    
    ## select gene sets to label for each clusters
    label = markers_aggregated_casted %>%  
      group_by( cluster ) %>% 
      top_n(wt = abs(ovarT), n = 40)
    label$geneCl <- paste0(label$geneSymbol,"_",label$cluster)
    
    markers_aggregated_casted$label2 <- markers_aggregated_casted$geneCl %in% label$geneCl 
    
    markers_aggregated_casted$label <- ifelse(markers_aggregated_casted$label2 == TRUE, ifelse(markers_aggregated_casted$melanT < 0, ifelse(markers_aggregated_casted$ovarT < 0, TRUE, FALSE), ifelse(markers_aggregated_casted$ovarT > 0, TRUE,FALSE)),FALSE)
    
    markers_aggregated_casted$significant <- ifelse(markers_aggregated_casted$label == TRUE, "Highlight", ifelse(markers_aggregated_casted$melanT < 0, ifelse(markers_aggregated_casted$ovarT < 0, "Decreased BOTH","Decreased Melanoma"), ifelse(markers_aggregated_casted$ovarT > 0, "Increased BOTH","Increased Melanoma")))
    
    ggplot(markers_aggregated_casted %>% arrange(label), aes(x = ovarT, y = melanT)) +facet_wrap(~cluster) + geom_hline(yintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + geom_point(aes(color = significant)) + theme_bw(base_size = 12)  + xlab("Ovarian CD8 TILs expression KO/WT") + ylab("Melanoma CD8 TILs expression KO/WT") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
      data = subset(markers_aggregated_casted, label == TRUE),
      aes(label = geneSymbol),
      size = 5,
      box.padding = unit(3, "lines"),
      point.padding = unit(0.25, "lines"),
      max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","grey","black","#7ac5cdff","grey" )) + ggtitle(paste0("KO vs WT RNA expression in CD8 TILs"))
ggsave(paste0(reclustered.path,"/10_volcano_Clusters_melanTvsovarT_topGenes.svg"), width = 14, height = 10)
    




### load and plot markers combined for melanoma and ovarian models, combined for ALL clusters merged
markers_aggregated <- as.data.frame(read_delim(paste0(markers.path,"/markers_melanT-ovarT_KOvsWT_all.txt"), delim = '\t', col_names = T))
markers_aggregated$rank <- sign(as.numeric(markers_aggregated$log2FoldChange)) * (-log10(as.numeric(markers_aggregated$pvalue)))

markers_aggregated_casted = dcast(markers_aggregated, geneSymbol+cluster+description~disease, mean,  value.var = c("rank")) 



markers_aggregated_casted <- na.omit(markers_aggregated_casted) #remove NAs

markers_aggregated_casted <- markers_aggregated_casted[order(-markers_aggregated_casted$melanT),]

markers_aggregated_casted$geneCl <- paste0(markers_aggregated_casted$geneSymbol,"_",markers_aggregated_casted$cluster)

markers_aggregated_casted$sumAbsVal <- as.numeric(paste(abs(markers_aggregated_casted$melanT) + abs(markers_aggregated_casted$ovarT)))


## select gene sets to label for each clusters
label = markers_aggregated_casted %>%  
  group_by( cluster ) %>% 
  top_n(wt = sumAbsVal, n = 70)
label$geneCl <- paste0(label$geneSymbol,"_",label$cluster)


markers_aggregated_casted$label2 <- markers_aggregated_casted$geneCl %in% label$geneCl 

markers_aggregated_casted$label <- ifelse(markers_aggregated_casted$label2 == TRUE, ifelse(markers_aggregated_casted$melanT < 0, ifelse(markers_aggregated_casted$ovarT < 0, TRUE, FALSE), ifelse(markers_aggregated_casted$ovarT > 0, TRUE,FALSE)),FALSE)

markers_aggregated_casted$significant <- ifelse(markers_aggregated_casted$label == TRUE, "Highlight", ifelse(markers_aggregated_casted$melanT < 0, ifelse(markers_aggregated_casted$ovarT < 0, "Decreased BOTH","Decreased Melanoma"), ifelse(markers_aggregated_casted$ovarT > 0, "Increased BOTH","Increased Melanoma")))

ggplot(markers_aggregated_casted %>% arrange(label), aes(x = ovarT, y = melanT)) +facet_wrap(~cluster) + geom_hline(yintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + geom_point(aes(color = significant)) + theme_bw(base_size = 12)  + xlab("Ovarian CD8 TILs expression KO/WT") + ylab("Melanoma CD8 TILs expression KO/WT") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
  data = subset(markers_aggregated_casted, label == TRUE),
  aes(label = geneSymbol),
  size = 7,
  box.padding = unit(3, "lines"),
  point.padding = unit(0.25, "lines"),
  max.overlaps = Inf) +  scale_color_manual(values = c("#f08080ff","grey","black","#7ac5cdff","grey" )) + ggtitle(paste0("KO vs WT RNA expression in CD8 TILs"))
ggsave(paste0(reclustered.path,"/10_volcano_ALL_melanTvsovarT_topGenes.svg"), width = 10, height = 7)
