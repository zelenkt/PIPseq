#!/usr/bin/env Rscript

#can be run directly from command line by: Rscript 1-Single_Cell_RNAseq_Analysis.R

####---- User Input ----####
## Set local Github repository as working directory
setwd("/Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN")


# Specify which cancer should be processed and adjust the metadata provided

# analysis of mouse recipient data is done sequentially up until line 800 - first for samples that were originally combined with human cells and second for samples that were only mouse recipient cells (below)
organism <- "Mouse"
mouseHuman_sample_origin <- TRUE
input_path <- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/2_recipients_mouse/Input/pipseeker_v3.3"
output_path<- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/2_recipients_mouse/Output"
original_datasets <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen_huTILs_p84","nt_recip_dnmt3_tbetcre_huTILs_p41","wt_recip_bcl11b_cd4ercre_huTILs_p19")
Project_Names <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen","nt_recip_dnmt3_tbetcre","wt_recip_bcl11b_cd4ercre")
Project_Names_all <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen","nt_recip_dnmt3_tbetcre","wt_recip_bcl11b_cd4ercre","ko_recip_bcl11b_cd4ercre","ko_recip_dnmt3_tbetcre","wt_recip_dnmt3_tbetcre")
condition <- c("NT","NT","WT")
disease <- c("BC_ovarT","DT_ovarT","BC_ovarT")
sample <- c("BCN","DTN","BCW")
# conditions for differential analysis - it will be calculated as condition1/condition2 so KO should be first
condition1 <- "NTBCovarT"
#condition2 <- "NTBCovarT"
condition2 <- "WTBCovarT"

# organism <- "Mouse"
# mouseHuman_sample_origin <- FALSE
# input_path <- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/2_recipients_mouse/Input/pipseeker_v3.3"
# output_path<- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/2_recipients_mouse/Output"
# Project_Names <- c("ko_recip_bcl11b_cd4ercre","ko_recip_dnmt3_tbetcre","wt_recip_dnmt3_tbetcre")
# Project_Names_all <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen","nt_recip_dnmt3_tbetcre","wt_recip_bcl11b_cd4ercre","ko_recip_bcl11b_cd4ercre","ko_recip_dnmt3_tbetcre","wt_recip_dnmt3_tbetcre")
# original_datasets <- NULL
# condition <- c("KO","KO","WT")
# disease <- c("BC_ovarT","DT_ovarT","DT_ovarT")
# sample <- c("BCK","DTK","DTW")
# # conditions for differential analysis - it will be calculated as condition1/condition2 so KO should be first
# condition1 <- "NTBCovarT"
# #condition2 <- "NTBCovarT"
# condition2 <- "WTBCovarT"







# Specify what variables should be regressed out
# here they say that we should not regress out nFeatures https://github.com/satijalab/seurat/issues/5667 , vars.to.regress=regress.out
# Also, UMI counts is not necessary set in the vars.to.regress, because it is set in the regularized negative binomial model. https://github.com/satijalab/seurat/issues/1739 ..it is not necessary but it can still be done
regress.out.list <- list(c("percent.mt"),c("nFeature_RNA","percent.mt"),c("nFeature_RNA","percent.mt","S.Score","G2M.Score"),c("nFeature_RNA","percent.mt","S.Score","G2M.Score"),c("percent.mt","S.Score","G2M.Score","condition"))
regress.out <- regress.out.list[[2]]





# Specify whether the integrative analysis will be run and which condition should be used as a reference (if any)
integration <- TRUE

# Select a reference for eference-based integration. Note that we do not recommended it for standard-size datasets.
reference <- NULL
# reference <- "WT" 



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
# resolution <- c(0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,10)
resolution <- c(0.1,0.2,0.3,0.4,0.5,0.7,0.85,1,1.2,1.35,1.5,1.75,2,2.5)



# Setup working resolution for downstream analyses
working_res <- "integrated_snn_res.0.7"


# gene.list <- list(genes=c("Ikzf1","Ets1","Irf1","Batf","Rbpj","Atxn1","Tcf7","Il7r","Ccr7","Cx3cr1","Slamf6","Icos","Cd80","Sell","Cxcr5","Id3","Lamp1","Mki67","Ifng","Il2","Il2ra","Tnf","Gzma","Gzmb","Gzmc","Gzmk","Prf1","Prdm1","Id2","Tox","Pdcd1", "Ctla4","Entpd1","Havcr2","Tigit", "Lag3"))
# 

# gene.list <- list(genes=rev(c("Tcf7", "Id3","Slamf6", "Ccr7", "Sell", "Myb",  "Cxcr5", "Bach2","Bcl6", "Il7r", "Foxo1", "Lef1","Zeb1", "Stat3", "Cd28","Bcl2","Tox","Stat4","Zeb2","Id2","Prdm1", "Batf","Prf1","Gzmb","Gzmc","Gzmf","Nfatc2","Hif1a","Ccl3","Ccl4","Ccl5","Irf4", "Pdcd1","Ctla4", "Tigit", "Lag3", "Havcr2","Entpd1", "Egr2","Nr4a1","Nr4a2")))

gene.list <- list(genes=rev(c("Tcf7", "Id3","Bach2","Bcl6","Slamf6", "Ccr7", "Sell","Il7r","S100a4","Bcl2","Mki67","Gzmb","Gzmc","Gzmf","Prf1","Batf","Tox","Ctla4","Pdcd1","Havcr2","Entpd1","Nr4a1","Nr4a2","Ccl3","Cxcr6")))


annot_wt = c("Cxcr6","Entpd1","Zeb2","Havcr2", "Ctla4","Tox" ,"Prdm1","Id2","Gzmb","Prf1","Mki67","Slamf6","Bach2","Lef1","Tcf7")
stemness = rev(c("Tcf7", "Bach2", "Id3", "Myb", "Foxo1", "Bcl6", "Lef1", "Zeb1", "Stat3", "Slamf6", "Ccr7", "Sell", "Il7r", "Cxcr5", "Cd28","Bcl2"))
exhaustion = c("Gzmf","Gzmc","Gzmb","Prf1","Id2","Ccl5", "Ccl4", "Ccl3", "Entpd1", "Ctla4", "Tigit", "Lag3", "Havcr2", "Pdcd1","Prdm1", "Batf", "Irf4", "Stat4", "Nfatc2", "Nr4a2", "Nr4a1", "Hif1a", "Egr2","Zeb2", "Tox2", "Tox")

fcer1g.signature =c("Cd247","Lat","Ide","Xcl1","Ccr2","Cxcr3","Napsa","Adamts9","Serpinb6a","B2m","Ncr1","Klrb1c","Il2rb","Cd7","Clec4d","Clec4e","Tyrobp","Fcer1g")

fcer1g.signature.long =c("Entpd1","Lag3", "Havcr2","Lat","Vegfc", "Chl1","Galnt18","Enpp6", "Ide","Xcl1","Ccr2","Cxcr3","Napsa","Adamts9","Serpinb6a","B2m","Ncr1","Klrb1c","Il2rb","Cd7","Tyrobp","Fcer1g")








# Setup gene sets to calculate signatures
g1 <- list(Tex =c("Havcr2", "Lag3", "Pdcd1", "Tigit", "Ctla4"),
           Tstem=c("Tcf7", "Lef1", "Sell", "Il7r"),
           Teff=c("Gzmb", "Prf1", "Fasl") 
)
g2 <- list(Tpex =c("Tcf7","Slamf6","Cxcr5","Xcl1","Havcr2-","Entpd1-","Cx3cr1-","Cd101-"),
           Tinex=c("Tcf7-","Slamf6-","Cxcr5-","Xcl1-","Cx3cr1","Cd101-","Cd69-"),
           Ttex=c("Tcf7-","Slamf6-","Cxcr5-","Xcl1-","Havcr2","Entpd1","Cx3cr1-","Cd101","Cd69") 
)
# based on CD8+ T cells in the cancer-immunity cycle from Giles (Wherry) 2023, Immunity; Fig. 1c
#CD8.Tcell.signature = c("Cd2","Cd8a","Cd8b1","Cd3e","Cd3d","Trac"),
#Myeloid_signature = c("Spi1","Fcer1g","Csf1r"),

gene.signatures.list <- list(g1=g1,g2=g2)



####---- Load Packages ----####


##### for comparative analysis of WT and KO samples addittional packages are needed:
# install.packages('BiocManager')
# install.packages('metap')
# install.packages("sctransform")
# BiocManager::install(c("glmGamPoi","celldex","multtest","SingleR")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# devtools::install_github('satijalab/seurat-data')

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix", "fields","SeuratDisk","tibble","SeuratData","SingleR","SingleCellExperiment", "ggplot2", "cowplot", "ggrepel", "UCell", "ggpubr","stringr","SeuratWrappers","reticulate","reshape2","reshape","tidyverse","scCustomize","Nebulosa")
invisible(lapply(packages, library, character.only = TRUE))

set.seed(seed)

# # perhaps not necessary
# library(SignatuR)
# library(DiagrammeR)
# library(GSEABase)



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

  
  # Retrieve save seurat objects
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


if (organism == "MouseHuman") {
  print("Analyzing dataset containing a mixture of mouse and human cells...")
    
  # Process the raw counts for all the samples and perform QC filtering
  for (i in 1:length(Project_Names)){
  
    Project_Name <- Project_Names[i]
    original_dataset <- original_datasets[i]
    
    output_folder <-  paste0(output_path,"/",Project_Name,"/Single_Cell_RNAseq_Output")
    if (!file.exists(output_folder)) {
      dir.create(output_folder,recursive = T)
    }
    
    barnyard_folder <-  paste0(output_path,"/extract_barnyard_plots_and_species_barcodes")
    if (!file.exists(barnyard_folder)) {
      dir.create(barnyard_folder,recursive = T)
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
    
  
    seurat_obj <- CreateSeuratObject(counts = counts, min.cells=3)
    SaveSeuratRds(seurat_obj, file = paste0(output_folder,"/",Project_Name,".Rds"))
    
    # save QC graph before filtering 
    VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    ggsave2(paste0(output_path,"/",Project_Name,"/QC_prefilt_",Project_Name,".pdf"), width = 7, height = 5)
    
    # Filter out cells with less than 100 nFeatures
    seurat_obj <- seurat_obj[, seurat_obj[["nFeature_RNA"]] >= 300]
    
    # save QC graph basic nFeature filtering 
    VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    ggsave2(paste0(output_path,"/",Project_Name,"/QC_post_nFeat_filt_",Project_Name,".pdf"), width = 7, height = 5)
    
    
    raw_counts <- as.data.frame(LayerData(seurat_obj, layer = "counts", assay="RNA"))
    
    # filter out genes expressed in less than 1% of cells
    genes.percent.expression <- rowMeans(raw_counts>0 )*100   
    genes.filter <- names(genes.percent.expression[genes.percent.expression>1])  
    raw_counts.sub <- as.data.frame(raw_counts[genes.filter,])
    
    rm(raw_counts)
    
    raw_counts.sub$species <- ifelse(startsWith(row.names(raw_counts.sub),"hg-") == TRUE, "human",ifelse(startsWith(row.names(raw_counts.sub),"mm-") == TRUE, "mouse","check")) 
    
    # save this and try to analyze it in bash either directly 
    # write.table(raw_counts.sub, file = paste0(output_folder,"/raw_counts_sub_withZeros_",Project_Name,".txt"),sep = "\t", row.names = TRUE, col.names = TRUE)
    
    # replace 0 by NA
    raw_counts.sub[raw_counts.sub == 0] <- NA
    
    # save this and try to analyze it in bash either directly 
    write.table(raw_counts.sub, file = paste0(barnyard_folder,"/raw_counts_sub_withNAs_",Project_Name,".txt"),sep = "\t", row.names = TRUE, col.names = TRUE)
    
  }
  
  stopifnot(organism == "MouseHuman",error=stop("Combined mouse-human samples is used, program extracts mouse and human barcodes that need to be processed in bash, then re-run the script with mouse- or human-specific reference, now quitting...."))
  
  
}



# Process the raw counts for all the samples and perform QC filtering
for (i in 1:length(Project_Names)){
    # i=1
    Project_Name <- Project_Names[i]
    original_dataset <- original_datasets[i]
    
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
    
  print(paste0("Processing dataset containing ",organism," cells"))
  
  # Create a Seurat object and perform QC
  seurat_obj <- CreateSeuratObject(counts = counts, min.cells=3)
  
  
  if (mouseHuman_sample_origin == TRUE) {
    print(paste0("Sample originates from a mixed mouse-human cell dataset: ",original_dataset,", performing isolation of ",organism,"-specific cells..."))
    
    # isolate cell id for mouse and human cells separately
    cast_sum_melt_raw_counts.sub <- read.table(paste0(input_path,"/",Project_Name,"/barnyard_metrics_",original_dataset,".txt"), header = T, col.names=c("barcode","hg_transcripts", "mm_transcripts", "percHuman", "percMouse", "annotation"))
    
    if (organism == "Mouse") { 
      selected_cells <- cast_sum_melt_raw_counts.sub[cast_sum_melt_raw_counts.sub$annotation == "mouse", ]
    } else if (organism == "Human") { 
      selected_cells <- cast_sum_melt_raw_counts.sub[cast_sum_melt_raw_counts.sub$annotation == "human", ]
    } else {
      
      stop("Please check what organism you used")

      # print("Please check what organism you used")
    }
    
    seurat.obj_combined <- seurat_obj
    index<- colnames(seurat.obj_combined) %in% selected_cells$barcode
    seurat_obj <- seurat.obj_combined[,index]
    
    # # this shows gene names
    # rownames(seurat_obj)
    
    write.table(selected_cells$barcode, file = paste0(output_folder,"/",organism,"_cells_barcodes_",Project_Name,".txt"),sep = "\t", row.names = FALSE, col.names = FALSE, quote=F)
  }

  
  seurat_obj <- qc.seurat(seurat_obj, organism, 100, 30000, 300)
  
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
    # seurat_obj <- ScaleData(seurat_obj, vars.to.regress=regress.out)
    seurat_obj <- CellCycleScoring(seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, verbose = FALSE, seed.use = seed)
    
    
    # here they say that we should not regress out nFeatures https://github.com/satijalab/seurat/issues/5667 , vars.to.regress=regress.out
    # Also, UMI counts is not necessary set in the vars.to.regress, because it is set in the regularized negative binomial model. https://github.com/satijalab/seurat/issues/1739 ..it is not necessary but it can still be done
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



#########################################################
# # This was run after I loaded the mixed mouse-human data mapped to mouse genome, filtered with basic qc function but not removing the human data...it was then processed with SCTransform and PCA, UMAP, doublets as regularly - in these subsequent steps I extracted the information on doublets and compared the assumed / calculated doublets with the real mouse/human doublet rate
# # check the distribution and ratio of real and assumed / calculated doublets
# doublets <- as.data.frame(seurat_obj$doublets)
# colnames(doublets) <- "doublets" 
# doublets$barcode <- substring(row.names(doublets),5)
# merged.doublets <- merge(cast_sum_melt_raw_counts.sub,doublets)
# 
# backup <- cast_sum_melt_raw_counts.sub
# 
# cast_sum_melt_raw_counts.sub <- merged.doublets
# 
# multiplet_threshold <- 85
# 
# cast_sum_melt_raw_counts.sub$annotation <- ifelse(cast_sum_melt_raw_counts.sub$percHuman > multiplet_threshold, "human", ifelse(cast_sum_melt_raw_counts.sub$percMouse > multiplet_threshold, "mouse", "multiplets"))
# 
# sum.percent <- paste(paste0("human=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["human"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%"),"\n",paste0("mouse=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["mouse"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%"), "\n",paste0("multiplets=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["multiplets"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%"),"\n",paste0("assumed doublets=",round(100*(table(cast_sum_melt_raw_counts.sub$doublets)["Doublet"] / (table(cast_sum_melt_raw_counts.sub$doublets)[["Doublet"]]+table(cast_sum_melt_raw_counts.sub$doublets)[["Singlet"]])), digits = 1),"%"))
# 
# cast_sum_melt_raw_counts.sub$annotation <- ifelse(cast_sum_melt_raw_counts.sub$doublets == "Doublet" , "assumed doublets",  ifelse(cast_sum_melt_raw_counts.sub$percHuman > multiplet_threshold, "human", ifelse(cast_sum_melt_raw_counts.sub$percMouse > multiplet_threshold, "mouse", "multiplets")))
# 
# # cast_sum_melt_raw_counts.sub$annotation <- factor(cast_sum_melt_raw_counts.sub$annotation, levels=c("assumed doublets","mouse","human","multiplets"))
# 
# 
# cast_sum_melt_raw_counts.sub$annotation <- factor(cast_sum_melt_raw_counts.sub$annotation, levels=c("mouse","human","multiplets","assumed doublets"))
# 
# 
# ggplot(cast_sum_melt_raw_counts.sub %>% arrange(annotation), aes(x = hg_transcripts, y = mm_transcripts))  + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed")  + geom_point(aes(color = annotation)) + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey","red")) + ggtitle(paste0("Barnyard metrics ",Project_Name)) + scale_x_continuous(limits = c(0,10000), breaks = c(0,2000,6000,10000)) +scale_y_continuous(limits = c(0,10000), breaks = seq(0,10000, by = 2000)) + annotate(geom="text", x=8000, y=9000, label=sum.percent ,color="black",size=5)
# ggsave2(paste0(output_path,"/",Project_Name,"/barnyard_doublets_",Project_Name,".png"), width = 9, height = 5.5)
# ggsave2(paste0(output_path,"/",Project_Name,"/barnyard_doublets_",Project_Name,".svg"), width = 9, height = 5.5)
# 
# ggplot(cast_sum_melt_raw_counts.sub %>% arrange(annotation), aes(x = log10(hg_transcripts), y =log10(mm_transcripts))) + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed") + geom_point(aes(color = annotation)) + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey","red")) + ggtitle(paste0("Barnyard metrics ",Project_Name))+ scale_x_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) +scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) + annotate(geom="text", x=1, y=1, label=sum.percent ,color="black",size=5)
# 
# ggsave2(paste0(output_path,"/",Project_Name,"/barnyard_doublets_log10_",Project_Name,".png"), width = 9, height = 5.5)
# ggsave2(paste0(output_path,"/",Project_Name,"/barnyard_doublets_log10_",Project_Name,".svg"), width = 9, height = 5.5)
############################################################################







############### from this point onwards, all mouse data can be analyzed together
if (organism == "Mouse") {
Project_Names <- Project_Names_all
}  else if (organism == "Human") {
  Project_Names <- Project_Names
} else {
  stop("Please check what organism you used")
  # print("Please check what organism you used")
}


for (i in 1:length(Project_Names)){
  # i=1
  Project_Name <- Project_Names[i]
  original_dataset <- original_datasets[i]
  working_res <- "SCT_snn_res.0.5"
  working_assay <- "RNA"
  
  feat_lymphocyte <- c("Thy1","Cd3e","Trac","Cd8a","Cd4","Cd19")
  feat_macrophage <- c("Itgam","Adgre1","Itgax","Ly6c1","Ly6g","Cxcr2","Cd33","Cd9","Cxcr4")
  feat_nk <- c("Klrg1","Klrb1","Klrk1","Klrd1", "Cd27", "Ncr1", "Ly6a", "Ncam1","Fcgr3")
  
  genes2 <- c("Thy1","Cd3e","Cd3d","Trac","Cd4","Cd8a","Cd8b1","Cd160","Klrk1","Klrd1","Cd19","Ms4a1","Cd22","Itgam","Adgre1","Itgax","Cxcr2","Ly6c1","Ly6g","Cd33", "Ccr5","Cr1l","Cxcr4","Ccr1","Cd14","Fcgr3","Cd9","Ptprc","Cpe","Mme","Cma1","Tlr7","Il10ra","Ptger2","Ifnar1","Itga4","Cd36","Fcgr2b")
              
           
  g6 <- list(	CD4=c("Cd4","Cd3e","Cd3d","Trac","Cd8a-","Cd8b1-"),
              CD8=c("Cd8a","Cd8b1","Cd3e","Cd3d","Trac","Cd4-"),
              Bcell=c("Cd19","Ms4a1"),
              NK=c("Klrg1","Klrb1","Klrk1","Klrd1", "Cd27", "Ncr1", "Ly6a", "Ncam1","Cd3e-","Cd3d-"),
              monoMacro=c("Itgam"),
              actMacro=c("Itgam","Adgre1"),
              dendritic=c("Itgax"),
              neutro=c("Ly6c1","Ly6g","Cxcr2","Cd33","Cd9","Cxcr4")
  )
  gene.signatures.list <- list(g6=g6)
  
  
  
  output_folder <-  paste0(output_path,"/",Project_Name,"/Single_Cell_RNAseq_Output")
  seurat_obj <- LoadSeuratRds(paste0(output_folder,"/",Project_Name,".Rds"))
  

  
  DefaultAssay(seurat_obj) <-  working_assay
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  print(table(seurat_obj$sample))
  
  print(table(seurat_obj$doublets))
  

  DimPlot(object = seurat_obj, group.by = "doublets", label = T, seed = seed, repel = TRUE)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_doublets_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_doublets_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  
  DimPlot(object = seurat_obj, group.by = "ImmGen.main", label = T, seed = seed, repel = TRUE)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_ImmGen_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_ImmGen_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  # DimPlot(object = seurat_obj, group.by = "ImmGen.fine", label = T)
  
  DimPlot(object = seurat_obj, group.by = "mouseRNAseq.main", label = T, seed = seed, repel = TRUE)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_mouseMain_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_mouseMain_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  DimPlot(object = seurat_obj, group.by = working_res, label = T, seed = seed, repel = TRUE)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_clusters_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_clusters_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  DimPlot(object = seurat_obj, group.by = "Phase", label = T, seed = seed, repel = TRUE)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_phase_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/umap_phase_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  # DotPlot(seurat_obj, features = c("Thy1","Cd3e","Cd3d","Trac","Cd4","Cd8a","Cd8b1","Cd160","Klrk1","Klrd1","Cd19","Ms4a1","Ccr3", "Cd33", "Il5ra","Ccr5", "Cxcr2","Enpp3","Fcer1a","Cr1l", "Il2ra", "Ptgdr2","Ptprc","Ctsg","Fut4","Rnase2a","Tlr7","Cd14","Elane","Fcgr3","Mme","Mpo","Cma1","Cpe","Cxcr4","Il10ra","Ptger2","Ccr1","Ifnar1","Itga4","Cd9","Cd22","Cd36","Cd40lg","Fcgr2b"),group.by = working_res, assay=working_assay, scale.min = 0,scale.max =50 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  # ggsave2(paste0(output_path,"/",Project_Name,"/dotPlot_Genes1_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 7, height = 9)
  # ggsave2(paste0(output_path,"/",Project_Name,"/dotPlot_Genes1_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 7, height = 9)
  
  DotPlot(seurat_obj, features = genes2 ,group.by = working_res, assay=working_assay, scale.min = 0,scale.max =100 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000") 
  ggsave2(paste0(output_path,"/",Project_Name,"/dotPlot_Genes2_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 7, height = 7)
  ggsave2(paste0(output_path,"/",Project_Name,"/dotPlot_Genes2_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 7, height = 7)
  
  # Shamima:
  # Macrophages + Monocytes: Cd11b=itgam...from this activated macrophages have F4/80=Adgre1 
  # Dendritic cell Cd11c=itgax
  # Neutro Ly6c, Ly6g, Gr1
  # 
  # 
  # CD16 = Fcgr3 also known as FcγRIII, is a cluster of differentiation molecule found on the surface of natural killer cells, neutrophils, monocytes, macrophages, and certain T cells
  # 
  # 
  # granulo
  # "Ccr3", "Cd33", "Il5ra"
  # 
  # nm
  # "Ccr5", "Cxcr2"
  # 
  # mb
  # "Enpp3","Fcer1a"
  # 
  # eb
  # "Cr1", "Il2ra", "Ptgdr2","Ptprc" # Cr1 replaced by Cr1l
  # 
  # ne
  # "Cd65","Ctsg","Fut4","Rnase3","Tlr7" # Rnase3 replaced by Rnase2a
  # 
  # neutro
  # "Cd14","Elane","Fcgr3","Mme","Mpo" # "Fcgr3a"= Fcgr3
  # 
  # mast
  # "Cma1","Cpe","Cxcr4","Il10ra","Ptger2"
  # 
  # eosin
  # "Ccr1","Fcar","Ifnar1","Itga4","Siglec8"
  # 
  # baso
  # "Cd9","Cd22","Cd36","Cd40lg","Fcgr2b"
  
  
  # table(seurat_obj$SCT_snn_res.0.1) # total 5875, with 112 T cells
  
  tryCatch({
    
    signatures <- AddModuleScore_UCell(seurat_obj, features = gene.signatures.list[[1]])
    signature.names <- paste0(names(gene.signatures.list[[1]]), "_UCell")
    
    # generate violin plots
    # VlnPlot(signatures, features = signature.names, group.by = working_res)
    Stacked_VlnPlot(seurat_object = signatures, features = signature.names, x_lab_rotate = TRUE, group.by = working_res)
    ggsave2(paste0(output_path,"/",Project_Name,"/vln_signatures_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
    ggsave2(paste0(output_path,"/",Project_Name,"/vln_signatures_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  
  tryCatch({
    

  FeaturePlot_scCustom(seurat_object = seurat_obj, features = feat_lymphocyte,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_Lymphocyte_",working_assay,"_",Project_Name,".svg"), width = 12, height = 10)
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_Lymphocyte_",working_assay,"_",Project_Name,".png"), width = 12, height = 10)
  
  Plot_Density_Custom(seurat_object = seurat_obj, features = feat_lymphocyte,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_Lymphocyte_",working_assay,"_",Project_Name,".svg"), width = 12, height = 6.7)
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_Lymphocyte_",working_assay,"_",Project_Name,".png"), width = 12, height = 6.7)
  
  FeaturePlot_scCustom(seurat_object = seurat_obj, features = feat_macrophage,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_Macrophage_",working_assay,"_",Project_Name,".svg"), width = 12, height = 10)
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_Macrophage_",working_assay,"_",Project_Name,".png"), width = 12, height = 10)
  
  Plot_Density_Custom(seurat_object = seurat_obj, features = feat_macrophage,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_Macrophage_",working_assay,"_",Project_Name,".svg"), width = 12, height = 10)
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_Macrophage_",working_assay,"_",Project_Name,".png"), width = 12, height = 10)
  
  FeaturePlot_scCustom(seurat_object = seurat_obj, features = feat_nk,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_NK_",working_assay,"_",Project_Name,".svg"), width = 12, height = 10)
  ggsave2(paste0(output_path,"/",Project_Name,"/featPlot_NK_",working_assay,"_",Project_Name,".png"), width = 12, height = 10)
  
  Plot_Density_Custom(seurat_object = seurat_obj, features = feat_nk,reduction = "umap")
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_NK_",working_assay,"_",Project_Name,".svg"), width = 12, height = 10)
  ggsave2(paste0(output_path,"/",Project_Name,"/densPlot_NK_",working_assay,"_",Project_Name,".png"), width = 12, height = 10)
  # Plot_Density_Custom(seurat_object = seurat_obj, features = c("Ncam1"),reduction = "umap")
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  tryCatch({
    
  # Create Plots
  Stacked_VlnPlot(seurat_object = seurat_obj, features = feat_lymphocyte, x_lab_rotate = TRUE, group.by = working_res)
  # , colors_use = c("#33A02C","#F781BF","cyan3","dodgerblue1","#E6AB02","purple"))
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_Lymphocyte_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 4.3)
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_Lymphocyte_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 4.3)
  
  Stacked_VlnPlot(seurat_object = seurat_obj, features = feat_macrophage, x_lab_rotate = TRUE, group.by = working_res)
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_Macrophage_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_Macrophage_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
  
  Stacked_VlnPlot(seurat_object = seurat_obj, features = feat_nk, x_lab_rotate = TRUE, group.by = working_res)
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_NK_",working_assay,"_",working_res,"_",Project_Name,".svg"), width = 8, height = 6)
  ggsave2(paste0(output_path,"/",Project_Name,"/vln_NK_",working_assay,"_",working_res,"_",Project_Name,".png"), width = 8, height = 6)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}








# c("Ptprc","Cd11b","Cd11c","Cd3e","Cd3d","Trac","Cd4","Cd8a","Cd8b1","Cd19","Ms4a1","Cd160","Klrd1","Klrk1")
# 
# g6 <- list(	immune=c("Ptprc"),
#             CD4=c("Cd4","Cd3e","Cd3d","Trac","Cd8a-","Cd8b1-"),
#             CD8=c("Cd8a","Cd8b1","Cd3e","Cd3d","Trac","Cd4-"),
#             Bcell=c("Cd19","Ms4a1"),
#             NK=c("Klrb1","Cd160","Gnly","Ncam1"),
#             myeloid=c("Cd11c","Cd11b"),
#             monocytes=c("Ly6c2","Ly6c1","Mgst1","Ccl9","Fcgr1","Cxcl10","Nr4a1","Fn1","Klf13","Itgal","Cd300a","Clec4a1","Lyz2","Chil3","S100a4","Vcan","Cd300e","Ace","Fcgr4","Il1b","S100a8","S100a9"),
#             cDC2=c("Cd74","Ciita","Csf1r","Cst3","H2-Aa","H2-Ab1","H2-DMa","H2-Eb1","H2-Oa","Irf4","Itgax","Siglecg","Sirpa"),
#             cDC1=c("B2m","Cd226","Clec9a","Dpp4","Gcsam","Ifi205","Il12b","Irf8","Itgae","Itgax","Nlrc5","Rab43","Swap70","Tap1","Tap2","Tapbp","Tlr3","Xcr1","Zdhhc2"),
#             macrophages=c("Cd68","Cd163","Mrc1","C1qc","C1qb","Apoe","Mafb","Ccr5","Ms4a6b","Maf","Cx3cr1","Trem2","Mgl2","Clec4b1","Lrp1","Spp1","Vegfa","Ssp1","Mmp12"),
#             stemness=c("Il7r","Tcf7", "Lef1", "Bach2", "Bcl6","Slamf6","Ccr7", "Dpp4"), 
#             effector=c("Tbx21","Id2","Prdm1","Fasl", "Gzmb", "Prf1"),
#             exhaustion =c("Havcr2", "Lag3", "Pdcd1", "Tigit", "Ctla4", "Nr4a1", "Nr4a2","Nr4a3","Egr2")
# )
# 
# 
# signatures <- AddModuleScore_UCell(seurat_obj, features = gene.signatures.list[[7]])
# signature.names <- paste0(names(gene.signatures.list[[7]]), "_UCell")
# 
# # generate violin plots
# VlnPlot(signatures, features = signature.names, group.by = working_res)
#         
# ggsave2(paste0(signatures.path,"/signatures_",names(gene.signatures.list[i]),"_vln_split_disCond_",working_res,".svg"), width = 20, height = 5)


# 
# # load the table created with HPC and test different plotting options
# cast_sum_melt_raw_counts.sub <- read.table(paste0(output_folder,"/barnyard_metrics_",Project_Name,".txt"), header = TRUE)
# 
# table(cast_sum_melt_raw_counts.sub$annotation)[1]
# 
# sum.percent <- paste(paste0("human=",round(100*(table(cast_sum_melt_raw_counts.sub$label)["human"] / (table(cast_sum_melt_raw_counts.sub$label)[["human"]]+table(cast_sum_melt_raw_counts.sub$label)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$label)[["multiplets"]])), digits = 1),"%"),"\n",paste0("mouse=",round(100*(table(cast_sum_melt_raw_counts.sub$label)["mouse"] / (table(cast_sum_melt_raw_counts.sub$label)[["human"]]+table(cast_sum_melt_raw_counts.sub$label)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$label)[["multiplets"]])), digits = 1),"%"), "\n",paste0("multiplets=",round(100*(table(cast_sum_melt_raw_counts.sub$label)["multiplets"] / (table(cast_sum_melt_raw_counts.sub$label)[["human"]]+table(cast_sum_melt_raw_counts.sub$label)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$label)[["multiplets"]])), digits = 1),"%") )
# 
# 
# ggplot(cast_sum_melt_raw_counts.sub, aes(x = human, y = mouse))  + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed")  + geom_point(aes(color = annotation)) + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey")) + ggtitle(paste0("Barnyard metrics ",Project_Name)) + annotate(geom="text", x=8000, y=9000, label=sum.percent ,color="black",size=5)
# 
# ggsave(paste0(figures.path,"/volcano_only_cl",i,"_",names(gene.list[g]),"_allMarkers_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
# 
# 
# ggplot(cast_sum_melt_raw_counts.sub, aes(x = log10(human), y =log10(mouse))) + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed")  + geom_point(aes(color = annotation))  + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey")) + ggtitle(paste0("Barnyard metrics ",Project_Name))+ scale_x_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) +scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) + annotate(geom="text", x=1, y=1, label=sum.percent ,color="black",size=5)
# 
# write.table(cast_sum_melt_raw_counts.sub, file = paste0(output_folder,"/barnyard_metrics_",Project_Name,".txt"),sep = "\t", row.names = FALSE, col.names = TRUE)

  
  