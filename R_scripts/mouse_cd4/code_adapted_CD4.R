#!/usr/bin/env Rscript

#can be run directly from command line by: Rscript 1-Single_Cell_RNAseq_Analysis.R

####---- User Input ----####
## Set local Github repository as working directory
setwd("/Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN")


## Specify organism
#organism= "Human"
organism <- "Mouse"

# Specify the input and output folder paths
input_path <- "2-Single_Cell_RNAseq_Pipeline/Input"
output_path<- "2-Single_Cell_RNAseq_Pipeline/Output"

# Specify what variables should be regressed out
regress.out.list <- list(c("nFeature_RNA","percent.mt"),c("nFeature_RNA","percent.mt","S.Score","G2M.Score"))
regress.out <- regress.out.list[[1]]


# Specify which cancer should be processed and adjust the metadata provided

Project_Names <- c("ko_ovarian_mouse_tumor_cd4_ot2_rep1","wt_ovarian_mouse_tumor_cd4_ot2_rep1")
condition <- c("KO","WT")
disease <- c("ovarT4","ovarT4")
sample <- c("OT4K1","OT4W1")
# conditions for differential analysis - it will be calculated as condition1/condition2 so KO should be first
condition1 <- "KOovarT4"
condition2 <- "WTovarT4"

# Project_Names <- c("ko_ovarian_mouse_tumor_cd8_ot1_rep1","ko_ovarian_mouse_tumor_cd8_ot1_rep2","wt_ovarian_mouse_tumor_cd8_ot1_rep1","wt_ovarian_mouse_tumor_cd8_ot1_rep2")
# condition <- c("KO","KO","WT","WT")
# disease <- c("ovarT","ovarT","ovarT","ovarT")
# sample <- c("OTK1","OTK2","OTW1","OTW2")
# # conditions for differential analysis - it will be calculated as condition1/condition2 so KO should be first
# condition1 <- "KOovarT"
# condition2 <- "WTovarT"

# Project_Names <- c("ko_melanoma_mouse_tumor_cd8_pmel_rep1","ko_melanoma_mouse_tumor_cd8_pmel_rep2","wt_melanoma_mouse_tumor_cd8_pmel_rep1","wt_melanoma_mouse_tumor_cd8_pmel_rep2")
# condition <- c("KO","KO","WT","WT")
# disease <- c("melanT","melanT","melanT","melanT")
# sample <- c("MTK1","MTK2","MTW1","MTW2")
# condition1 <- "KOmelanT"
# condition2 <- "WTmelanT"

# Setup extension for labeling
extension <- "ovarT4_integrated_noRefUsed"


# Specify whether the integrative analysis will be run and which condition should be used as a reference (if any)
integration <- TRUE
# reference <- "WT"
reference <- NULL

### Specify additional parameters
# number of variable features for integration
varFeatures <- 3000

# lower and upper filtering thresholds for nFeatures and upper threshold for mitochondrial content
MINnFeature <- 400
MAXnFeature <- 20000
mtPerc <- 20


# Number of randomly sampled cells
num_cells <- 1000


# Define the file names for count file and/or for raw matrix data from CellRanger / PIPseeker
Count_file_name <- "Count_File.txt.gz"
Matrix_name <- "matrix.mtx.gz"
Barcodes_name <- "barcodes.tsv.gz"
Features_name <- "features.tsv.gz"

# Define the file name for clinical/meta data (OPTIONAL)
meta_name <- "Meta_File_Example_GSE116256_AML921A-D0.txt.gz"


seed <- 42

# Setup range of resultions for cluster generation
resolution <- c(0.1,0.2,0.3,0.4,0.5,0.7,1,1.5,2)

# Setup working resolution for downstream analyses
working_res <- "integrated_snn_res.1"



# Setup gene sets to calculate signatures
g0 <- list(Th1a=c("Nkg7","Gnly","Gzmh","Ccl5","Ccl4","Gzmb","Gzma","Prf1","Ifng","Gzmm","Klrg1","Zeb2","Cx3cr1","Tnf","Il2rg","Ly6e","Tbx21","Gbp5","Il32","Runx3","Ccl3"),
           Th1b=c("Cxcr4","Pecam1","Pdcd1","Il2rb","Eomes","Klrg1","Cd27","Tnf","Gzmm","Ifng","Nkg7","Gzma","Ccl4","Ccl5","Gzmk"),
           Th1c=c("Gbp5","Samd3","Klrg1","Hopx","Il32","Ifng","Gzmm","Cxcr3","Ccl5","Gzma","Gzmk"),
           Tr1=c("Id2","Irf4","Ly6e","Icos","Gzmm","Cxcr6","Batf","Prdm1","Cxcr3","Tigit","Il21","Cd27","Gzmk","Il2rb","Gzma","Il32","Cd38","Il10","Cd59","Ctla4","Lag3"),
           Treg=c("Icos","Gata3","Ccr4","Il10ra","Prdm1","Il2ra","Batf","Tigit","Cd27","Il32","Sell","Ikzf2","Ctla4","Foxp3"),
           Tcm=c("Icos","Cd27","Areg","Il6st","Sell","Ccd7","Lef1","Tcf7"),
           Th2=c("Il13","Lmo4","Socs3","Ccr10","Tnfsf10","Lgals1","Gata3"),
           Tfhcirc=c("Il7r","Ccr6","Ltb","Aqp3","Tnfsf13b","Cmtm6","Lgals3","Tob1","Zfp36","Bcl6")
)

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

g3 <- list(Tpex =c("Tcf7","Bcl6","Foxo1","Id3","Fli1","Myb","Eomes","Stat3","Bach2","Tox","Runx1"),
           Tpexdif=c("Id2","Batf","Irf4","Tbx21","Runx3","Zeb2","Irf7","Prdm1","Erg2","Hif1a"))
# from Stem-like exhausted and memory CD8+ T cells in cancer Gebhardt, 2023, Nat Rev Canc

g4 <- list(pA.eff=c("Ifng","Zeb2","Myc","Il2ra","Gzma","Gzmk"),
           pB.exh=c("Mt1","Mt2","Entpd1","Klra5","Klre1","Cd38"),
           pC.stem=c("Nt5e","Tcf7","Lef1","Id3","Sell","Slamf6"),
           pD.prol=c("Mki67","Birc5","Cdk1"),
           M2 =c("Zeb2","Id2"),
           M3=c("Nr4a3","Nr4a2","Ikzf1","Nfat5","Runx3"),
           M5 =c("Bach2","Bcl6","Tcf12","Jund"),
           M6 =c("Foxo1","Tox2","Fhl2","Glis1"),
           M7 =c("Tcf7","Myb","Ets1","Plagl1"),
           M8 =c("Tox","Nfatc2","Rbpj","Ikzf2","Klf13")
)
# M - modules from Hongbo paper 2023 Nature Fig1e

g5 <- list(Tpex1=c("Tcf7","Slamf6","Myb","Sell","Bach2","Havcr2-","Entpd1-","Cx3cr1-"),
           Tpex2=c("Tcf7","Slamf6","Myb","Sell","Bach2","Havcr2-","Entpd1-","Cx3cr1-","Mki67"),
           Tex1=c("Tcf7-","Slamf6-","Sell-","Havcr2","Pdcd1","Entpd1","Cx3cr1","Mki67","Ifng"),
           Tex2=c("Tcf7-","Slamf6-","Sell-","Havcr2","Pdcd1","Entpd1","Cd38","Cd244a","Cxcr6","Mki67-")
)

g6 <- list(Stem=c("Bcl2","Cd28","Cxcr5","Il7r","Sell","Ccr7","Slamf6","Stat3","Zeb1","Lef1","Bcl6","Foxo1","Myb","Id3","Bach2","Tcf7"),
           Exhaust=c("Gzmf","Gzmc","Gzmb","Ccl5", "Ccl4", "Ccl3", "Entpd1", "Ctla4", "Tigit", "Lag3", "Havcr2", "Pdcd1",
                     "Prdm1", "Batf", "Irf4", "Stat4", "Nfatc2", "Nr4a2", "Nr4a1", "Hif1a", "Egr2", "Id2",
                     "Zeb2", "Tox2", "Tox")
)
# M - modules from Hongbo paper 2023 Nature Fig1j


gene.signatures.list <- list(g0=g0,g1=g1,g2=g2,g3=g3,g4=g4,g5=g5,g6=g6)



####---- Load Packages ----####


##### for comparative analysis of WT and KO samples addittional packages are needed:
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
#install.packages("sctransform")
#BiocManager::install("glmGamPoi")

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix", "fields","SeuratDisk","tibble","SeuratData","SingleR","SingleCellExperiment", "ggplot2", "cowplot", "ggrepel", "UCell", "SAVER", "ggpubr")
invisible(lapply(packages, library, character.only = TRUE))

# perhaps not necessary
library(SignatuR)
library(DiagrammeR)
library(GSEABase)



#### Define functions


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
  
  # Write seurat H5 file 
  SaveH5Seurat(seurat_obj, filename = file.path(output_folder, paste0(Project_Name,"_h5friendly")),
               overwrite = TRUE, verbose = TRUE)
  
  # write out count data
  for (j in names(seurat_obj@assays)) {
    
    # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
    # scaled_counts <- as.data.frame(seurat_obj@assays[[j]]@scale.data)
    raw_counts <- as.data.frame(seurat_obj@assays[[j]]@counts)
    normalized_data <- as.data.frame(seurat_obj@assays[[j]]@data)
    
    # Add row names as a column for Alyssa Umap app
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
      seurat_obj <- LoadH5Seurat(paste0(output_path.list[i]))
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
      
    
    if (!file.exists(paste0(output_path,"/",extension))) {
      dir.create(paste0(output_path,"/",extension),recursive = T)
    }
    
    SaveH5Seurat(integrated.seurat.object, filename = file.path(output_path,extension, paste0("/",extension,"_seurat_obj")), overwrite = TRUE, verbose = TRUE)
    
    raw_counts <- as.data.frame(integrated.seurat.object@assays$RNA@counts)
    normalized_data <- as.data.frame(integrated.seurat.object@assays$RNA@data)
    metadata <- as.data.frame(integrated.seurat.object@meta.data)
    
    # Add row names as a column for the Umap app
    raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
    normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
    metadata <- tibble::rownames_to_column(metadata, var = "Cell")
    
    # Write raw_counts, normalized_data, and metadata to separate files
    write.table(raw_counts, file = paste0(output_path,"/",extension,"/",extension,"_raw_counts.txt"),sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(normalized_data, file = paste0(output_path,"/",extension,"/",extension,"_normalized_counts.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(metadata, file = paste0(output_path,"/",extension,"/",extension,"_metafile.txt"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    
    # save dimplots as pdf
    DimPlot(integrated.seurat.object, label = T, repel = T, group.by = "Phase") + ggtitle("Unsupervised clustering grouped by Phase")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_phase.svg"), width = 7, height = 5)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, group.by = "disease") + ggtitle("Unsupervised clustering grouped by Disease")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_disease.svg"), width = 7, height = 5)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, group.by = "condition") + ggtitle("Unsupervised clustering grouped by Condition")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_cond.svg"), width = 7, height = 5)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, group.by = "sample") + ggtitle("Unsupervised clustering grouped by Sample")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_sample.svg"), width = 7, height = 5)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, group.by = "conditionDisease") + ggtitle("Unsupervised clustering grouped by Condition-Disease")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_condDis.svg"), width = 7, height = 5)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "sample",group.by = "integrated_snn_res.0.5") + ggtitle("Unsupervised clustering split by Sample grouped by integrated_snn_res.0.5")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_sample_split_res005.svg"), width = 14, height = 4)
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "diseaseCondition",group.by = "integrated_snn_res.0.5") + ggtitle("Unsupervised clustering split by Disease-Condition grouped by integrated_snn_res.0.5")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_condDis_split_res005.svg"), width = 10, height = 4)
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "sample",group.by = "integrated_snn_res.0.3") + ggtitle("Unsupervised clustering split by Sample grouped by integrated_snn_res.0.3")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_sample_split_res003.svg"), width = 14, height = 4)
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "diseaseCondition",group.by = "integrated_snn_res.0.3") + ggtitle("Unsupervised clustering split by Disease-Condition grouped by integrated_snn_res.0.3")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_condDis_split_res003.svg"), width = 10, height = 4)
    
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "condition",group.by = "integrated_snn_res.0.5") + ggtitle("Unsupervised clustering split by Condition grouped by integrated_snn_res.0.5")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_cond_split_res005.svg"), width = 10, height = 4)
    
    DimPlot(integrated.seurat.object, label = T, repel = T, split.by = "condition",group.by = "integrated_snn_res.0.3") + ggtitle("Unsupervised clustering split by Condition grouped by integrated_snn_res.0.3")
    ggsave2(paste0(output_path,"/",extension,"/",extension,"_dimplot_cond_split_res003.svg"), width = 10, height = 4)
    
    return(integrated.seurat.object)
    
  } else {
    print("Integration analysis will be skipped....")
  }
}


####---- Run Script ----####
set.seed(seed)


# Process the raw counts for all the samples and perform QC filtering
for (i in 1:length(Project_Names)){
  i=2
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
    
  } else {
    print(paste0(Project_Name," raw counts were already processed, skipping..."))
  }
}

# integrate and save the data
output_path.list <- c()
for (i in 1:length(Project_Names)){
  path <- paste0(output_path,"/",Project_Names[i],"/Single_Cell_RNAseq_Output/",Project_Names[i],"_h5friendly.h5seurat")
  output_path.list <- append(output_path.list,path)
}

integrated.seurat.object <- integrate.save.data(output_path, output_path.list, extension, Project_Names, varFeatures, integration, reference, condition, resolution, seed)


######### generate more detailed graphs for expression and signatures
extension <- "ovarT4_integrated_noRefUsed_100nFeat_filt_WTKO"
integrated.path <- paste0(output_path,"/",extension)
scData <- LoadH5Seurat(paste0(integrated.path,"/",extension,"_seurat_obj.h5seurat"))

# scData <- LoadH5Seurat(paste0("2-Single_Cell_RNAseq_Pipeline/Output/wt_ovarian_mouse_tumor_cd4_ot2_rep1/Single_Cell_RNAseq_Output/wt_ovarian_mouse_tumor_cd4_ot2_rep1_h5friendly.h5seurat"))

# integrated.path = "2-Single_Cell_RNAseq_Pipeline/Output/wt_ovarian_mouse_tumor_cd4_ot2_rep1"
# working_res = "SCT_snn_res.1"

scData@meta.data$diseaseCondition <- paste0(scData@meta.data$disease,scData@meta.data$condition)


# Randomly sample cells
random_cells <- sample(1:ncol(scData), size = num_cells, replace = FALSE)
subset_seurat_obj <- scData[,random_cells]

# write out count data
for (j in names(subset_seurat_obj@assays)) {
  
  # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
  raw_counts <- as.data.frame(subset_seurat_obj@assays[[j]]@counts)
  normalized_data <- as.data.frame(subset_seurat_obj@assays[[j]]@data)
  metadata <- as.data.frame(subset_seurat_obj@meta.data)
  
  # Add row names as a column for Alyssa Umap app
  raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
  normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
  metadata <- tibble::rownames_to_column(metadata, var = "Cell")
  
  # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
  write.table(raw_counts, file = paste0(integrated.path, "/", extension,"_",num_cells,"_",j,"_raw_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(normalized_data, file = paste0(integrated.path, "/", extension,"_",num_cells,"_",j,"_normalized_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(metadata, file = paste0(integrated.path,"/",extension,"_",num_cells,"_",j,"_metafile.txt"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}


figures.path <- paste0(integrated.path,"/figures")
if (!file.exists(figures.path)) {
  dir.create(figures.path,recursive = T)
}

  working_assay <- "RNA"
  DefaultAssay(scData) <- "RNA"
  Idents(scData) <- working_res
  markers <- FindAllMarkers(scData, assay = working_assay, random.seed = seed) # only.pos = TRUE
  
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene))  %>% #filters out all ribosomal and mitochondrial genes
    slice_head(n = 10) %>%
    ungroup() -> top10
  top10 = top10[!duplicated(top10$gene),]
  
  
  write.table(markers, file = paste0(figures.path,"/",extension,"_allMarkers_combined_cond_RNA_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  DotPlot(scData, features = top10$gene, group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size')+coord_flip()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_allMarkers_no_split_noRiboMt_ClustMarkers_assay",working_assay,"_",working_res,".svg"), width = 7, height = 20)

  DotPlot(scData, features = top10$gene,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "condition") + coord_flip() +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_allMarkers_split_cond_noRiboMt_assay",working_assay,"_",working_res,".svg"), width = 10, height = 20)
  
  DotPlot(scData, features = top10$gene,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "sample") + coord_flip() +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_allMarkers_split_sample_noRiboMt_ClustMarkers_assay",working_assay,"_",working_res,".svg"), width = 24, height = 20)
  
  DotPlot(scData, features = top10$gene,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip() +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_allMarkers_split_condDis_noRiboMt_ClustMarkers_assay",working_assay,"_",working_res,".svg"), width = 14, height = 20)
  
  DotPlot(scData, features = gene.list,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size') + coord_flip()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_selectMarkers_no_split_assay",working_assay,"_",working_res,".svg"), width = 6, height = 9)
  
  DotPlot(scData, features = gene.list,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "condition") + coord_flip()  +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_selectMarkers_split_cond_assay",working_assay,"_",working_res,".svg"), width = 10, height = 9)
  
  DotPlot(scData, features = gene.list,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "sample") + coord_flip()  +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_selectMarkers_split_sample_assay",working_assay,"_",working_res,".svg"), width = 24, height = 9)

  DotPlot(scData, features = gene.list,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "conditionDisease") + coord_flip()  +RotatedAxis()
  ggsave2(paste0(figures.path,"/",extension,"_bubble_selectMarkers_split_condDis_assay",working_assay,"_",working_res,".svg"), width = 14, height = 9)
  
  VlnPlot(scData, features = tcell_features, split.by = "conditionDisease",group.by = working_res, pt.size = 0, combine = T, ncol = 2)
  ggsave2(paste0(figures.path,"/",extension,"_violin_tcell_splitCondDis_assay",working_assay,"_",working_res,".svg"), width = 16, height = 8)

  VlnPlot(scData, features = stem_features, split.by = "conditionDisease",group.by = working_res, pt.size = 0, combine = T, ncol = 3)
  ggsave2(paste0(figures.path,"/",extension,"_violin_stem_splitCondDis_assay",working_assay,"_",working_res,".svg"), width = 25, height = 16)
  
  VlnPlot(scData, features = effector_features, split.by = "conditionDisease",group.by = working_res, pt.size = 0, combine = T, ncol = 2)
  ggsave2(paste0(figures.path,"/",extension,"_violin_effector_splitCondDis_assay",working_assay,"_",working_res,".svg"), width = 16, height = 8)
  VlnPlot(scData, features = exhaustion_features, split.by = "conditionDisease",group.by = working_res, pt.size = 0, combine = T, ncol = 2)
  ggsave2(paste0(figures.path,"/",extension,"_violin_exhaust_splitCondDis_assay",working_assay,"_",working_res,".svg"), width = 16, height = 8)
  
    # fix the umap coordinates for feature and signature plots https://github.com/satijalab/seurat/issues/2507 https://satijalab.org/seurat/articles/essential_commands.html
  DefaultAssay(scData) <- "integrated"
  new.embeddings = Embeddings(scData, reduction = "umap")
  new_reduction <- CreateDimReducObject(embeddings = new.embeddings, key = "umap")
  scData[["umap"]] <- new_reduction

  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = working_res) + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_clusters_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_cond_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "Phase") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_phase_",working_res,".svg"), width = 7, height = 5)

  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "sample") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_sample_",working_res,".svg"), width = 7, height = 5)

  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = working_res, split.by = "diseaseCondition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-diseaseCondition_",working_res,".svg"), width = 12, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, split.by = "disease",group.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-disease_groupCond_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "conditionDisease") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_conditionDisease_",working_res,".svg"), width = 7, height = 5)

  DimPlot(scData, reduction= "umap", label = T, repel = T,group.by = working_res, split.by = "condition") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-condition_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T,group.by = working_res, split.by = "disease") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_split-disease_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "ImmGen.main") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_ImmGen.main_",working_res,".svg"), width = 7, height = 5)
  
  DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "mouseRNAseq.main") + ggtitle("Unsupervised clustering")
  ggsave2(paste0(figures.path,"/",extension,"_dimplot_mouseRNAseq.main_",working_res,".svg"), width = 7, height = 5)
  
  

  DefaultAssay(scData) <- "RNA"
  # FeaturePlot(scData, features = tcell_features, max.cutoff = "q80", cols=c("grey","red"),ncol=2, pt.size = 0.2 )
  # ggsave2(paste0(figures.path,"/",extension,"_featPlot_tcell_assay",working_assay,".svg"), width = 10, height = 10)
  # 
  # FeaturePlot(scData, features = stem_features, max.cutoff = "q05", cols=c("grey","red"),ncol=4, pt.size = 0.2 )
  # ggsave2(paste0(figures.path,"/",extension,"_featPlot_stem_assay",working_assay,".svg"), width = 14, height = 10)
  # 
  # FeaturePlot(scData, features = effector_features, max.cutoff = "q90", cols=c("grey","red"),ncol=3, pt.size = 0.2 )
  # ggsave2(paste0(figures.path,"/",extension,"_featPlot_effect_assay",working_assay,".svg"), width = 10, height = 7)
  # 
  # FeaturePlot(scData, features = exhaustion_features, max.cutoff = "q80", cols=c("grey","red"),ncol=3, pt.size = 0.2 )
  # ggsave2(paste0(figures.path,"/",extension,"_featPlot_exhaust_assay",working_assay,".svg"), width = 10, height = 7) 
  

  #marker plots
  for (i in 0:(length(table(scData@meta.data[[working_res]]))-1))  {
    tryCatch({
    DefaultAssay(scData) <- working_assay
    Idents(scData) <- working_res
    
    # finds all markers for each cluster
    markers <- FindMarkers(scData, assay = working_assay, ident.1 = i, verbose = FALSE, random.seed = seed)
    markers$gene <- rownames(markers)
    
    write.table(markers, file = paste0(figures.path,"/only_cl_",i,"_allMarkers_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
    
    markers$label <- markers$gene %in% gene.list
    
    markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
    
    
     ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
      data = subset(markers, label == TRUE),
      aes(label = gene),
      size = 7,
      box.padding = unit(3, "lines"),
      point.padding = unit(0.25, "lines"),
      max.overlaps = Inf) +  scale_color_manual(values = c("#7ac5cdff","black","#f08080ff" )) + ggtitle(paste0(disease," CD8 T cells, cluster ",i," vs all clusters; all samples together"))
    
    ggsave(paste0(figures.path,"/volcano_only_cl",i,"_allMarkers_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  

  
  # differential analysis between two conditions
  scData$conditionDisease.clust <- paste(scData$conditionDisease, scData$integrated_snn_res.1,sep = "_")

  markers.path <- paste0(figures.path,"/markers")
  if (!file.exists(markers.path)) {
    dir.create(markers.path,recursive = T)
  }
  
  
  
  DefaultAssay(scData) <- working_assay
  Idents(scData) <- "conditionDisease"
  # finds markers within each cluster differential between the two conditions / replicates
  
  markers <- FindMarkers(scData, assay = working_assay, ident.1 = paste0(condition1), ident.2 = paste0(condition2), verbose = T, random.seed = seed)
  
  markers$gene <- rownames(markers)
  
  write.table(markers, file = paste0(markers.path,"/allclusters_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  

  
  
  for (i in 0:(length(table(scData@meta.data[[working_res]]))-1))  {
    tryCatch({
      DefaultAssay(scData) <- working_assay
      Idents(scData) <- "conditionDisease.clust"
      # finds markers within each cluster differential between the two conditions / replicates
      
      markers <- FindMarkers(scData, assay = working_assay, ident.1 = paste0(condition1,"_",i), ident.2 = paste0(condition2,"_",i), verbose = T, random.seed = seed)
      
      markers$gene <- rownames(markers)
      
      write.table(markers, file = paste0(markers.path,"/cl_",i,"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      
      
      markers$label <- markers$gene %in% gene.list
      
      markers$significant <- ifelse(markers$label == TRUE, "Highlight", ifelse(markers$avg_log2FC < 0, "Decreased signal", "Increased Signal"))
      
      ggplot(markers %>% arrange(label), aes(x = avg_log2FC, y = -log10(p_val)))  + geom_point(aes(color = significant)) + geom_vline(xintercept = 0, col="#333333ff", linewidth=0.7, linetype="dashed") + theme_bw(base_size = 12)  + xlab("Log2 fold change") + ylab("-log10(P value)") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) + geom_text_repel(
        data = subset(markers, label == TRUE),
        aes(label = gene),
        size = 7,
        box.padding = unit(3, "lines"),
        point.padding = unit(0.25, "lines"),
        max.overlaps = Inf) +  scale_color_manual(values = c("#7ac5cdff","black","#f08080ff" )) + ggtitle(paste0(disease," CD8 T cells, cluster ",i," ",condition1,"/",condition2))
      
      ggsave(paste0(markers.path,"/volcano_cl",i,"_",condition1,"_vs_",condition2,"_assay",working_assay,"_",working_res,".svg"), width = 8, height = 4)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  


############## to subset only WT sample and use it to find markers within WT samples

Idents(scData) <- "condition"
scData.subset <- subset(scData, idents = c("WT"))

wt.clusters.anal.path <- paste0(integrated.path,"/wt_clusters_anal")
if (!file.exists(wt.clusters.anal.path)) {
  dir.create(wt.clusters.anal.path,recursive = T)
}
g1 <- list(Stemness=c("Bcl2","Cd28","Cxcr5","Il7r","Sell","Ccr7","Slamf6","Stat3","Zeb1","Lef1","Bcl6","Foxo1","Myb","Id3","Bach2","Tcf7"),
           Exhaustion=c("Gzmf","Gzmc","Gzmb","Ccl5", "Ccl4", "Ccl3", "Entpd1", "Ctla4", "Tigit", "Lag3", "Havcr2", "Pdcd1",
                     "Prdm1", "Batf", "Irf4", "Stat4", "Nfatc2", "Nr4a2", "Nr4a1", "Hif1a", "Egr2", "Id2",
                     "Zeb2", "Tox2", "Tox"),
           DA.combined = c("Cd69","Ccl5","Ccl4","Ccl3","Nt5e","Entpd1","Cd101","Ctla4","Tigit","Lag3","Havcr2","Pdcd1","Gzma","Gzmb","Ifng","Xcl1","Cx3cr1","Cd28","Cxcr5","Il7r","Sell","Ccr7","Slamf6","Tbx21","Prdm1","Batf","Irf4","Stat4","Nfatc2","Nr4a2","Nr4a1","Hif1a","Egr2","Id2","Zeb2","Tox2","Tox","Mki67","Ets1","Ikzf1","Fli1","Stat3","Zeb1","Eomes","Bcl6","Foxo1","Myb","Id3","Bach2","Tcf7"),
           DA.Tpex.WT = c("Tcf7",	"Bach2",	"Id3",	"Myb",	"Foxo1",	"Bcl6",	"Eomes",	"Zeb1",	"Fli1",	"Ikzf1",	"Ets1",	"Stat3"	,"Slamf6"	,"Ccr7"	,"Sell",	"Ifng",	"Il7r",	"Cd69"	,"Cxcr5",	"Pdcd1",	"Lag3",	"Havcr2",	"Entpd1"	,"Mki67",	"Xcl1"	,"Ccr6"	,"Cd28"),
           DA.Ttex.WT = c("Tox","Tox2",	"Prdm1",	"Tbx21","Zeb2",	"Id2","Pdcd1",	"Havcr2",	"Nr4a1",	"Stat4",	"Irf4",	"Batf",	"Egr2",	"Hif1a"	,"Entpd1",	"Nt5e",	"Cd69",	"Ctla4",	"Cd101",	"Tcf7",	"Gzma",	"Mki67",	"Ccl3",	"Ccl4",	"Ccl5",	"Layn",	"Tigit"),
           DA.Tinex.WT = c("Tox",	"Tox2",	"Prdm1",	"Tbx21",	"Zeb2",	"Id2",	"Pdcd1","Lag3","Layn","Gzmb","Cxcr1","Cd69","Havcr2",	"Entpd1"),
           hongbo.general = c("Tox",	"Tcf7",	"Slamf6",	"Sell",	"Havcr2",	"Pdcd1","Entpd1","Cd38","Cd244a","Ifng","Gzma","Gzmb","Itgax"	),
           bence.zeb = c("Klrd1",	"Klrk1",	"Klrc1",	"Klre1",	"Klrg1",	"Zeb2","S1pr5","Ctla4","Pdcd1","Lag3","Coro1a","Cxcr6","Cd101","Tigit","Cx3cr1"	),
           bence.general = c("Pdcd1","Slamf6","Cx3cr1","Ccr7","Sell","Lef1","Klrg1","Ly6c2","Klrb1c","Klrd1","S1pr1","Il7r","Sidt1","Gpr183","Id3","Lgals3","S100a4","Mki67","Gzma","Cd101","Entpd1","Top2a","Mif","S1pr5","Klrc1","Lag3","Ifng","Ccl4","Isg15", "Ifit1", "Isg20"),
           bence.ex = c("Bcl2","Bcl2a1b","Bcl2a1d","Itga4","Itgb7","Itgb1","Cd44","Cxcr6","Ccl3","Cd69","Cd103")
           )

# Stemness = c("Ifng","Xcl1","Cx3cr1","Cd28","Cxcr5","Il7r","Sell","Ccr7","Slamf6","Ets1","Ikzf1","Fli1","Stat3","Zeb1","Eomes","Bcl6","Foxo1","Myb","Id3","Bach2","Tcf7"),
# Exhaustion = c("Ccl5","Ccl4","Ccl3","Nt5e","Entpd1","Cd101","Ctla4","Tigit","Lag3","Havcr2","Pdcd1","Gzma","Prdm1","Batf","Irf4","Stat4","Nfatc2","Nr4a2","Nr4a1","Hif1a","Egr2","Id2","Zeb2","Tox2","Tox"),

g1.for.signature <- list(DA.Tpex.WT = c("Tcf7",	"Bach2",	"Id3",	"Myb",	"Foxo1",	"Bcl6",	"Eomes",	"Zeb1",	"Fli1",	"Ikzf1",	"Ets1",	"Stat3"	,"Slamf6"	,"Ccr7"	,"Sell",	"Ifng",	"Il7r",	"Cxcr5",	"Pdcd1-",	"Lag3-",	"Havcr2-",	"Entpd1-"	,"Mki67",	"Xcl1"	,"Ccr6"	,"Cd28"),
           DA.Ttex.WT = c("Tox",	"Tox2",	"Prdm1",	"Tbx21",	"Zeb2",	"Id2",	"Pdcd1",	"Havcr2",	"Nr4a1",	"Stat4",	"Irf4",	"Batf",	"Egr2",	"Hif1a"	,"Entpd1",	"Nt5e",	"Cd69",	"Ctla4",	"Cd101",	"Tcf7-",	"Gzma",	"Mki67-",	"Ccl3",	"Ccl4",	"Ccl5",	"Layn",	"Tigit"),
           DA.Tinex.WT = c("Tox",	"Tox2",	"Prdm1",	"Tbx21",	"Zeb2",	"Id2",	"Pdcd1","Lag3","Layn","Gzmb","Cxcr1","Cd69-"))


gene.signatures.list <- list(g1=g1.for.signature)


DefaultAssay(scData.subset) <- "RNA"
Idents(scData.subset) <- working_res
markers <- FindAllMarkers(scData.subset, assay = working_assay, random.seed = seed) # only.pos = TRUE

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.7) %>%  dplyr::filter(!grepl(("^Rp[ls]|^mt-|^MT-|^RP[LS]"), gene))  %>% #filters out all ribosomal and mitochondrial genes
  slice_head(n = 10) %>%
  ungroup() -> top10
top10 = top10[!duplicated(top10$gene),]

write.table(markers, file = paste0(wt.clusters.anal.path,"/",extension,"_allMarkers_combined_cond_RNA_",working_res,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

DotPlot(scData.subset, features = top10$gene, group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size')+coord_flip()
ggsave2(paste0(wt.clusters.anal.path,"/",extension,"_WTonly_bubble_allMarkers_no_split_noRiboMt_ClustMarkers_assay",working_assay,"_",working_res,".svg"), width = 7, height = 20)

DotPlot(scData.subset, features = top10$gene,group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "sample") + coord_flip() +RotatedAxis()
ggsave2(paste0(wt.clusters.anal.path,"/",extension,"_WTonly_bubble_allMarkers_split_sample_noRiboMt_ClustMarkers_assay",working_assay,"_",working_res,".svg"), width = 10, height = 20)


for (l in 1:(length(g1)))  {
  DotPlot(scData.subset, features = g1[l],group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size')  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave2(paste0(wt.clusters.anal.path,"/ovarT_integrated","_",names(g1[l]),"_WTonly_bubble_selectMarkers_no_split_assay",working_assay,"_",working_res,".svg"), width = 6, height = 11)
}

# generate the same bubble plot also for non-subsetted all the data
for (l in 1:(length(g1)))  {
  DotPlot(scData, features = g1[l],group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "condition") + coord_flip()   + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave2(paste0(wt.clusters.anal.path,"/ovarT_integrated_CombGenes_ClustIndiv","_",names(g1[l]),"_bubble_selectMarkers_split_Cond_assay",working_assay,"_",working_res,".svg"), width = 7, height = 11)
}



for (i in 1:(length(gene.signatures.list)))  {
  signatures <- AddModuleScore_UCell(scData.subset, features = gene.signatures.list[[i]])
  signature.names <- paste0(names(gene.signatures.list[[i]]), "_UCell")
  
  # generate violin plots
  VlnPlot(signatures, features = signature.names, group.by = working_res)
  ggsave2(paste0(wt.clusters.anal.path,"/signatures",i,"_vln_no-split_",working_res,".svg"), width = 14, height = 5)
  
  # generate feature plots
  FeaturePlot(signatures, features = signature.names,ncol = 3, order = T, pt.size=0.9)
  ggsave2(paste0(wt.clusters.anal.path,"/signatures",i,"_featplot_no-split_",working_res,".svg"), width = 16, height = 5)
}


tz_selected=c("Lag3","Tigit","Entpd1","Tox2","Tox","S100a8","Lcn2","Cd44","Gzmc","Gzmf","Gzmb","Gzmk","Ccl5","Cxcr6","Itga1","Cxcr5","Pdcd1","Bcl6","Ctla4","Nrp1","Tnfsf18","Il2ra","Foxp3","Il23a","Il22","Il17f","Il17a","Il10","Csf2","Prdm1","Rorc","Il13","Il5","Il4","Il6","Tgfb1","Gata3","Tnf","Ifng","Il1b","Tbx21","Il1r1","Ccr7","Sell","Tcf7","Lef1","Bcl11b","Satb1","Runx3","Zbtb7b","Zbtb32","Cd5l","Runx1","Cd4","Cd8b1","Cd8a","Cd3g","Cd3d","Cd3e","Trac","Cd2")

Tfh=c("Klf2","Foxp1","Foxo1","Id3","Prdm1","Stat5","Egr3","Egr2","Batf","Maf","Lef1","Tcf7","Stat4","Stat1","Stat3","Icos","Irf4","Pdcd1","Ascl2")

g0 <- list(Tgeneral=c("Cd4","Cd8b","Cd8a","Cd3g","Cd3f","Cd3e","Trbc","Trac"),
           Th1_tz=c("Tnf","Ifng","Il1b","Tbx21"),
           Th2_tz=c("Il13","Il5","Il4","Il6","Tgfb1","Gata3"),
           Th17_tz=c("Il23a","Il22","Il17f","Il17a","Il10","Csf2","Prdm1","Rorc"),
           Treg_tz=c("Ctla4","Nrp1","Tnfsf18","Il2ra","Foxp3"),
           Tfh_tz=c("Pdcd1","Bcl6","Cxcr5"),
           Trest_tz=c("Il1r1","Ccr7","Sell","Tcf7"),
           Trm_tz=c("Cxcr6","Itga1"),
           Tem_tz=c("Gzmb","Gzmk","Ccl5","Cxcr6","Itga1") #https://www.nature.com/articles/s41467-019-12464-3
)
  
  Th1a=c("Nkg7","Gnly","Gzmh","Ccl5","Ccl4","Gzmb","Gzma","Prf1","Ifng","Gzmm","Klrg1","Zeb2","Cx3cr1","Tnf","Il2rg","Ly6e","Tbx21","Gbp5","Il32","Runx3","Ccl3"),
           Th1b=c("Cxcr4","Pecam1","Pdcd1","Il2rb","Eomes","Klrg1","Cd27","Tnf","Gzmm","Ifng","Nkg7","Gzma","Ccl4","Ccl5","Gzmk"),
           Th1c=c("Gbp5","Samd3","Klrg1","Hopx","Il32","Ifng","Gzmm","Cxcr3","Ccl5","Gzma","Gzmk"),
           Tr1=c("Id2","Irf4","Ly6e","Icos","Gzmm","Cxcr6","Batf","Prdm1","Cxcr3","Tigit","Il21","Cd27","Gzmk","Il2rb","Gzma","Il32","Cd38","Il10","Cd59","Ctla4","Lag3"),
           Treg=c("Icos","Gata3","Ccr4","Il10ra","Prdm1","Il2ra","Batf","Tigit","Cd27","Il32","Sell","Ikzf2","Ctla4","Foxp3"),
           Tcm=c("Icos","Cd27","Areg","Il6st","Sell","Ccd7","Lef1","Tcf7"),
           Th2=c("Il13","Lmo4","Socs3","Ccr10","Tnfsf10","Lgals1","Gata3"),
           Tfhcirc=c("Il7r","Ccr6","Ltb","Aqp3","Tnfsf13b","Cmtm6","Lgals3","Tob1","Zfp36","Bcl6") #https://www.biorxiv.org/content/10.1101/2021.07.21.453277v1.full
)



DimPlot(scData, reduction= "umap", label = T, repel = T, group.by = "SCT_snn_res.1.5", split.by = "diseaseCondition") + ggtitle("Unsupervised clustering")



scData$ConditionClusterName <- paste(scData$condition, scData$integrated_snn_res.1,sep = "_")
Idents(scData) <- "ConditionClusterName"

scData$ClusterNameCondition <- paste(scData$integrated_snn_res.1,scData$condition,sep = "_")
Idents(scData) <- "ClusterNameCondition"

scData$ClusterNameCondition <- factor(scData$ClusterNameCondition , levels = c("0_KO","0_WT","1_KO","1_WT","2_KO","2_WT","3_KO","3_WT","4_KO","4_WT","5_KO","5_WT","6_KO","6_WT","7_KO","7_WT","8_KO","8_WT","9_KO","9_WT","10_KO","10_WT","11_KO","11_WT"))
  
DotPlot(scData, features = tz_selected,group.by = "ClusterNameCondition", assay=working_assay,scale.min = 0,scale.max =10 ,scale.by = 'size') + coord_flip()   + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_gradient2(low = "#005b96", mid = "#FFFFFF", high = "#990000")


ggsave2(paste0(integrated.path,"/bubble_400nFeat-filt_cd4_Markers_split_Cond_assay",working_assay,"_",working_res,".svg"), width = 4, height = 12)



DotPlot(scData.subset, features = c("Ccr7","Il7r","Klf6","Glmap5","Zbtb32","Ikzf2","Tox","Tox2","Id3","Icos","Slamf6","Il21","Izumo1r","Gzmk","Gzmb","Bhlhe40","Cxcr6","Id2","Selplg","Nkg7"),group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =30 ,scale.by = 'size',split.by = "condition") + coord_flip()   + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #https://iiif.elifesciences.org/lax/76339%2Felife-76339-fig3-v2.tif/full/,1500/0/default.jpg


VlnPlot(scData, features = "Cd8a", group.by = working_res,assay=working_assay, split.by = "condition")


i=1
gene.signatures.list[[i]]

signatures <- AddModuleScore_UCell(scData.subset, features = g0)
signature.names <- paste0(names(g0), "_UCell")

# generate violin plots
VlnPlot(signatures, features = signature.names, group.by = working_res)
ggsave2(paste0(wt.clusters.anal.path,"/signatures",i,"_vln_no-split_",working_res,".svg"), width = 14, height = 5)

# generate feature plots
FeaturePlot(signatures, features = signature.names,ncol = 3, order = T, pt.size=0.9)
ggsave2(paste0(wt.clusters.anal.path,"/signatures",i,"_featplot_no-split_",working_res,".svg"), width = 16, height = 5)
}






Idents(scData) <- working_res




#to subset only selected clusters and create dot plots
select.clusters <- c("3","4","5","6","10","11")

for (l in 1:2)  {

  scData.subset <- subset(scData, idents = select.clusters)
  
  DotPlot(scData.subset, features = g1[l],group.by = working_res, assay=working_assay,cols = "RdBu", scale.min = 0,scale.max =70 ,scale.by = 'size',split.by = "condition")  + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="darkgrey", stroke=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave2(paste0(figures.path,"/ovarT_integrated_Clselected_",names(g1[l]),"_",working_res,".svg"), width = 6, height = 6)
  
}


# plot violins with p values
scData$condition.clust <- paste(scData$condition, scData$integrated_snn_res.0.7,sep = "_")
scData.subset <- subset(scData, idents = select.clusters)

# to add p values to violin plots, use this function
violin_plots_p_val <- function(scData, features.to.plot, file_name, comparisons, order.comparisons){
  plot_case1 <- function(features.to.plot, y_max = NULL){
    VlnPlot(scData, features = features.to.plot, cols=c(rep("#ff8686ff",length(comparisons)),rep("white",length(comparisons))),
            pt.size = 0.1,
            group.by = "condition.clust",
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = comparisons, label = "p.format", step.increase = 0, tip.length = 0, method = "t.test", method.args=list(alternative = "two.sided", var.equal = F, paired=FALSE)) + scale_x_discrete(limits=order.comparisons)
    #method="t.test",method.args=list(alternative = "two.sided", var.equal = T, paired=FALSE)
      }
  plot_list <- list()
  y_max_list <- list()
  for (gene in features.to.plot) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, ".svg")
  ggsave(file_name, width = 6, height = 4)
}


# set feature to be plotted
features.to.violinplot.list <- c("Tcf7","Tox","Bach2","Id3","Xcl1","Slamf6","Ikzf1","Fli1","Zeb2","Bcl2","Mcl1")

comparisons <- list(c("KO_3","WT_3"),c("KO_4","WT_4"),c("KO_5","WT_5"),c("KO_6","WT_6"),c("KO_10","WT_10"),c("KO_11","WT_11"))
order.comparisons <- c("KO_5","WT_5","KO_6","WT_6","KO_3","WT_3","KO_11","WT_11","KO_4","WT_4","KO_10","WT_10")
for (w in 1:(length(features.to.violinplot.list)))  {
  features.to.plot <- features.to.violinplot.list[w]
  violin_plots_p_val(scData.subset, features.to.plot, file_name = paste0("violin_",features.to.plot), comparisons, order.comparisons)
}














############## Calculate and plot signature scores
signatures.path <- paste0(figures.path,"/signatures")
if (!file.exists(signatures.path)) {
  dir.create(signatures.path,recursive = T)
}

DefaultAssay(scData) <- "RNA"
Idents(scData) <- working_res

for (i in 1:(length(gene.signatures.list)))  {
  signatures <- AddModuleScore_UCell(scData, features = gene.signatures.list[[i]])
  signature.names <- paste0(names(gene.signatures.list[[i]]), "_UCell")
  
  # generate violin plots
  VlnPlot(signatures, features = signature.names, group.by = working_res)
  ggsave2(paste0(signatures.path,"/signatures",i,"_vln_no-split_",working_res,".svg"), width = 14, height = 5)

  VlnPlot(signatures, features = signature.names, group.by = working_res, split.by = "condition") # 
  ggsave2(paste0(signatures.path,"/signatures",i,"_vln_split_Cond_",working_res,".svg"), width = 20, height = 5)

  VlnPlot(signatures, features = signature.names, group.by = working_res, split.by = "diseaseCondition") # order: melanKO, melanWT, ovarKO, ovarWT
  ggsave2(paste0(signatures.path,"/signatures",i,"_vln_split_disCond_",working_res,".svg"), width = 25, height = 5)
  
  # generate feature plots
  FeaturePlot(signatures, features = signature.names,ncol = 3, order = T, pt.size=0.9)
  ggsave2(paste0(signatures.path,"/signatures",i,"_featplot_no-split_",working_res,".svg"), width = 14, height = 10)

  FeaturePlot(signatures, features = signature.names,split.by = "diseaseCondition", ncol = 3, order = T,pt.size=0.9)
  ggsave2(paste0(signatures.path,"/signatures",i,"_featplot_split-disCond_",working_res,".svg"), width = 17, height = 25)
  
  FeaturePlot(signatures, features = signature.names,split.by = "condition", ncol = 3, order = T,pt.size=0.9)
  ggsave2(paste0(signatures.path,"/signatures",i,"_featplot_split-cond_",working_res,".svg"), width = 10, height = 20)
  
}






















##### Velocity analysis

library(pagoda2)
library(Seurat)
library(SeuratWrappers)
library(reticulate)
library(velocyto.R)
ldat <- read.loom.matrices("/Users/4472414/Projects/Rexamples/Trajectory_Analysis/frozen_pbmc_donor_a_possorted_genome_bam_RJU0T.loom")


# Gather the spliced and unspliced estimates
emat <- ldat$spliced
nmat <- ldat$unspliced

myData <- readRDS("/Users/4472414/Projects/Rexamples/Trajectory_Analysis/seurat_obj.rds")
seurat_obj = myData

# take embedding from the Seurat data object
# NOTE: This assumes you have a seurat data object loaded
# into R memory prior to using this script. STOP and rerun seurat
# pipeline if you do not have this loaded. In my case, my seurat object is simply myData
emb <- myData@reductions$tsne@cell.embeddings

# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='/Users/4472414/Projects/Rexamples/Trajectory_Analysis/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
out_data_dir = "/Users/4472414/Projects/Rexamples/Trajectory_Analysis/"
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='/Users/4472414/Projects/Rexamples/Trajectory_Analysis/pca.csv', quote=F, row.names=F)
write.csv(seurat_obj@reductions$tsne@cell.embeddings, file='/Users/4472414/Projects/Rexamples/Trajectory_Analysis/tsne.csv', quote=F, row.names=F)
write.csv(seurat_obj@reductions$umap@cell.embeddings, file='/Users/4472414/Projects/Rexamples/Trajectory_Analysis/umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='/Users/4472414/Projects/Rexamples/Trajectory_Analysis/gene_names.csv',
  quote=F,row.names=F,col.names=F
)