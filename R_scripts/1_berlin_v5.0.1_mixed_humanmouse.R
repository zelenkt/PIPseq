#!/usr/bin/env Rscript

#can be run directly from command line by: Rscript 1-Single_Cell_RNAseq_Analysis.R

####---- User Input ----####
## Set local Github repository as working directory
setwd("/Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN")


# This script will only pre-process the mixed samples with combined human and mouse scRNAseq data

organism <- "MouseHuman"

# Specify the input and output folder paths
input_path <- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/1_mixed_human_mouse_samples/Input/pipseeker_v3.3"
output_path<- "2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/1_mixed_human_mouse_samples/Output"

Project_Names <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen","nt_recip_dnmt3_tbetcre","wt_recip_bcl11b_cd4ercre")







# Define the file names for count file and/or for raw matrix data from CellRanger / PIPseeker
Count_file_name <- "Count_File.txt.gz"
Matrix_name <- "matrix.mtx.gz"
Barcodes_name <- "barcodes.tsv.gz"
Features_name <- "features.tsv.gz"

# Define the file name for clinical/meta data (OPTIONAL)
meta_name <- "Meta_File_Example_GSE116256_AML921A-D0.txt.gz"


seed <- 42


####---- Load Packages ----####

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix", "fields","SeuratDisk","tibble","SeuratData","SingleR","SingleCellExperiment", "ggplot2", "cowplot", "ggrepel", "UCell", "ggpubr","stringr","SeuratWrappers","reticulate","reshape2","reshape","tidyverse","scCustomize","Nebulosa")
invisible(lapply(packages, library, character.only = TRUE))

set.seed(seed)




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


