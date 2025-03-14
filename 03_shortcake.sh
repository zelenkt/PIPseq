#!/bin/bash

#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=15GB
#SBATCH -t 24:00:00
#SBATCH --qos=medium
#SBATCH --job-name=Shortcake
#SBATCH --output=Shortcake.%j.out
#SBATCH --error=Shortcake.%j.err

# run as
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/03_shortcake.sh

# https://mediawiki.moffitt.org/index.php/Running_Basic_Singularity_Containers
# module load singularity # this is not needed

# THIS IS COMPREHENSIVE TOOL to analyze scRNAseq data. Contains all packages one can think of - I will use it to overcome troubles with certificates with scvelo and velocyto 
# https://github.com/rnakato/ShortCake
# cd /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker
# singularity pull docker://rnakato/shortcake:<version>
# they also provide already built sif singularity files on their dropbox so I downloaded them

# there is full and light version and there are multiple python environments that can be checked with micromamba
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif micromamba env list 

# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_light.3.0.0.sif micromamba env list 

# singularity shell /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_light.3.0.0.sif
# the scvelo and velocyto are both in the base environment





# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep1/results_v3.3.0_PIP3M8OTW1_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTW1_10-15-24.txt -o velocyto_wt_ovarian_mouse_tumor_cd8_ot1_rep1_7-23-24 starsolo_sorted_renamedCB_OTW1.bam ../../gencode.vM29.primary_assembly.annotation.gtf

# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep1/results_v3.3.0_PIP3M8OTK1_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTK1_10-15-24.txt -o velocyto_ko_ovarian_mouse_tumor_cd8_ot1_rep1_10-15-24 starsolo_sorted_renamedCB_OTK1.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep2/results_v3.3.0_PIP3M8OTK2_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTK2_10-15-24.txt -o velocyto_ko_ovarian_mouse_tumor_cd8_ot1_rep2_7-23-24 starsolo_sorted_renamedCB_OTK2.bam ../../gencode.vM29.primary_assembly.annotation.gtf

cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep2/results_v3.3.0_PIP3M8OTW2_no_annot
singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTW2_10-15-24.txt -o velocyto_wt_ovarian_mouse_tumor_cd8_ot1_rep2_7-23-24 starsolo_sorted_renamedCB_OTW2.bam ../../gencode.vM29.primary_assembly.annotation.gtf








# original attempts from 7-22-24
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep1/results_v3.3.0_PIP3M8OTK1_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTK1_7-22-24.txt -o velocyto_ko_ovarian_mouse_tumor_cd8_ot1_rep1_7-23-24 starsolo_sorted_renamedCB_OTK1.bam ../../gencode.vM29.primary_assembly.annotation.gtf

# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep2/results_v3.3.0_PIP3M8OTK2_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTK2_7-22-24.txt -o velocyto_ko_ovarian_mouse_tumor_cd8_ot1_rep2_7-23-24 starsolo_sorted_renamedCB_OTK2.bam ../../gencode.vM29.primary_assembly.annotation.gtf

# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep2/results_v3.3.0_PIP3M8OTW2_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTW2_7-22-24.txt -o velocyto_wt_ovarian_mouse_tumor_cd8_ot1_rep2_7-23-24 starsolo_sorted_renamedCB_OTW2.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep1/results_v3.3.0_PIP3M8OTW1_no_annot
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b barcodes.filt.FULL_BARCODE.OTW1_7-22-24.txt -o velocyto_wt_ovarian_mouse_tumor_cd8_ot1_rep1_7-23-24 starsolo_sorted_renamedCB_OTW1.bam ../../gencode.vM29.primary_assembly.annotation.gtf

















# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b raw_matrix/barcodes.tsv -o TZ4_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 starsolo_sorted.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# # here I let the program to generate the sorted bam
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b filtered_matrix/sensitivity_1/barcodes.tsv -o TZ5_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 starsolo_sorted.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# here I used my sorted bam generated by samtools....the results are the same when further inspecting the created loom files so either way should work
# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b filtered_matrix/sensitivity_1/barcodes.tsv -o TZ6_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 starsolo_sorted.bam ../../gencode.vM29.primary_assembly.annotation.gtf






# if all attemps with velocyty fail, we can try https://github.com/minoda-lab/universc universial pipeline which supports pipseq and which modifies it to be able to process it with cell-ranger which should unite the format and make it compatible with other platforms



# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b raw_matrix/barcodes.tsv -o TZ3_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 test.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b OTW1.filt.barcodes.7-17-24.txt -o TZ2_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 test.bam ../../gencode.vM29.primary_assembly.annotation.gtf


# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_full.3.0.0.sif run_env.sh shortcake_default velocyto run -b OTW1.filt.barcodes.7-17-24.txt -o TZ_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 starsolo_sorted.bam ../../gencode.vM29.primary_assembly.annotation.gtf

# singularity exec /share/lab_avram/HPC_Cluster/user/tomas/shortcake-docker/shortcake_light.3.0.0.sif velocyto run -b OTW1.filt.barcodes.7-17-24.txt -o TZ_velocyto_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1 starsolo_sorted.bam ../../gencode.vM29.primary_assembly.annotation.gtf



# velocyto run -b velocyto_barcodes/wt_ovarian_mouse_tumor_cd8_ot1_rep1/barcodes.tsv -o TZ_chr1_TEST_wt_ovarian_mouse_tumor_cd8_ot1_rep1_output wt_ovarian_mouse_tumor_cd8_ot1_rep1/results_PIP3M8OTW1_no_annotation/starsolo_out.chr1.bam gencode.vM31.primary_assembly.annotation.gtf   

# I will run velocyto with filtered barcodes only - I need to extract them from the seurat file first
# t <- seurat_obj@meta.data %>% 
#   rownames_to_column("barcodes") %>%
#   select(barcodes)
# library(stringr)
# z <- as.data.frame(str_split_fixed(t$barcodes, '_', 2))
# colnames(z) <- c("id","barcode")
# table(z$id)
# OTW1.filt.barcodes <- subset(z, id=="OTW1")[,2]

# write.table(OTW1.filt.barcodes, file = paste0("OTW1.filt.barcodes.7-17-24.txt"), sep = "\t", row.names = FALSE, col.names = F,quote=F)

# OTW1.filt.barcodes.7-17-24.txt





# velocyto run \
#     --bcfile /velocyte/C2Routput/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \ # filtered barcode file in gzipped format
#     --outputfolder /velocyte/C2Routput \ # location for the output to go
#     --samtools-threads 18 \ # optional, doesn't increase the speed much
#     --samtools-memory 8000 \ # optional, doesn't increase the speed much
#     -vvv \ # highest verbosity so I can debug if needed
#     /velocyte/C2Routput/outs/possorted_genome_bam.bam \ # BAM file path
#     /velocyte/mm10-transgene-premrna/genes/genes.gtf  # GTF file path - GTF file that was used for creating the cell ranger output



date

