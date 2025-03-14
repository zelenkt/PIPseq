#!/bin/bash

#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=15GB
#SBATCH -t 24:00:00
#SBATCH --qos=medium
#SBATCH --job-name=Pipseeker
#SBATCH --output=Pipseeker.%j.out
#SBATCH --error=Pipseeker.%j.err



# run as
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/pipseeker_docker.sh
# https://mediawiki.moffitt.org/index.php/Running_Basic_Singularity_Containers
# module load singularity # this is not needed


# cd /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/
# singularity pull docker://public.ecr.aws/w3e1n2j6/fluent-pipseeker:3.3.0

# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_3_P11HM4ODWN_UDILibraryIndexmix4_UDILibraryIndexmix4_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq/P11HM4ODWN/
# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_1_PIP8M4ODW_UDILibraryIndexmix1_UDILibraryIndexmix1_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/PIP8M4ODW/
# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_2_P9HM4OC4WR_UDILibraryIndexmix2_UDILibraryIndexmix2_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P9HM4OC4WR/
# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_4_P11M4ODKR_UDILibraryIndexmix5_UDILibraryIndexmix5_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P11M4ODKR/
# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_5_P12HM4OC4WN_UDILibraryIndexmix6_UDILibraryIndexmix6_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P12HM4OC4WN/
# ln -s /share/lab_avram/basespace_Moffitt/15_bulk_sc_rnaseq_ovarian_tbet_dnmt1_10-18-2024/NS3936-DAvram/NS-3936a-DAvram-25B-22HMFTLT4-Lane7/NS-3936_6_P13M4OC4KR_UDILibraryIndexmix7_UDILibraryIndexmix7_L007/* /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P13M4OC4KR/





# GRCm39 (GENCODE vM29 2022.04, Ensembl 106)
# GRCh38.p13 (GENCODE v40 2022.04, Ensembl 106)
# index_path="/share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCh38-2022.04"  # for human


# our first PIP-seq runs were prepared with V4 chemistry. Starting with the recipient cells with T10 we will have V chemistry. The chemistry needs to be specified when running the pipeline as there are differences how the reads were assembled.


# # Please choose from: full, cells, barcode, feature, buildmapref, buildannotref, merge, extract
# use this flag --run-barnyard when using combination of mouse and human data


# these sample were only from mouse (select the proper index_path above)
# sample="PIP8M4ODW"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/PIP8M4ODW

# sample="P11M4ODKR"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P11M4ODKR

# sample="P13M4OC4KR"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P13M4OC4KR


# # to perform the main full run with mouse only data
# index_path="/share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCm39-2022.04"  # for mouse
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif full --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot --save-svg --dpi 400 --random-seed 42 --resume-last-run --skip-version-check

# # this option is not compatible with the full run so I need to separately run the barcode option after running the full run - this is required for some downstream analyses such as for scVelo
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif barcode --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot --random-seed 42 --sorted-bam --resume-last-run --skip-version-check




# these sample were from mouse combined with human cells so they have to be run with the --run-barnyard flag (select the proper index_path above)
# first we need to run it together in order to determine which cells are human and which mouse - after we rerun these samples with either human or mouse mapping and then filter the barcodes of interest
# sample="P11HM4ODWN" 
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P11HM4ODWN

# sample="P9HM4OC4WR"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P9HM4OC4WR

# sample="P12HM4OC4WN"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P12HM4OC4WN

# # to perform the main full run with combined mouse+human data
# index_path="/share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCh38-and-GRCm39-2022.04" # for mouse+human

# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif full --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot --save-svg --dpi 400 --random-seed 42 --resume-last-run --skip-version-check --run-barnyard

# # this option is not compatible with the full run so I need to separately run the barcode option after running the full run
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif barcode --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot --random-seed 42 --sorted-bam --resume-last-run --skip-version-check






# re-run barnyard (mouse+human) samples with separately human and mouse annotations
sample="P11HM4ODWN"
cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P11HM4ODWN

# sample="P9HM4OC4WR"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P9HM4OC4WR

# sample="P12HM4OC4WN"
# cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/P12HM4OC4WN

# index_path="/share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCh38-2022.04"  # for human
# # testing an option with --retain-barcoded-fastqs that could be used for further filtering of mouse vs human cells (maybe)
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif full --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot_human_mapping --save-svg --dpi 400 --random-seed 42 --resume-last-run --skip-version-check

# # this option is not compatible with the full run so I need to separately run the barcode option after running the full run
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif barcode --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot_human_mapping --random-seed 42 --sorted-bam --resume-last-run --skip-version-check



index_path="/share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCm39-2022.04"  # for mouse
# testing an option with --retain-barcoded-fastqs that could be used for further filtering of mouse vs human cells (maybe)
singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif full --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot_mouse_mapping --save-svg --dpi 400 --random-seed 42 --resume-last-run --skip-version-check

# this option is not compatible with the full run so I need to separately run the barcode option after running the full run
singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif barcode --chemistry V --fastq $sample --star-index-path $index_path --output-path results_v3.3.0_"$sample"_no_annot_mouse_mapping --random-seed 42 --sorted-bam --resume-last-run --skip-version-check

























# FOLLOWING ARE COMMANDS USED FOR PREVIOUS PIPSEQ RUNS

# sample="PIP3M8OTK1"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep1

# sample="PIP3M8OTK2"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep2


# sample="PIP3M8OTW2"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep2


# sample="PIP3M8OTW1"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_ovarian_mouse_tumor_cd8_ot1_rep1



# sample="PIP2M8MW1"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_melanoma_mouse_tumor_cd8_pmel_rep1

# sample="PIP2M8MW2"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_melanoma_mouse_tumor_cd8_pmel_rep2


# sample="PIP2M8MK1"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_melanoma_mouse_tumor_cd8_pmel_rep1


# sample="PIP4M8MK2"
# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_melanoma_mouse_tumor_cd8_pmel_rep2


# # Please choose from: full, cells, barcode, feature, buildmapref, buildannotref, merge, extract

# to perform the main full run
# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif full --chemistry v4 --fastq $sample --star-index-path /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCm39-2022.04 --output-path results_v3.3.0_"$sample"_no_annot --save-svg --dpi 400 --random-seed 42 --resume-last-run --skip-version-check



# The BAM file produced by STAR during the PIPseeker analysis is not sorted by genome position. Some users may want to have a sorted BAM compatible with various downstream applications, such as velocyto and scvelo. Running PIPseeker barcode with the --sorted-bam flag will generate a sorted BAM file (not used by PIPseeker), that can serve as an input to third party applications. The barcode and MI information will be incorporated as tags in the BAM file. Note that in PIPseq V, the MI will be a concatenation of the binning index and the IMI for each read. Using it in place of a UMI in the above applications may lead to inaccurate results for certain cell types.


# this option is not compatible with the full run so I need to separately run the barcode option after running the full run

# singularity run /share/lab_avram/HPC_Cluster/user/tomas/pipseeker-docker/fluent-pipseeker_3.3.0.sif barcode --chemistry v4 --fastq $sample --star-index-path /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/references/pipseeker-gex-reference-GRCm39-2022.04 --output-path results_v3.3.0_"$sample"_no_annot --random-seed 42 --sorted-bam --resume-last-run --skip-version-check


# it was always stuck in one original step so I'm trying these options to skip the original checks, the first version-check skipping solved the problem
# --skip-version-check
# --skip-preflight-check











# # to be able to distinguish different samples upon merging into a single loom file, append sample name prefix to the cell barcodes, both in the bam file and in the barcode.tsv. The python script requires argparse, re, pysam so this needs to be loaded - I have them in tobias env
# ml purge
# ml Anaconda3/2024.02-1

# source activate /home/4476125/tobias_tz


# # sample="MTK2"

# # cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_melanoma_mouse_tumor_cd8_pmel_rep2/results_v3.3.0_PIP4M8MK2_no_annot


# sample="MTW2"

# cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/wt_melanoma_mouse_tumor_cd8_pmel_rep2/results_v3.3.0_PIP2M8MW2_no_annot

# # append to beginning - USED
# python /share/lab_avram/HPC_Cluster/user/tomas/bin/bamtagregex.py starsolo_sorted.bam starsolo_sorted_renamedCB_"$sample".bam --tag CB --pattern ^ --replace "${sample}_"

# # # # append to the end - NOT USED
# # # python bamtagregex.py starsolo_sorted.chr19.bam zrenamed_starsolo_sorted.chr19.bam --tag CB --pattern $ --replace "test"

# ml SAMtools/1.16.1-GCC-11.3.0

# samtools index starsolo_sorted_renamedCB_"$sample".bam
# samtools sort -t CB -O BAM -o cellsorted_starsolo_sorted_renamedCB_"$sample".bam starsolo_sorted_renamedCB_"$sample".bam
# samtools index cellsorted_starsolo_sorted_renamedCB_"$sample".bam




















date