#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=80GB
#SBATCH --qos=medium
#SBATCH --time=72:00:00
#SBATCH --job-name=melt
#SBATCH --output=melt.%j.out
#SBATCH --error=melt.%j.err


# # to run
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/02_convert_barcodes_meltData.sh



##### LOAD MODULES #####
# ml Anaconda3/2024.02-1
# conda create --prefix /home/4476125/CondaR421 python=3.6
# source activate /home/4476125/CondaR421
# conda remove -p /home/4476125/CondaR421 --all



ml R/4.4.1-gfbf-2023b
export R_LIBS_USER='/home/4476125/Rmelt'

# modify the path to the directory containing raw_counts_sub_withNAs_nt_recip_bcl11b_cd4ercre_noTamoxifen_huTILs_p84.txt file for each mixed human-mouse sample
cd /share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/BERLIN/2_recipients_mouse_human_TILs_PIPseq8-13_10-18-2024/1_mixed_human_mouse_samples/Output/extract_barnyard_plots_and_species_barcodes

# you may need to adapt the sample names also in the R script here:
cp /share/lab_avram/HPC_Cluster/user/tomas/bin/scRNAseq_pipeline_pipseeker_v3.3.0/meltData.R ./

time Rscript meltData.R

rm meltData.R

date

