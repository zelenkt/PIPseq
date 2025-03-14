#newly installed version of trust4 on HPC
ml TRUST4/1.0.13-r474
grabnode

# or run:
sbatch ~/bin/trust4.sh

cd /share/lab_avram/HPC_Cluster/user/tomas/10_PIPseq2-4_scRNAseq_CUTRUN_CD8_melanoma_ovarian_mouse_10-19-2023/ko_ovarian_mouse_tumor_cd8_ot1_rep1/results_v3.3.0_PIP3M8OTK1_no_annot/



results_v3.3.0_PIP3M8OTK1_no_annot




# when using bam files:
run-trust4 -f ~/lib/GRCm38_bcrtcr.fa --ref ~/lib/mouse_IMGT+C.fa -b starsolo_out.bam 

# when using the raw sequences:
run-trust4 -f ~/lib/GRCm38_bcrtcr.fa --ref ~/lib/mouse_IMGT+C.fa -1 S01_0_2_2_009212022_S1_L001_R1_001.fastq.gz -2 S01_0_2_2_009212022_S1_L001_R2_001.fastq.gz

# for human
-f ~/lib/hg38_bcrtcr.fa --ref ~/lib/human_IMGT+C.fa


#note that the version of TRUST4 on HPC is old and not supporting standardized airr format https://github.com/liulab-dfci/TRUST4/issues/96
#the AIRR format can be used with many other programs https://docs.airr-community.org/en/latest/resources/rearrangement_support.html
#so I had to rerun it from my mac using newer version using the following command with flag --stage 3:

#however I noticed that this way some of the information - particularly info on the  TCR sequence is lost so it is better to extract all the information from scratch
conda activate tcr
run-trust4 -f /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/GRCm38_bcrtcr.fa --ref /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/mouse_IMGT+C.fa -b starsolo_out.bam

https://www.researchgate.net/post/How-can-I-find-the-sequence-of-the-transgenic-T-cell-receptor-of-the-OT-1-mouse
https://www.informatics.jax.org/allele/MGI:3054907#TransgeneDescription
here they describe what OT1 is...however based on our results, I believe our OT1 sequence might be different
https://www.jax.org/strain/003831
https://www.addgene.org/52111/ here is some sequence




on HPC I can run the following
grabnode

ml TRUST4
run-trust4 -f ~/lib/GRCm38_bcrtcr.fa --ref ~/lib/mouse_IMGT+C.fa -b starsolo_out.bam 

run-trust4 -f /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/GRCm38_bcrtcr.fa --ref /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/mouse_IMGT+C.fa -b starsolo_out.bam  --stage 3 # this is to be rerun on HPC output to generate additional airre file using the newer local version of trust4
