#!/usr/bin/env Rscript

suppressMessages(
if  (!require("reshape")) {
    install.packages("reshape", dependencies = TRUE)
    library(reshape)
    }
)

suppressMessages(
if  (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE)
    library(reshape2)
    }
)

suppressMessages(
if  (!require("tidyverse")) {
    install.packages("tidyverse", dependencies = TRUE)
    library(tidyverse)
    }
)

suppressMessages(
if  (!require("plyr")) {
    install.packages("plyr", dependencies = TRUE)
    library(plyr)
    }
)

suppressMessages(
if  (!require("dplyr")) {
    install.packages("dplyr", dependencies = TRUE)
    library(dplyr)
    }
)

suppressMessages(
if  (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE)
    library(ggplot2)
    }
)



# setwd("/share/lab_avram/HPC_Cluster/user/tomas/15_PIPseq8-13_10-18-2024/test")
# file_to_analyze <- "subset1000_raw_counts_sub_withNAs_nt_recip_bcl11b_cd4ercre_noTamoxifen_huTILs_p84.txt"

Project_Names <- c("nt_recip_bcl11b_cd4ercre_noTamoxifen_huTILs_p84","nt_recip_dnmt3_tbetcre_huTILs_p41","wt_recip_bcl11b_cd4ercre_huTILs_p19")

multiplet_threshold <- 85  # % of human or mouse transcript needed to annotate a cell as human or mouse


for (i in 1:length(Project_Names)){
    
  Project_Name <- Project_Names[i]

  file_to_analyze <- paste0("raw_counts_sub_withNAs_",Project_Name,".txt")

  raw_counts.sub <- read.table(paste0(file_to_analyze))

  melt_raw_counts.sub <- melt(raw_counts.sub, id = "species")  %>% drop_na(value)

  # write.table(melt_raw_counts.sub, file = paste0("melt_raw_counts_sub_",Project_Name,".txt"),sep = "\t", row.names = TRUE, col.names = TRUE)



  # function needed
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE, .drop=TRUE) {
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                  .fun = function(xx, col) {
                    c(N    = length2(xx[[col]], na.rm=na.rm),
                      mean = mean   (xx[[col]], na.rm=na.rm)
                    )
                  },
                  measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    return(datac)
  }

  sum_melt_raw_counts.sub <- summarySE(melt_raw_counts.sub,measurevar = "value",groupvars = c("species","variable"))
  # sum_melt_raw_counts.sub <- sum_melt_raw_counts.sub[,-4]``
  cast_sum_melt_raw_counts.sub = dcast(sum_melt_raw_counts.sub, variable~species, sum, value.var="N") 

  cast_sum_melt_raw_counts.sub$percHuman <- cast_sum_melt_raw_counts.sub$human/(cast_sum_melt_raw_counts.sub$human+cast_sum_melt_raw_counts.sub$mouse)*100
  cast_sum_melt_raw_counts.sub$percMouse <- cast_sum_melt_raw_counts.sub$mouse/(cast_sum_melt_raw_counts.sub$human+cast_sum_melt_raw_counts.sub$mouse)*100

  cast_sum_melt_raw_counts.sub$annotation <- ifelse(cast_sum_melt_raw_counts.sub$percHuman > multiplet_threshold, "human", ifelse(cast_sum_melt_raw_counts.sub$percMouse > multiplet_threshold, "mouse", "multiplets"))

  write.table(cast_sum_melt_raw_counts.sub, file = paste0("barnyard_metrics_",Project_Name,".txt"),sep = "\t", row.names = FALSE, col.names = TRUE)

  sum.percent <- paste(paste0("human=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["human"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%"),"\n",paste0("mouse=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["mouse"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%"), "\n",paste0("multiplets=",round(100*(table(cast_sum_melt_raw_counts.sub$annotation)["multiplets"] / (table(cast_sum_melt_raw_counts.sub$annotation)[["human"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["mouse"]]+table(cast_sum_melt_raw_counts.sub$annotation)[["multiplets"]])), digits = 1),"%") )


  fit <- ggplot(cast_sum_melt_raw_counts.sub, aes(x = human, y = mouse))  + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed")  + geom_point(aes(color = annotation)) + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey")) + ggtitle(paste0("Barnyard metrics ",Project_Name)) + scale_x_continuous(limits = c(0,10000), breaks = c(0,2000,6000,10000)) +scale_y_continuous(limits = c(0,10000), breaks = seq(0,10000, by = 2000)) + annotate(geom="text", x=8000, y=9000, label=sum.percent ,color="black",size=5)
  svg(paste0("barnyard_dotplot_",Project_Name,".svg"), width = 6, height = 4, family ="Arial")
  print(fit)
  dev.off()


  fit <- ggplot(cast_sum_melt_raw_counts.sub, aes(x = log10(human), y =log10(mouse))) + geom_abline(intercept = 0, col="grey", linewidth=0.7, linetype="dashed") + geom_point(aes(color = annotation)) + theme_bw(base_size = 12)  + xlab("# Human transcripts") + ylab("# Mouse transcripts") + theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(color="black",size=20) ,axis.title.y = element_text(color="black",size=20), axis.text.x = element_text(color="black",size=20), axis.text.y = element_text(color="black",size=20)) + theme(legend.text = element_text(color="black",size = 18), legend.title = element_text(color="black",size = 18)) +  scale_color_manual(values = c("#7ac5cdff","#f08080ff","darkgrey")) + ggtitle(paste0("Barnyard metrics ",Project_Name))+ scale_x_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) +scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5, by = 1)) + annotate(geom="text", x=1, y=1, label=sum.percent ,color="black",size=5)
  svg(paste0("barnyard_log10_dotplot_",Project_Name,".svg"), width = 6, height = 4, family ="Arial")
  print(fit)
  dev.off()


}