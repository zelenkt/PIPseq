
##manual investigation of tcrs - I used the trust4 report file from which I filtered out out-of-frame and non-completed tcrs...and I only used TCRa or TCRb receptors
library(ggplot2)
library(dplyr)

ovar_wo_outframe <- read_delim("ovar_tcr_modif_reports_without_out-frame.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

set.seed(42) # just to make it reproducible

palette <- sample(colorRampPalette(RColorBrewer::brewer.pal(12,name = 'Paired'))(208))


ggplot(ovar_wo_outframe, aes(x = "", y=perc_inFrame_sum, fill = reorder(CDR3aa, perc_inFrame_sum)))+geom_bar(stat = "identity") +coord_polar("y")+ facet_wrap(TCR~sample2,ncol=4) + theme_void() +theme(legend.position="none") +geom_text(aes(label = scales::percent(round(perc_inFrame_sum,2))), position = position_stack(vjust = 0.5)) + scale_fill_manual(values = palette)

ggsave("ovar_wo_outframe_wPerc.svg", width = 7, height = 4.5)


ggplot(ovar_wo_outframe, aes(x = "", y=perc_inFrame_sum, fill = reorder(CDR3aa, perc_inFrame_sum)))+geom_bar(stat = "identity") +coord_polar("y")+ facet_wrap(TCR~sample2,ncol=4) + theme_void() +theme(legend.position="none") + scale_fill_manual(values = palette)

ggsave("ovar_wo_outframe_woPerc.svg", width = 7, height = 4.5)





# from here below is for testing






ggplot(ovar_wo_outframe, aes(x = "", y=perc_inFrame_sum, fill = reorder(CDR3aa, perc_inFrame_sum)))  +geom_bar(stat = "identity", width=1, color="white") +coord_polar("y")+ facet_wrap(TCR~sample2,ncol=4) + theme_void() +theme(legend.position="none")+geom_col()


ovar_wo_outframe %>%
  group_by(sample2) %>%
  summarise(volume = sum(Cnt)) %>%
  mutate(share=volume/sum(volume)) %>%
  ggplot(aes(x="", y= perc_inFrame_sum, fill=reorder(CDR3aa, perc_inFrame_sum))) +
  geom_col() +
  geom_text(aes(label = scales::percent(round(share,3))), position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") + 
  theme_void()

df %>%
  group_by(Make) %>%
  summarise(volume = sum(Cnt)) %>%
  mutate(share=volume/sum(volume)) %>%
  ggplot(aes(x="", y= share, fill=reorder(Make, volume))) +
  geom_col() +
  geom_text(aes(label = scales::percent(round(share,3))), position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") + 
  theme_void()




ggplot(ovar_wo_outframe, aes(x = sample, y=sum, fill = sample))+geom_bar(stat = "identity", color="black") + facet_wrap(~TCR, scales = "free") 

+theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442", "#6A51A3", "#CC79A7")) + ggtitle("Full-length human TCR, read counts")+ ylab("Read count")


ggsave("tcrab_human_read_counts.svg", width = 7, height = 4.5)












TCRaTCRbMouse <- read_delim("TCRaTCRbMouse.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

TCRaTCRbHumanIncomplete <- read_delim("TCRaTCRbHumanIncomplete.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

TCRaTCRbMouseIncomplete <- read_delim("TCRaTCRbMouseIncomplete.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)






ggplot(TCRaTCRbHuman, aes(x = patient, y=counts, fill = sample))+geom_bar(stat = "identity", color="black") + facet_wrap(~tcr2, scales = "free") +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442", "#6A51A3", "#CC79A7")) + ggtitle("Full-length human TCR, read counts")+ ylab("Read count")
ggsave("tcrab_human_read_counts.svg", width = 7, height = 4.5)

ggplot(TCRaTCRbHuman, aes(x = patient, fill = sample))+geom_bar() + facet_wrap(~tcr2) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442", "#6A51A3", "#CC79A7"))+ ggtitle("Full-length human TCR, clonotype counts")+ ylab("Clonotype count")
ggsave("tcrab_human_clonotype_counts.svg",  width = 6.5, height = 4.5)


ggplot(TCRaTCRbHumanIncomplete, aes(x = patient, y=counts, fill = sample))+geom_bar(stat = "identity", color="black") + facet_wrap(~tcr2, scales = "free") +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442", "#6A51A3", "#CC79A7")) + ggtitle("Human TCRs including incomplete, read counts")+ ylab("Read count")
ggsave("tcrab_human_read_counts_incl_incomplete.svg", width = 7, height = 4.5)

ggplot(TCRaTCRbHumanIncomplete, aes(x = patient, fill = sample))+geom_bar() + facet_wrap(~tcr2) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442", "#6A51A3", "#CC79A7"))+ ggtitle("Human TCRs including incomplete, clonotype counts")+ ylab("Clonotype count")
ggsave("tcrab_human_clonotype_counts_incl_incomplete.svg",  width = 6.5, height = 4.5)

#mouse
ggplot(TCRaTCRbMouse, aes(x = patient, y=counts, fill = sample))+geom_bar(stat = "identity", color="black") + facet_wrap(~tcr1, scales = "free") +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442","#999999","#D55E00", "#6A51A3", "#CC79A7")) + ggtitle("Full-length mouse TCR, read counts")+ ylab("Read count")
ggsave("tcrab_mouse_read_counts.svg", width = 10, height = 5)


ggplot(TCRaTCRbMouse, aes(x = patient, fill = sample))+geom_bar() + facet_wrap(~tcr1) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442","#999999","#D55E00", "#6A51A3", "#CC79A7")) + ggtitle("Full-length mouse TCR, clonotype counts")+ ylab("Clonotype count")
ggsave("tcrab_mouse_clonotype_counts.svg",  width = 9, height = 5)


ggplot(TCRaTCRbMouseIncomplete, aes(x = patient, y=counts, fill = sample))+geom_bar(stat = "identity", color="black") + facet_wrap(~tcr1, scales = "free") +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442","#999999","#D55E00", "#6A51A3", "#CC79A7")) + ggtitle("Mouse TCRs including incomplete, read counts")+ ylab("Read count")
ggsave("tcrab_mouse_read_counts_incl_incomplete.svg", width = 10, height = 5)


ggplot(TCRaTCRbMouseIncomplete, aes(x = patient, fill = sample))+geom_bar() + facet_wrap(~tcr1) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + scale_fill_manual(values = c("#0072B2", "#56B4E9","darkgreen" , "#009E73", "#E69F00",  "#F0E442","#999999","#D55E00", "#6A51A3", "#CC79A7")) + ggtitle("Mouse TCRs including incomplete, clonotype counts")+ ylab("Clonotype count")
ggsave("tcrab_mouse_clonotype_counts_incl_incomplete.svg",  width = 9, height = 5)


"#D55E00",
"#999999",
  c("#DADAEB", "#9E9AC8", "#6A51A3")






library(abdiv)
library('rlist')
library(tidyr)

x <- c(1, 1,1, 1, 1, 1)
y <- c(0.2,0.2,0.244,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
## to normalize the shannon index to be comparable between samples I have to divide it by the count of occurrences
shannon(x) / log(richness(x))
shannon(y) / log(richness(y))

### for mouse, manual shannon index
TCRaTCRbMouseWide <- subset(TCRaTCRbMouseIncomplete, tcr2=="TRB")
shannon_names <- c()
shannon_results <- c()
TCRaTCRbMouseWide <- TCRaTCRbMouseWide %>% spread(key = sample, value = counts)
TCRaTCRbMouseWide <- TCRaTCRbMouseWide[,8:17]
for (i in 1:length(colnames(TCRaTCRbMouseWide))){
test <- as.vector(TCRaTCRbMouseWide[,i] %>% drop_na())
shannon <- shannon(test[[1]])/ log(richness(test[[1]]))
shannon_names <- list.append(shannon_names, names(test))
shannon_results <- list.append(shannon_results, shannon)
}
shannon_mouse_final <- as.data.frame(shannon_results,shannon_names) %>% cbind(shannon_names, c("peritoneum","peritoneum","peritoneum","peritoneum","tumor", "tumor","tumor","tumor","tumor","tumor"), c("KO","KO","WT","WT","KO", "KO","KO","WT","WT","WT"))
colnames(shannon_mouse_final)[2] = "sample"
colnames(shannon_mouse_final)[3] = "type"
colnames(shannon_mouse_final)[4] = "gen"
shannon_mouse_final[is.na(shannon_mouse_final)] <- 0
ggplot(shannon_mouse_final, aes(x = type,y=shannon_results, group = gen))+geom_jitter (aes(col =gen), size=3, width = 0.05, height=0) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + ylim(0,1) + ylab("Normalized Shannon index")+ ggtitle("Mouse TCRb including incomplete")

+ ggtitle("Full-length mouse TCRa") + ggtitle("Mouse TCRs including incomplete")
ggsave("tcrb_shannon_mouse_TCR_incl_incomplete.svg",  width = 5, height = 4)



### for humans, manual shannon index
TCRaTCRbHumanWide <- subset(TCRaTCRbHumanIncomplete, tcr2=="TRA")
shannon_names <- c()
shannon_results <- c()
TCRaTCRbHumanWide <- TCRaTCRbHumanWide %>% spread(key = sample, value = counts)
TCRaTCRbHumanWide <- TCRaTCRbHumanWide[,8:15]
for (i in 1:length(colnames(TCRaTCRbHumanWide))){
  test <- as.vector(TCRaTCRbHumanWide[,i] %>% drop_na())
  shannon <- shannon(test[[1]])/ log(richness(test[[1]]))
  shannon_names <- list.append(shannon_names, names(test))
  shannon_results <- list.append(shannon_results, shannon)
}
shannon_human_final <- as.data.frame(shannon_results,shannon_names) %>% cbind(shannon_names, c("patient 19","patient 19","patient 19","patient 19","patient 41", "patient 41","patient 41","patient 41"), c("KO","KO","WT","WT","KO", "KO","WT","WT"))
colnames(shannon_human_final)[2] = "sample"
colnames(shannon_human_final)[3] = "type"
colnames(shannon_human_final)[4] = "gen"
shannon_human_final[is.na(shannon_human_final)] <- 0
ggplot(shannon_human_final, aes(x = type,y=shannon_results, group = gen))+geom_jitter (aes(col =gen), size=3, width = 0.05, height=0) +theme_bw() + theme(legend.position = "right",strip.text.x = element_text(size=20,  family="Arial"), strip.background = element_rect(colour="white", fill="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank() ,axis.title.y = element_text(color="black",size=20,  family="Arial"), axis.text.x = element_text(color="black",size=20,  family="Arial"), axis.text.y = element_text(color="black",size=20,  family="Arial")) + theme(legend.text = element_text(color="black",size = 18,  family="Arial"), legend.title = element_text(color="black",size = 18,  family="Arial")) + ylim(0,1) + ylab("Normalized Shannon index")+ ggtitle("Human TCRa including incomplete")
+ ggtitle("Full-length human TCRb")

ggsave("tcra_shannon_human_TCR_incl_incomplete.svg",  width = 5, height = 4)

