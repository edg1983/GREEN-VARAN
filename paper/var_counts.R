library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
source("Read_VCF_stats.R")

# Variables -----------------
data_folder <- "data/var_counts/"
x_axis_text_45 <- theme(axis.text.x = element_text(angle=45, hjust=1))
chr_order <- paste("chr",c(1:22,"X","Y","M"), sep = "")
data_tranches <- c(
  "gnomAD/1000G AF < 0.01\ncohort AF < 0.1",
  "overlap GREEN-DB region",
  "overlap TFBS/DNase/UCNE",
  "FATHMM/ncER/ReMM >= FDR50 threshold",
  "region constraint >= 0.7"
)

# Summarize variants count per individual ------------------
data <- NULL
data_tranches_labels <- c()
group_label <- ""
files <- list.files(data_folder,pattern = ".txt")
for (f in files) {
  step <- as.numeric(gsub("_.+","",f))
  df <- readVCFstat(paste0(data_folder, f))$PSC
  
  if (step == 0) {
    df$group <- "rare"
    group_label <- paste0("rare: ",data_tranches[step+1]) 
    df$group_label <- group_label
  } else if (step == 1) {
    df$group <- paste0("level",step)
    group_label <- paste0("level",step,": ",data_tranches[step+1])
    df$group_label <- group_label
  } else {
    df$group <- paste0("level",step)
    group_label <- paste0("level",step,": previous step +\n",data_tranches[step+1])
    df$group_label <- group_label
  }
  data <- rbind(data,df)  
  data_tranches_labels <- c(data_tranches_labels, group_label)
}
data$group <- factor(data$group,levels=c("rare",paste0("level",1:4)))
data$group_label <- factor(data$group_label,levels=data_tranches_labels)
data$n_vars <- data$nNonRefHom + data$nHets + data$nIndels

# Plot ----------------
p <- ggplot(data, aes(x=group, y=n_vars, color=group_label)) + 
  geom_quasirandom(bandwidth=5, size=1, shape=1, alpha=.5, varwidth=TRUE) + 
  scale_y_log10() + 
  theme_bw() + theme(legend.position="right", legend.direction = "vertical", legend.key.size = unit(22,"point")) + 
  x_axis_text_45 + labs(y="N vars",x="Variant tranches", color="Tranches")
ggsave(p, filename = "plots/Supp Figure XX - Trio analysis var counts per sample.png", height=6, width=8, dpi=150, device="png")


