library(ROCR)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Variables ----------------
scores <- c("CADD","DANN","ExPECTO","LinSight","NCBoost","ReMM","FIRE",
            "GenoCanyon","GenoSkylinePlus","GWAVA_1","GWAVA_2","GWAVA_3",
            "EIGEN_PHRED","EIGEN_PC_PHRED","FATHMM_MKL_NC",
            "FATHMM_XF","ncER_perc","FINSURF")
data_folder <- "data/scores_comparison/"

# Load data -----------------------
#Curated TP and TN vars from NCBoost paper
TN <- read.csv(paste0(data_folder, "TN_set.18scores.csv"), header=T, sep=" ", stringsAsFactors = F, na.strings = '.')
TP <- read.csv(paste0(data_folder, "TP_set.18scores.csv"), header=T, sep=" ", stringsAsFactors = F, na.strings = '.')
colnames(TN) <- gsub("X\\.+[0-9]+\\.","",colnames(TN))
colnames(TP) <- gsub("X\\.+[0-9]+\\.","",colnames(TP))
TN$Class <- "TN"
TP$Class <- "TP"
scores_df <- rbind(TP,TN)
scores_df$Class <- factor(scores_df$Class, levels=c("TP","TN"))
scores_df <- scores_df[,c("CHROM","POS","REF","ALT",scores,"Class")]
scores_df$ExPECTO <- apply(scores_df,1, function(x) max(abs(as.numeric(unlist(strsplit(x["ExPECTO"], split = ","))))))

# Performance evaluation -----------------------
#OPM is a summary metric developed to capture the overall performance of a classifier 
#OPM defined in: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0117380
#((PPV+NPV)(sens+spec)(acc+nMCC)) / 8
#nMCC = (1+MCC) / 2
preds <- list()
for (s in scores) {
  df <- scores_df %>% select(s,Class) %>% filter(!is.na(!!as.name(s)))
  message(s," ",nrow(df))
  preds[[s]] <- prediction(df[,s], df$Class)
}

ROC_curves <- list()
PrecRecall_curves <- list()
AUC_values <- list()
results <- list()

for (s in scores) {
  ROC_curves[[s]] <- performance(preds[[s]], measure="tpr", x.measure="fpr")
  PrecRecall_curves[[s]] <- performance(preds[[s]], measure="prec", x.measure="rec")
  
  AUC_values[[s]] <- performance(preds[[s]], measure="auc")@y.values[[1]]
  
  F_values <- performance(preds[[s]], measure="f")@y.values[[1]]
  acc_values <- performance(preds[[s]], measure="acc")@y.values[[1]]
  MCC_values <- performance(preds[[s]], measure="mat")@y.values[[1]]
  NPV_values <- performance(preds[[s]], measure="npv")@y.values[[1]]
  PPV_values <- performance(preds[[s]], measure="ppv")@y.values[[1]]
  tpr_values <- performance(preds[[s]], measure="tpr")@y.values[[1]]
  tnr_values <- performance(preds[[s]], measure="tnr")@y.values[[1]]
  fdr_values <- performance(preds[[s]], measure="pcfall")@y.values[[1]]
  OPM_values <- NULL
  
  for (i in 1:length(ROC_curves[[s]]@alpha.values[[1]])) {
    nMCC <- (1+MCC_values[i]) / 2
    OPM_values <- c(OPM_values, ((PPV_values[i]+NPV_values[i])*(tpr_values[i]+tnr_values[i])*(acc_values[i]+nMCC)) / 8)
  }
  
  results[[s]] <- data.frame(Cutoff=ROC_curves[[s]]@alpha.values[[1]],OPM=OPM_values,ACC=acc_values,F_value=F_values,MCC=MCC_values,TPR=tpr_values,TNR=tnr_values,PPV=PPV_values,NPV=NPV_values,FDR=fdr_values, stringsAsFactors = F)
  results[[s]] <- results[[s]][!is.infinite(results[[s]]$Cutoff),]
}

ROC_tidy <- data.frame(Score=character(), x=numeric(), y=numeric(), AUC=numeric(), stringsAsFactors = F)
for (n in names(ROC_curves)) {
  df <- data.frame(x=ROC_curves[[n]]@x.values[[1]], y=ROC_curves[[n]]@y.values[[1]], AUC=round(AUC_values[[n]],2), stringsAsFactors = F)
  df$Score <- paste0(n," (AUC ",round(AUC_values[[n]],2),")")
  ROC_tidy <- rbind(ROC_tidy,df)
}

PrecRecall_tidy <- data.frame(Score=character(), x=numeric(), y=numeric(), stringsAsFactors = F)
for (n in names(ROC_curves)) {
  df <- data.frame(x=PrecRecall_curves[[n]]@x.values[[1]], y=PrecRecall_curves[[n]]@y.values[[1]], stringsAsFactors = F)
  df$Score <- paste0(n," (AUC ",round(AUC_values[[n]],2),")")
  PrecRecall_tidy <- rbind(PrecRecall_tidy,df)
}
PrecRecall_tidy$Score <- factor(PrecRecall_tidy$Score, 
                                levels=levels(reorder(ROC_tidy$Score, -(ROC_tidy$AUC))))

# Optimal thresholds based on ACC, TPR and FDR ----------------------
Optimal_cutoff <- list()
tpr90_cutoff <- list()
fdr50_cutoff <- list()
for (s in scores) {
  idx <- which.max(results[[s]]$ACC)
  Optimal_cutoff[[s]] <- results[[s]][idx,]
  results[[s]] <- results[[s]][order(results[[s]]$TPR, decreasing = T),]
  idx <- sum(results[[s]]$TPR >= 0.9)
  tpr90_cutoff[[s]] <- results[[s]][idx,]
  results[[s]] <- results[[s]][order(results[[s]]$FDR,results[[s]]$TPR, decreasing=T),]
  idx <- sum(results[[s]]$FDR > 0.5, na.rm = T)
  if (idx < nrow(results[[s]])) { idx <- idx + 1}
  fdr50_cutoff[[s]] <- results[[s]][idx,]
}

# PLOTS -------------------
## Performance comparison plot -------------------
# ROC curve plot
panel_a <- ggplot(ROC_tidy, aes(x=x, y=y, color=reorder(Score,-AUC))) + 
  geom_line() + labs(x="FPR", y="TPR", color="Score") + theme_bw() + 
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.key.size = unit(5,"mm"),)

# Precision-Recall curve plot
panel_b <- ggplot(PrecRecall_tidy %>% filter(Score != "FATHMM_MKL_COD"), aes(x=x, y=y, color=Score)) + 
  geom_line() + labs(x="Recall", y="Precision") + theme_bw() + theme(legend.position = "none")

figure <- panel_a / panel_b
figure <- figure + plot_annotation(tag_levels = 'A') + plot_layout(guides="collect")
  ggsave(figure, filename = "plots/Supp Figure XX - Score performances.png", device="png", width=8, height = 8, dpi=150)

## Cutoff plot ----------------
#Sensitivity / specificity plot with cutoffs
tpr_tnr_tidy <- data.frame(Score=character(), CutOff=numeric(), Thr=numeric(), TPR=numeric(), TNR=numeric(), stringsAsFactors = F)
tpr_tnr_points <- data.frame(Score=character(), Thr=character(), x=numeric(), y=numeric(), label=character(), stringsAsFactors = F)
for (s in scores) {
  df <- results[[s]] %>% select(Cutoff,TPR,TNR) %>% mutate(Score=s)
  tpr_tnr_tidy <- rbind(tpr_tnr_tidy,df)

  label <- paste0("TNR=",round(Optimal_cutoff[[s]]$TNR[[1]],2),"\nTPR=",round(Optimal_cutoff[[s]]$TPR[[1]],2))
  tpr_tnr_points[nrow(tpr_tnr_points)+1,] <-c(s,"max ACC",Optimal_cutoff[[s]]$TNR[[1]],Optimal_cutoff[[s]]$TPR[[1]], label)

  label <- paste0("TNR=",round(tpr90_cutoff[[s]]$TNR[[1]],2),"\nTPR=",round(tpr90_cutoff[[s]]$TPR[[1]],2))
  tpr_tnr_points[nrow(tpr_tnr_points)+1,] <- c(s,"TPR 0.9",tpr90_cutoff[[s]]$TNR[[1]],tpr90_cutoff[[s]]$TPR[[1]], label)
  
  if(length(fdr50_cutoff[[s]]$TNR)>0){
    label <- paste0("TNR=",round(fdr50_cutoff[[s]]$TNR[[1]],2),"\nTPR=",round(fdr50_cutoff[[s]]$TPR[[1]],2))
    tpr_tnr_points[nrow(tpr_tnr_points)+1,] <- c(s,"FDR 0.5",fdr50_cutoff[[s]]$TNR[[1]],fdr50_cutoff[[s]]$TPR[[1]], label)
  }
}
tpr_tnr_points$x <- as.numeric(tpr_tnr_points$x)
tpr_tnr_points$y <- as.numeric(tpr_tnr_points$y)

figure <- ggplot(tpr_tnr_tidy, aes(x=TNR, y=TPR)) + 
  geom_line() + 
  geom_point(data=tpr_tnr_points, aes(x=x,y=y,color=Thr)) + 
  geom_label_repel(data=tpr_tnr_points, aes(x=x,y=y,color=Thr,label=label),size=2.5,show.legend = F,label.padding = 0.1) + 
  facet_wrap(~Score,scales="free_x") + scale_color_brewer(palette="Set1") +
  labs(color="Threshold") +
  theme_bw()
ggsave(figure, filename = "plots/Supp Figure XX - Sens_spec_with_cutoff.png", height=9, width=12, dpi=150)

# Save result tables ------------------
out_thresholds <- NULL
for (s in scores) {
  df <- Optimal_cutoff[[s]]
  df$Threshold <- "max_ACC"
  df$Score <- s
  out_thresholds <- rbind(out_thresholds,df)
  df <- tpr90_cutoff[[s]]
  df$Threshold <- "tpr90"
  df$Score <- s
  out_thresholds <- rbind(out_thresholds,df)
  df <- fdr50_cutoff[[s]]
  df$Threshold <- "fdr50"
  df$Score <- s
  out_thresholds <- rbind(out_thresholds,df)
}

out_performance <- data.frame(Score=character(), AUC=numeric(), max_F1=numeric(), max_MCC=numeric(), max_ACC=numeric(), max_OPM=numeric(), stringsAsFactors = F)
for (s in scores) {
  out_performance[nrow(out_performance)+1,] <- c(s,max(AUC_values[[s]],na.rm = T),max(results[[s]]$F_value,na.rm = T),max(results[[s]]$MCC,na.rm = T),max(results[[s]]$ACC,na.rm = T),max(results[[s]]$OPM,na.rm = T))
}

write.table(out_thresholds, file = "results/Suggested_thresholds.tsv", sep="\t", row.names = F, quote = F)
write.table(out_performance, file = "results/Scores_performances.tsv", sep="\t", row.names = F, quote = F)
save.image("ROC_curves.RData")

