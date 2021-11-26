library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Variables -----------------
data_folder <- "data/spike_in_vars/"
gado_data <- "data/GADO/"
greendb <- "/well/gel/HICF2/HICF2_hg38_remap/annotations/v1/RegulatoryRegions/GRCh38_regions.annotated.tsv.gz"

x_axis_text_45 <- theme(axis.text.x = element_text(angle=45, hjust=1))
small_text <- theme(axis.text.x = element_text(angle=45, hjust=1, size=8), 
                    axis.title.x = element_text(size=8),
                    axis.text.y = element_text(size=8),
                    axis.title.y = element_text(size=8),
                    axis.line = element_line(color="black"),
                    strip.text = element_text(size=8),
                    legend.text = element_text(size=8),
                    legend.title = element_text(size=8, face="bold"),
                    legend.key.size = unit(3,"mm"))
chr_order <- paste("chr",c(1:22,"X","Y","M"), sep = "")
NC_tranches <- list(
  level1 = dplyr::expr(inGREENDB == 1),
  level2 = dplyr::expr(TFBS == 1 | DNase == 1 | UCNE == 1),
  level3 = dplyr::expr(FATHMM_MKLNC_score >= 0.9 | ReMM_score >= 0.96 | ncER_score >= 96),
  level4 = dplyr::expr(Reg_constraint >= 0.70)
)

# Load data -------------------
spikein_vars <- fread(paste0(data_folder, "GRCh38_genomizer_49_var.anno.csv"), sep=" ", header=T, na.strings = c("NA","."))
colnames(spikein_vars) <- gsub("#*\\[[0-9]+\\]","",colnames(spikein_vars))
spikein_vars <- spikein_vars[, vid := paste(CHROM,POS,REF,ALT, sep="_")] %>% replace_na(list(TFBS=0,DNase=0,UCNE=0,max_AF_all=0))
spikein_vars <- spikein_vars[max_AF_all < 0.01]
greendb_df <- fread(cmd=paste0("zcat ",greendb), header=T, sep="\t")
greendb_df <- as.data.table(greendb_df %>% separate_rows(controlled_genes,sep=","))
spikein_vars <- as.data.table(spikein_vars %>% separate_rows(Reg_id, sep=","))
spikein_vars <- merge(spikein_vars, greendb_df[,.(regionID,controlled_genes)], by.x="Reg_id",by.y="regionID", all.x=T)
spikein_vars <- merge(spikein_vars, greendb_df[,.(regionID,closestProt_symbol)], by.x="Reg_id",by.y="regionID", all.x=T)

# Evaluate if Spike-in variants are captured -----------------
# 45 rare (AF < 0.01) used for evaluation
spikein_vars <- spikein_vars[, samegene := ifelse(disease_gene == controlled_genes | disease_gene == closestProt_symbol, 1, 0)][, samegene := ifelse(all(is.na(samegene)),0,max(samegene, na.rm=T)), by=vid]
spikein_vars <- spikein_vars[, inGREENDB := ifelse(Reg_id == '.', 0, 1)]
spikein_vars_simple <- as.data.table(spikein_vars %>%
  select(vid,disease_gene,ID,inGREENDB,samegene,TFBS,DNase,UCNE,ReMM_score,ncER_score,FATHMM_MKLNC_score,Reg_constraint) %>%
  distinct()) %>% mutate(ID = paste0("O",gsub("_V[0-9]+","",ID)))
tot_vars <- nrow(spikein_vars_simple)
spikein_counts <- data.frame(group=character(),n=numeric(), pct=numeric())
df <- spikein_vars_simple
for (l in names(NC_tranches)) {
  df <- df[eval(NC_tranches[[l]])]
  n <- nrow(df)
  spikein_counts[nrow(spikein_counts)+1,] <- c(l,n,n/tot_vars)
}
spikein_counts$n <- as.numeric(spikein_counts$n)
spikein_counts$pct <- as.numeric(spikein_counts$pct)


## GREENDB regions associated to spike in genes -------------
df1 <- greendb_df[controlled_genes %in% spikein_vars_simple$disease_gene] %>% 
  select(regionID, controlled_genes,stop,start) %>% setnames("controlled_genes","gene")
df2 <- greendb_df[closestProt_symbol %in% spikein_vars_simple$disease_gene] %>% 
  select(regionID, closestProt_symbol,stop,start) %>% setnames("closestProt_symbol","gene")
regions_per_spikein_gene <- rbind(df1,df2) %>% mutate(size=stop-start) %>% select(regionID, size, gene) %>%
  distinct() %>% group_by(gene) %>% summarize(n_regions=n(), bp=sum(size))
genes_dist <- greendb_df %>% filter(controlled_genes != "") %>% 
  select(regionID,stop,start,controlled_genes) %>% 
  mutate(size=stop-start) %>% select(regionID, size, controlled_genes) %>%
  group_by(controlled_genes) %>%
  summarize(n_regions=n(), bp=sum(size)) %>% rename(gene = controlled_genes)
genes_dist$group <- "all_genes"
regions_per_spikein_gene$group <- "spikein_genes"
regions_per_spikein_gene <- rbind(genes_dist, regions_per_spikein_gene)

## Prioritization in NA12878 assuming diseases gene is known --------------------
#For each variant in the 44 validation set, we applied the 4 prioritization levels described above
#At each level, we tested if the variant is captured and how many other non-coding variants were selected affecting the same gene based.

#Based on these counts we computed:
#- TPR (N validation vars captured / total validation vars)
#- PPV (N validation vars captured / total non-coding vars captured)
#- FDR 

## Load NA12878_vars ---------------
NA12878_vars <- fread(paste0(data_folder, "NA12878.anno.rare.GQ10.csv"), sep=" ", header=T, na.strings = c("NA","."))
colnames(NA12878_vars) <- gsub("#*\\[[0-9]+\\]","",colnames(NA12878_vars))
setnames(NA12878_vars, c("NA12878:GT","NA12878:GQ"), c("GT","GQ"))
NA12878_vars <- NA12878_vars[, vid := paste(CHROM,POS,REF,ALT, sep="_")] %>% replace_na(list(TFBS=0,DNase=0,UCNE=0,max_AF_all=0))
NA12878_vars <- NA12878_vars[max_AF_all < 0.01]
NA12878_vars <- as.data.table(NA12878_vars %>% separate_rows(Reg_id, sep=","))
NA12878_vars <- merge(NA12878_vars, greendb_df[,.(regionID,controlled_genes)], by.x="Reg_id",by.y="regionID", all.x=T)
NA12878_vars <- merge(NA12878_vars, greendb_df[,.(regionID,closestProt_symbol)], by.x="Reg_id",by.y="regionID", all.x=T)
NA12878_vars <- NA12878_vars[, inGREENDB := ifelse(Reg_id == '.', 0, 1)]
NA12878_vars_hom <- NA12878_vars[GT == "1/1"]
NA12878_vars_het <- NA12878_vars[GT == "0/1" | GT == "1/2"]

## Count prioritized vars in known_disease_genes ---------------
df <- copy(spikein_vars_simple)
df <- df[, is_validation := 1]
col_names <- c("vid","inGREENDB","TFBS","DNase","UCNE","ReMM_score","ncER_score","FATHMM_MKLNC_score","Reg_constraint","is_validation")
hom_df <- NA12878_vars_hom[, vid := paste(CHROM,POS,REF,ALT, sep="_")][, is_validation := 0]
het_df <- NA12878_vars_het[, vid := paste(CHROM,POS,REF,ALT, sep="_")][, is_validation := 0]
prioritized_var_counts <- data.frame(level=character(),GT=character(), n=numeric(),validation_captured=character())
for (i in 1:nrow(df)) {
  hom_vars <- hom_df[controlled_genes == df$disease_gene[i] | closestProt_symbol == df$disease_gene[i]]
  het_vars <- het_df[controlled_genes == df$disease_gene[i] | closestProt_symbol == df$disease_gene[i]]
  hom_vars <- rbind(hom_vars[,..col_names], df[i,..col_names])
  het_vars <- rbind(het_vars[,..col_names], df[i,..col_names])
  for (l in names(NC_tranches)) {
  hom_vars <- hom_vars[eval(NC_tranches[[l]])]
  het_vars <- het_vars[eval(NC_tranches[[l]])]
  if (1 %in% hom_vars$is_validation) {hom_known_captured <- "Y"} else {hom_known_captured <- "N"}
  if (1 %in% het_vars$is_validation) {het_known_captured <- "Y"} else {het_known_captured <- "N"}
  n_hom <- length(unique(hom_vars[, vid]))
  n_het <- length(unique(het_vars[, vid]))
  prioritized_var_counts[nrow(prioritized_var_counts)+1,] <- c(l,"hom",n_hom,hom_known_captured)
  prioritized_var_counts[nrow(prioritized_var_counts)+1,] <- c(l,"het",n_het,het_known_captured)
  }
}
prioritized_var_counts$n <- as.numeric(prioritized_var_counts$n)

# Prioritization in NA12878 based on HPO profiles -------------------
#1. Across the 29 OMIM conditions associated with validation vars, we used HPOs from HPO official repository to perform GADO predictions
#2. For diseases associated with more than 5 HPOs, 5 terms were randomly selected
#3. We performed gene ranking predictions using GADO CLI 1.0.1 and July 2020 prediction matrix
#4. All diseases genes associated to validation variants were present in the genes evaluated by GADO, leaving a total of 45 variants
#5. For each variant in this set, we applied the 4 prioritization levels described above
#6. At each level, we tested if the variant is captured and how many other non-coding variants were selected overall or considering only genes above the 90th, 95th and 99th percentile in GADO score (strong predicted association with phenotype)

## Read HPOs and create GADO input profiles ----------------
HPOs <- read.csv(paste0(gado_data, "HPO_term_for_diseases.txt"), sep="\t", header=F)
HPOs$V3 <- gsub(":","",HPOs$V3)
HPO_per_disease <- HPOs %>% group_by(V3) %>% summarize(N_hpo=n()) 

#The following is code to generate simulated HPO profiles
#When more than 5 HPOs are present for a specific disease, it randomly subsamples to 5
#for (d in HPO_per_disease$V3) {
#  if (length(HPOs$V1[HPOs$V3 == d]) > 5) {
#    rnd_HPOs <- sample(HPOs$V1[HPOs$V3 == d],size=5,replace = F)
#  } else {
#    rnd_HPOs <- HPOs$V1[HPOs$V3 == d]
#  }
#  write(paste(c(d, rnd_HPOs), collapse="\t"), file=paste0(gado_data, "GADO_input.txt"), append=T)
#}


## Read GADO_prioritization results ---------------
files <- list.files(paste0(gado_data, "GADO_predictions"), pattern = "txt")
GADO_predictions <- list()
for (f in files) {
  mim <- gsub("\\.txt","", f)
  GADO_predictions[[mim]] <- fread(paste0(gado_data, "GADO_predictions/",f), header=T, sep="\t")
  GADO_predictions[[mim]] <- GADO_predictions[[mim]][, pct := 1-(Rank/.N)]
}
GADO_genes <- GADO_predictions[[1]]$Hgnc

## GADO ranking for the causative genes ------------------
df <- spikein_vars_simple[disease_gene %in% GADO_genes & ID %in% names(GADO_predictions)] #45 vars
HPO_plot_df <- HPO_per_disease %>% mutate(disease_gene=0, GADO_rank =0, GADO_percentile=0)
for (i in 1:nrow(df)) {
  HPO_plot_df$disease_gene[HPO_plot_df$V3 == df$ID[i]] = df$disease_gene[i]
  HPO_plot_df$GADO_rank[HPO_plot_df$V3 == df$ID[i]]= GADO_predictions[[df$ID[i]]]$Rank[GADO_predictions[[df$ID[i]]]$Hgnc == df$disease_gene[i]]
  HPO_plot_df$GADO_percentile[HPO_plot_df$V3 == df$ID[i]] = GADO_predictions[[df$ID[i]]]$pct[GADO_predictions[[df$ID[i]]]$Hgnc == df$disease_gene[i]]
}
HPO_plot_df <- HPO_plot_df %>% filter(disease_gene != "0") %>% 
  select(V3,disease_gene,N_hpo, GADO_rank, GADO_percentile) %>%
  gather(key = "metric", value="value", N_hpo:GADO_percentile) %>% 
  mutate(disease_gene = paste0("(",disease_gene,")"), V3 = paste(V3,disease_gene,sep=""))
HPO_plot_df$metric <- factor(HPO_plot_df$metric, levels=c("N_hpo","GADO_rank","GADO_percentile"))
order_V3 <- HPO_plot_df %>% filter(metric == "N_hpo") %>% arrange(desc(value)) %>% pull(V3)
HPO_plot_df$V3 <- factor(HPO_plot_df$V3, levels=order_V3)

## Candidate vars using GADO prioritization -------------------
df <- copy(spikein_vars_simple)
df <- df[, is_validation := 1][, closestProt_symbol := ""]
setnames(df, "disease_gene","controlled_genes")
col_names <- c("vid","inGREENDB","TFBS","DNase","UCNE","ReMM_score","ncER_score","FATHMM_MKLNC_score","Reg_constraint","is_validation","controlled_genes", "closestProt_symbol")
hom_df <- NA12878_vars_hom[, vid := paste(CHROM,POS,REF,ALT, sep="_")][, is_validation := 0]
het_df <- NA12878_vars_het[, vid := paste(CHROM,POS,REF,ALT, sep="_")][, is_validation := 0]
GADO_prioritized_vars <- data.frame(level=character(),GT=character(), GADO_pct=character(), tot_genes= numeric(), n_genes=numeric(), n_vars=numeric())
for (i in 1:nrow(df)) {
  for (p in c(0,0.9,0.95,0.99)) {
    top_genes <- unique(GADO_predictions[[df$ID[i]]][pct >= p & Hgnc != "", Hgnc])
  
    hom_vars <- hom_df[controlled_genes %in% top_genes | closestProt_symbol %in% top_genes]
    het_vars <- het_df[controlled_genes %in% top_genes | closestProt_symbol %in% top_genes]
    hom_vars <- rbind(hom_vars[,..col_names], df[i,..col_names])
    het_vars <- rbind(het_vars[,..col_names], df[i,..col_names])
    for (l in names(NC_tranches)) {
      hom_vars <- hom_vars[eval(NC_tranches[[l]])]
      het_vars <- het_vars[eval(NC_tranches[[l]])]
      n_genes_hom <- length(unique(c(hom_vars[,controlled_genes], hom_vars[, closestProt_symbol], df[i,controlled_genes])))
      n_genes_het <- length(unique(c(het_vars[,controlled_genes], het_vars[, closestProt_symbol], df[i,controlled_genes])))
      if (1 %in% hom_vars$is_validation) {hom_known_captured <- "Y"} else {hom_known_captured <- "N"}
      if (1 %in% het_vars$is_validation) {het_known_captured <- "Y"} else {het_known_captured <- "N"}
      n_hom <- length(unique(hom_vars[, vid]))
      n_het <- length(unique(het_vars[, vid]))
      GADO_prioritized_vars[nrow(GADO_prioritized_vars)+1,] <- c(l,"hom",p,length(top_genes), n_genes_hom, n_hom)
      GADO_prioritized_vars[nrow(GADO_prioritized_vars)+1,] <- c(l,"het",p,length(top_genes), n_genes_het, n_het)
    }
  }
}
GADO_prioritized_vars$n_vars <- as.numeric(GADO_prioritized_vars$n_vars)
GADO_prioritized_vars$n_genes <- as.numeric(GADO_prioritized_vars$n_genes)

# PLOTS ----------------

## N regions and bp associated to validations vars
p1 <- ggplot(regions_per_spikein_gene, aes(x=n_regions, fill=group)) + geom_density(alpha=.5)
p2 <- ggplot(regions_per_spikein_gene, aes(x=bp, fill=group)) + geom_density(alpha=.5)
p <- p1 / p2

## Supp figure - Summary of candidate diseases and HPO profiles -----------------
p1 <- ggplot(HPO_plot_df %>% filter(metric=="N_hpo"), aes(x=gsub("OMIM","",V3), y=value)) + geom_bar(stat="identity") + 
  theme_bw() + small_text + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(y="N") +
  scale_y_continuous(expand = c(0.03,0), limits=c(0,50)) +
  facet_wrap(~metric)
p2 <- ggplot(HPO_plot_df %>% filter(metric=="GADO_rank"), aes(x=gsub("OMIM","",V3), y=value)) + geom_bar(stat="identity") + 
  theme_bw() + small_text + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(y="Rank") +
  scale_y_continuous(expand = c(0.03, 0)) +
  facet_wrap(~metric)
p3 <- ggplot(HPO_plot_df %>% filter(metric=="GADO_percentile"), aes(x=gsub("OMIM","",V3), y=value)) + 
  geom_bar(stat="identity") + 
  theme_bw() + small_text + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) + 
  labs(y="Pct", x ="disease OMIM (gene)") +
  scale_y_continuous(expand=c(0.15,0)) + coord_cartesian(ylim=c(0.9,1.0)) +
  facet_wrap(~metric)

p_HPOs <- p1 / p2 / p3 
ggsave(plot = p_HPOs, filename = "plots/Supp Figure XX - HPO profiles.png", height=7, width=7, dpi=150, device = "png")

## Figure 5 - Results for validations vars -----------------------
#Validation vars captured and N candidates when disease gene is known
p2 <- ggplot(prioritized_var_counts,aes(x=n, fill=validation_captured)) + geom_histogram(bins=max(prioritized_var_counts$n)) +
  scale_x_continuous(breaks=c(0:5,7,seq(10,25,5))) + scale_y_sqrt(breaks=c(0,1,10,25,50,100)) + 
  scale_fill_brewer(palette="Set1") + facet_grid(level~GT, scales="free_x") + 
  theme_bw() + small_text + theme(panel.grid.minor = element_blank()) + 
  labs(x="N candidate vars selected", y="N validation examples", fill="true var\ncaptured")
p1 <- ggplot(spikein_counts, aes(x="Validation\nvars", y=pct, label=n)) + 
  geom_bar(stat="identity") + 
  geom_label(x=1, y=0.5, size=2.5, label.padding = unit(1,"mm")) + 
  theme_bw() + small_text + 
  theme(axis.ticks.x = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5, size=8),
        axis.title.x = element_blank()) + 
  facet_wrap(~group, nrow=4) +
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.5)) +
  labs(y="% captured")
p_knowngene <- patchworkGrob(p1 + p2 + plot_layout(widths=c(1,12)))

#GADO prioritized vars
GADO_prioritized_vars$GT <- factor(GADO_prioritized_vars$GT, levels=c("het","hom"), labels =c("dominant","recessive"))
p_GADO <- ggplot(
  GADO_prioritized_vars %>%
    gather(key="variable",value="n",n_genes:n_vars)) +
  geom_pointrange(
    aes(x=as.factor(GADO_pct),color=variable, y=n), 
    position=position_dodge(0.5), 
    stat = "summary", 
    fun.min="min", fun="median", fun.max="max", 
    size=0.6, fatten=1.5) +
  scale_y_log10(limits=c(1,1000000), breaks=c(10,1000,100000)) + scale_color_brewer(palette="Set1") +
  facet_grid(level~GT) +
  theme_bw() + small_text + 
  theme(legend.position="right", panel.grid.minor = element_blank()) + 
  labs(x="GADO percentile threshold", y="N candidates")

#p_GADO <- patchworkGrob(p3 + p4 + plot_layout(widths=c(1,8)))

#Make paper figure 
layout <- "A
          B"
p <- wrap_plots(p_knowngene, p_GADO) + plot_annotation(tag_levels= 'A') + plot_layout(design = layout)
ggsave(plot = p, filename = "plots/Figure_5_color.pdf", height=160, width=190, dpi=350, unit="mm", device = "pdf")
