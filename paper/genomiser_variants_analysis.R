library(dplyr)
library(tidyr)
library(ggplot2)

findRank <- function(x, y) {
  rank_value <- sum(as.numeric(y) >= as.numeric(x)) + 1
  return(rank_value)
}

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

# READ DATA --------------
var_disease_genes <- read.table("data/genomiser_comparison/45_rarevars_disease_genes.tsv", header=F, sep="\t")
colnames(var_disease_genes) <- c("vstring","vid","disease_gene")
var_disease_genes$MIM <- gsub("_V[0-9]+","",var_disease_genes$vid)

results_dir <- "data/genomiser_comparison/output_49var/"
validated_vars_genomiser <- NULL
for (vid in unique(var_disease_genes$vid)) {
  df <- read.csv(paste0(results_dir,vid,".variants.tsv"), sep="\t", header=T)
  validated_vars_genomiser <- rbind(
    validated_vars_genomiser,
    df %>% 
      mutate(vid=vid) %>%
      select(vid,FUNCTIONAL_CLASS,EXOMISER_GENE,EXOMISER_GENE_COMBINED_SCORE,EXOMISER_GENE_VARIANT_SCORE,EXOMISER_GENE_PHENO_SCORE)
  )
}

results_dir <- "data/genomiser_comparison/output_NA12878/"
NA12878_vars_genomiser <- list(recessive=list(), dominant=list())
for (mim in unique(var_disease_genes$MIM)) {
  for (i in c("AR","XR")) {
    df <- read.csv(paste0(results_dir,mim,"_",i,".variants.tsv"), sep="\t", header=T)
    NA12878_vars_genomiser$recessive[[mim]] <- rbind(
      NA12878_vars_genomiser$recessive[[mim]],
      df %>%
        filter(GENOTYPE == "1/1") %>%
        mutate(vid=vid) %>%
        select(vid,EXOMISER_GENE,EXOMISER_GENE_COMBINED_SCORE,EXOMISER_GENE_VARIANT_SCORE,EXOMISER_GENE_PHENO_SCORE)
    )
  }
  for (i in c("AD","XD")) {
    df <- read.csv(paste0(results_dir,mim,"_",i,".variants.tsv"), sep="\t", header=T)
    NA12878_vars_genomiser$dominant[[mim]] <- rbind(
      NA12878_vars_genomiser$dominant[[mim]],
      df %>% 
        filter(FILTER == "PASS") %>%
        mutate(vid=vid) %>%
        select(vid,EXOMISER_GENE,EXOMISER_GENE_COMBINED_SCORE,EXOMISER_GENE_VARIANT_SCORE,EXOMISER_GENE_PHENO_SCORE)
    )
  }
}

# MERGE VAR INFO AND CHECK RIGHT GENE ---------
validated_vars_genomiser <- merge(var_disease_genes, validated_vars_genomiser, by="vid")
validated_vars_genomiser$correct_gene <- ifelse(validated_vars_genomiser$disease_gene == validated_vars_genomiser$EXOMISER_GENE,1,0)

# COMPUTE RANKING --------------
validated_vars_genomiser$recessive_combined_rank <- apply(
  validated_vars_genomiser,1,
  function(x) findRank(
    x["EXOMISER_GENE_COMBINED_SCORE"], 
    NA12878_vars_genomiser$recessive[[x["MIM"]]]$EXOMISER_GENE_COMBINED_SCORE
  )
)
validated_vars_genomiser$recessive_variant_rank <- apply(
  validated_vars_genomiser,1,
  function(x) findRank(
    x["EXOMISER_GENE_VARIANT_SCORE"], 
    NA12878_vars_genomiser$recessive[[x["MIM"]]]$EXOMISER_GENE_VARIANT_SCORE
    )
  )
validated_vars_genomiser$dominant_combined_rank <- apply(
  validated_vars_genomiser,1,
  function(x) findRank(
    x["EXOMISER_GENE_COMBINED_SCORE"], 
    NA12878_vars_genomiser$dominant[[x["MIM"]]]$EXOMISER_GENE_COMBINED_SCORE
  )
)
validated_vars_genomiser$dominant_variant_rank <- apply(
  validated_vars_genomiser,1,
  function(x) findRank(
    x["EXOMISER_GENE_VARIANT_SCORE"], 
    NA12878_vars_genomiser$dominant[[x["MIM"]]]$EXOMISER_GENE_VARIANT_SCORE
  )
)

# COUNT VARS FOR THE KNOWN GENE --------------
validated_vars_genomiser$recessive_vars_for_gene <- apply(
  validated_vars_genomiser,1,
  function(x) nrow(
    NA12878_vars_genomiser$recessive[[x["MIM"]]] %>% filter(EXOMISER_GENE == x["disease_gene"])
  ) + ifelse(x["correct_gene"] == 1, 1, 0)
)
validated_vars_genomiser$dominant_vars_for_gene <- apply(
  validated_vars_genomiser,1,
  function(x) nrow(
    NA12878_vars_genomiser$dominant[[x["MIM"]]] %>% filter(EXOMISER_GENE == x["disease_gene"])
  ) + ifelse(x["correct_gene"] == 1, 1, 0)
)

# ANALYSIS ---------------
## N variants with correct genes captured by GREEN-VARAN Level1 and Level3 vs Genomiser all and Genomiser TOP 10
N_correct_gene <- data.frame(
  method=character(),
  group=character(),
  N=numeric()
)
N_correct_gene <- N_correct_gene %>% 
  add_row(method="Genomiser", group="All_variants", N=sum(validated_vars_genomiser$correct_gene)) %>% 
  add_row(method="Genomiser",group="TOP10", N=sum(validated_vars_genomiser$correct_gene[validated_vars_genomiser$recessive_combined_rank < 10])) %>%
  add_row(method="GREEN-VARAN",group="All_variants", N=40) %>% 
  add_row(method="GREEN-VARAN",group="Level3", N=32)

p_truevars_captured <- ggplot(N_correct_gene, aes(x=group, y=N)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~method, scales="free_x") +
  scale_y_continuous(expand=c(0,1)) +
  theme_bw() + small_text +
  labs(y="N validated variants")

## Reproduce the histogram of N variants selected when disease gene is known adding the Genomiser data
## Here we consider all variants associated by Genomiser to the relevant disease gene since we assume only this gene is investigated
prioritized_var_counts <- readRDS("GREEN-DB_validated_vars_counts.R4.1.RDS")
prioritized_var_counts <- prioritized_var_counts %>% filter(level=="level1") %>% mutate(method="GREEN-VARAN") %>% select(-validation_captured,-level) 
genomiser_var_counts_rec <- data.frame(
  GT="hom",
  n=validated_vars_genomiser$recessive_vars_for_gene,
  method="Genomiser")
genomiser_var_counts_dom <- data.frame(
  GT="het",
  n=validated_vars_genomiser$dominant_vars_for_gene,
  method="Genomiser")
prioritized_var_counts <- rbind(prioritized_var_counts, genomiser_var_counts_rec, genomiser_var_counts_dom)

p_hist_knowngene <- ggplot(prioritized_var_counts,aes(x=n)) + geom_histogram(bins=max(prioritized_var_counts$n)) +
  scale_x_continuous(breaks=c(0:5,7,seq(10,25,5))) + scale_y_sqrt(breaks=c(0,1,10,25,50,100)) + 
  scale_fill_brewer(palette="Set1") + facet_grid(GT~method) + 
  theme_bw() + small_text + theme(panel.grid.minor = element_blank()) + 
  labs(x="N candidate vars selected", y="N validation examples")


## For NA12878 simulation: 
GADO_prioritized_vars <- readRDS("GADO_prioritized_vars.R4.1.RDS")
# N candidates resulting from GREEN-VARAN Levels without GADO VS Genomiser ranking variant only
N_candidates_var_only <- GADO_prioritized_vars %>% filter(GADO_pct == 0) %>% 
  select(level,GT,n_genes) %>%
  mutate(method="GREEN-VARAN") %>%
  mutate(GT=case_when(GT=="het" ~ "dominant", GT=="hom" ~ "recessive"))
genomiser_hom <- data.frame(method="Genomiser", level="variant_score", GT="recessive", n_genes=validated_vars_genomiser$recessive_variant_rank[validated_vars_genomiser$correct_gene==1])
genomiser_het <- data.frame(method="Genomiser", level="variant_score", GT="dominant", n_genes=validated_vars_genomiser$dominant_variant_rank[validated_vars_genomiser$correct_gene==1])
N_candidates_var_only <- rbind(N_candidates_var_only, genomiser_hom, genomiser_het)

p_NA12878_vars <- ggplot(N_candidates_var_only,aes(x=level,y=n_genes)) + geom_jitter(height=0, size=1) +
  facet_grid(GT~method, scales="free_x") + 
  scale_y_sqrt(breaks=c(0,1000,10000,25000,50000,75000,100000)) +
  theme_bw() + small_text + theme(panel.grid.minor = element_blank()) + 
  labs(x="Prioritization class", y="N candidate genes")

# Ranking of candidate vars at Level 3 + GADO VS Genomiser combined score ranking
GADO_rank_plot <- GADO_prioritized_vars %>% filter((level=="level1" & GADO_pct==0) | (level=="level3" & GADO_pct %in% c(0,0.99))) %>% 
  mutate(GADO_pct = paste0("GV_",level,"\nGADO_",GADO_pct)) %>%
  select(GT,GADO_pct,candidate_var_rank) %>%
  mutate(rank_tranche = cut(candidate_var_rank, breaks=c(0,5,10,20,50,max(candidate_var_rank, na.rm = T)), labels=c("top5","top10","top20","top50","above50"))) %>%
  group_by(GT,GADO_pct,rank_tranche) %>% summarise(N=n()) %>% #filter(GADO_pct != "GADO_0") %>%
  mutate(GT=case_when(GT=="het" ~ "dominant", GT=="hom" ~ "recessive"))
genomiser_rank <- validated_vars_genomiser %>% filter(correct_gene == 1) %>%
  select(dominant_combined_rank, recessive_combined_rank) %>%
  gather(key="GT", value="var_rank", dominant_combined_rank:recessive_combined_rank) %>%
  mutate(GT=gsub("_combined_rank","",GT), GADO_pct="Genomiser\nCombined_score") %>%
  mutate(rank_tranche = cut(var_rank, breaks=c(0,5,10,20,50,max(var_rank, na.rm = T)), labels=c("top5","top10","top20","top50","above50"))) %>%
  group_by(GT,GADO_pct,rank_tranche) %>% summarise(N=n())

GADO_rank_plot <- rbind(GADO_rank_plot, genomiser_rank)
GADO_rank_plot$rank_tranche <- factor(GADO_rank_plot$rank_tranche, levels = c("above50","top50","top20","top10","top5"))

p_ranks <- ggplot(GADO_rank_plot %>% filter(!is.na(rank_tranche)) , aes(x=GADO_pct, fill=rank_tranche, y=N)) + geom_bar(stat="identity") + 
  facet_wrap(~GT, scales="free") + 
  theme_bw() + small_text + scale_fill_brewer(palette = "Set1") +
  labs(x="Rank method", y="N variants", fill="ranking")

#Make paper figure 
layout <- "AABBB
          CCDDD"
p <- wrap_plots(p_truevars_captured, p_hist_knowngene, p_NA12878_vars, p_ranks) + plot_annotation(tag_levels= 'A') + plot_layout(design = layout)
ggsave(plot = p, filename = "plots/Supp Figure XX - Genomiser comparison validated vars.png", height=160, width=250, dpi=150, unit="mm", device = "png")


## For distant enhancer vars: N captured by GREEN-VARAN by level VS N captured by Genomiser with/without correct gene
N_distant_vars <- data.frame(
  method=character(),
  correct_gene=character(),
  N=numeric()
)
N_distant_vars <- N_distant_vars %>% 
  add_row(method="Genomiser", correct_gene="N", N=15) %>% 
  add_row(method="GREEN-VARAN\nLevel2",correct_gene="Y", N=18) %>% 
  add_row(method="GREEN-VARAN\nLevel3",correct_gene="Y", N=13) %>% 
  add_row(method="GREEN-VARAN\nLevel4",correct_gene="Y", N=2)

p_distantvars <- ggplot(N_distant_vars , aes(x=method, fill=correct_gene, y=N)) + geom_bar(stat="identity") + 
  theme_bw() + small_text + scale_fill_brewer(palette = "Set1") +
  labs(x="Method", y="N variants", fill="correct gene")
ggsave(plot = p_distantvars, filename = "plots/Supp Figure XX - Genomiser comparison distant vars.png", height=5, width=5,dpi=150, device = "png")


# SAVE DATA -----------
saveRDS(validated_vars_genomiser, file = "genomiser_results_R4.1.RDS")
