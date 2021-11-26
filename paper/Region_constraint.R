library(data.table)
library(tidyverse)
library(ggsignif)
library(patchwork)

# auxiliary functions ------------
blom =
  function(x, c = 3/8) {
    r = rank(x)
    N = length(x)
    qnorm((r - c) / (N - 2 * c + 1))
  }
logit =
  function(x) {
    log(x) - log1p(-x)
  }
binarise =
  function(x) {
    if_else(x == 0, 0, 1)
  }

# Plot formats ----------------
small_text <- theme(axis.text.x = element_text(angle=45, hjust=1, size=8), 
                    axis.title.x = element_text(size=8),
                    axis.text.y = element_text(size=8),
                    axis.title.y = element_text(size=8),
                    axis.line = element_line(color="black"), 
                    plot.title = element_text(size=10),
                    legend.text = element_text(size=8),
                    legend.title = element_text(size=8, face="bold"),
                    legend.key.size = unit(3,"mm"))
# gene sets -------------
gene_lists <- list(GO_BP="/well/gel/HICF2/ref/MSigDB/c5.bp.v7.1.symbols.gmt",
                   GO_MF="/well/gel/HICF2/ref/MSigDB/c5.mf.v7.1.symbols.gmt",
                   GO_CC="/well/gel/HICF2/ref/MSigDB/c5.cc.v7.1.symbols.gmt",
                   PATHWAYS="/well/gel/HICF2/ref/MSigDB/c2.cp.v7.1.symbols.gmt",
                   REACTOME="/well/gel/HICF2/ref/MSigDB/c2.cp.reactome.v7.1.symbols.gmt",
                   ESSENTIAL_MOUSE="/well/gel/HICF2/ref/gene_lists/gene_lists/lists/mgi_essential.tsv",
                   ESSENTIAL_CORE="/well/gel/HICF2/ref/gene_lists/gene_lists/lists/core_essentials_hart.tsv",
                   HAPLO_SEVERE="/well/gel/HICF2/ref/gene_lists/gene_lists/lists/haploinsufficiency_severe_curated_2016.tsv",
                   CLINVAR_PATH="/well/gel/HICF2/ref/ClinVar/gene_pathogenic_20190704.list")

gene_groups <- list()
for (n in names(gene_lists)) {
  inputfile <- gene_lists[[n]]
  pathways <- list()
  if (endsWith(inputfile,".gmt")) {
    con  <- file(inputfile, open = "r")
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
      myvector <- (strsplit(oneLine, "\t"))
      pathways[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
    }
    close(con)
  } else {
    pathways[[n]] <- scan(inputfile, what="",sep="\n")
  }
  gene_groups[[n]] <- pathways
}

# Load data ----------------------------
data_folder <- "data/region_constraint/"
genes_df <- fread(cmd="zcat data/regions/genes.tsv.gz", header=T, sep="\t")
genes_df <- unique(genes_df[ , .(regionID,gene_symbol)])
tissues <- fread(cmd="zcat data/regions/tissues.tsv.gz", header=T, sep="\t")
tissues <- unique(tissues[ , .(regionID,cell_or_tissue)])

GC <- fread(paste0(data_folder, "GRCh38_regions_GCcontent.simple.tsv"), header=T)
var_data <- fread(paste0(data_folder, "GRCh38_regions_varcounts.simple.tsv"), header=F)
SegDup <- fread(paste0(data_folder,"GRCh38_regions_SegDups.simple.cov"), header=F)
LCR <- fread(paste0(data_folder,"GRCh38_regions_LCR.simple.cov"), header=F)
exons_overlap <- fread(paste0(data_folder,"GRCh38_GREEN-DB.exonsOverlap.tsv"), header=F) %>% 
  group_by(V4) %>% mutate(tot_overlap=sum(V24),size=V3-V2) %>% select(V4, tot_overlap, size) %>%
  distinct() %>% mutate(exonoverlap = tot_overlap / size)
gnomAD <- fread(paste0(data_folder,"gnomad.v2.1.1.lof_metrics.by_gene.txt"), header=T)

df <- merge(var_data,GC[, c("4_usercol","12_pct_gc","19_seq_len")], by.x="V4", by.y="4_usercol")
colnames(df) <- c("regionID","chrom","start","end","std_type","DB_source","N_var","gc_pct","seq_len")
df <- merge(df, SegDup[,c("V4","V7")], by.x="regionID",by.y="V4")
colnames(df)[10] <- "SegDup"
df <- merge(df, LCR[,c("V4","V7")], by.x="regionID",by.y="V4")
colnames(df)[11] <- "LCR"
df <- merge(df, exons_overlap[,c("V4","exonoverlap")], by.x="regionID",by.y="V4", all.x=T) %>%
  replace_na(list(exonoverlap = 0))

df_clean <- df %>% filter(SegDup < 0.5, LCR < 0.5, chrom != "chrY")

# prepare data for the model --------------------
DATA =
  df_clean %>%
  as_tibble() %>%
  select(regionID, N_var, gc_pct, seq_len, SegDup, LCR, exonoverlap)
DATA2 =
  DATA %>%
  mutate(across(-c(regionID, LCR, SegDup, exonoverlap), blom)) %>%
  mutate(LCR = binarise(exonoverlap)) %>%
  mutate(LCR = binarise(LCR)) %>%
  mutate(SegDup = binarise(SegDup))

## Eventually look at data
DATA %>%
  pivot_longer(cols = -regionID, names_to = 'VAR', values_to = 'VALUE') %>%
  ggplot() +
  geom_histogram(aes(x = VALUE), bins = 100) +
  facet_wrap(~VAR, scales = 'free')
DATA %>%
  sample_n(10000) %>%
  pivot_longer(cols = -c(regionID, N_var),
               names_to = 'VAR', values_to = 'VALUE') %>%
  ggplot() +
  geom_point(aes(x = VALUE, y = N_var), size = 0.25, alpha = 0.5) +
  facet_wrap(~VAR, scales = 'free')

DATA2 %>%
  pivot_longer(cols = -regionID, names_to = 'VAR', values_to = 'VALUE') %>%
  ggplot() +
  geom_histogram(aes(x = VALUE), bins = 100) +
  facet_wrap(~VAR, scales = 'free')
DATA2 %>%
  sample_n(10000) %>%
  pivot_longer(cols = -c(regionID, N_var),
               names_to = 'VAR', values_to = 'VALUE') %>%
  ggplot() +
  geom_point(aes(x = VALUE, y = N_var), size = 0.25, alpha = 0.5) +
  facet_wrap(~VAR, scales = 'free')

# linear regression --------------
RES =
  DATA2 %>%               ## <- alternatively, change this to DATA
  as.data.frame() %>%
  column_to_rownames(var = 'regionID') %>%
  lm(N_var ~ ., data = .)

RES_DF = broom::augment(RES)
residuals <- RES_DF[order(RES_DF$.resid, decreasing=T),]
residuals$order <- seq(1:nrow(residuals))
residuals$pct <- residuals$order / nrow(residuals)

#Look at residuals
summary(RES)
RES_DF %>%
  ggplot() +
  geom_histogram(aes(x = .resid), bins = 100, color = 'black', fill = 'white')

# Find constrained regions -------------------------------
df_clean <- merge(df_clean,residuals[,c(".rownames",".resid","pct")],by.x="regionID",by.y=".rownames", all.x=T) %>%
  arrange(desc(pct)) %>% 
  mutate(group = case_when(pct >= 0.99 ~ "constrained_regions", TRUE ~ "other_regions"))
df_clean$percentiles <- cut(df_clean$pct, 
                            breaks = c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,.99,1), 
                            labels = c("<0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-0.99",">0.99"), 
                            include.lowest = T)

constrained_NCR <- df_clean %>% filter(pct >= 0.99) %>% pull(regionID)
#ggplot(aes(x=group, y=oe_lof)) + geom_violin() + geom_boxplot(width=0.2)


# Correlation with PhyloP conservation ---------------
phyloP <- fread("data/regions/GRCh38_regions_scores.phyloP_pct.fix.tsv", header=T, sep="\t")
phyloP <- phyloP[ , .(regionID,pct_P100_min_2)][, group := fcase(
  regionID %in% constrained_NCR, "constrained_regions",
  default = "other_regions"
)]

# GSEA for genes connected to constrained regions -------------------------------
## set genes of intereset and background ---------------------
### Genes of interest are those controlled by constrained regions
genes <- unique(genes_df$gene_symbol[genes_df$regionID %in% constrained_NCR]) #4579

### Load background genes (all controlled genes in the database)
backgenes <- unique(genes_df$gene_symbol)

## Compute hypergeometric test enrichment -----------------
gene_groups_enrich <- list()
for (n in names(gene_groups)) {
  message("Enrichment for ", n)
  results <- data.frame(category=character(), pvalue=numeric(), results_in_path=numeric(), backgenes_in_path=numeric(), tot_results=numeric(), tot_backgenes=numeric(), stringsAsFactors=FALSE)
  for (h in 1:length(gene_groups[[n]])) {
    successes <- length(which(genes %in% gene_groups[[n]][[h]]))
    tot_path <- length(which(backgenes %in% gene_groups[[n]][[h]]))
    not_path <- length(backgenes) - tot_path
    results[nrow(results)+1,] <- c(names(gene_groups[[n]][h]), phyper(successes-1, tot_path, not_path, length(genes), lower.tail=FALSE), successes, tot_path, length(genes), length(backgenes) )
  }
  results$backgenes_in_path <- as.numeric(results$backgenes_in_path)
  results$results_in_path <- as.numeric(results$results_in_path)
  results$pvalue <- as.numeric(results$pvalue)
  results$BH <- p.adjust(results$pvalue, method = "BH")
  results$BY <- p.adjust(results$pvalue, method = "BY")
  results <- results[order(results$pvalue),]
  gene_groups_enrich[[n]] <- results
}

## Merge results across all groups and re-compute FDR -----------------
GSEA_constrained_NCR <- NULL
for (n in names(gene_groups_enrich)) {
  df <- gene_groups_enrich[[n]]
  df$group <- n
  GSEA_constrained_NCR <- rbind(GSEA_constrained_NCR, df)
}
GSEA_constrained_NCR <- GSEA_constrained_NCR %>% 
  select(group,category,pvalue,results_in_path) %>%
  mutate(BH = p.adjust(pvalue, method="BH"),
         BY = p.adjust(pvalue, method="BY")) %>%
  arrange(BY)
genes_constrained_NCR <- genes

## Save to table -------------------
write.table(GSEA_constrained_NCR, file="results/GSEA_constrained_NCR.tsv", sep="\t", row.names=F, quote=F)

# Analysis of N_tissues for constrained regions ----------------
tissues_counts <- tissues %>% group_by(regionID) %>% summarise(N_tissues=n()) %>% 
  mutate(group = case_when(regionID %in% constrained_NCR ~ "99th_pct_constrained",
                           TRUE ~ "Other regions"))

# Analysis of N_controlled_genes for constrained regions ------------
genes_per_region <- genes_df %>% group_by(regionID) %>% summarise(N_genes=n()) %>% 
  mutate(group = case_when(regionID %in% constrained_NCR ~ "99th_pct_constrained",
                           TRUE ~ "Other regions"))

# Chi square test for tissue / gene specific constrained regions ----------
constrained_1_tissue <- length(tissues_counts$regionID[tissues_counts$group=="99th_pct_constrained" & tissues_counts$N_tissues == 1])
constrained_multi_tissue <- length(tissues_counts$regionID[tissues_counts$group=="99th_pct_constrained"]) - constrained_1_tissue
other_1_tissue <- length(tissues_counts$regionID[tissues_counts$group=="Other regions" & tissues_counts$N_tissues == 1])
other_multi_tissue <- length(tissues_counts$regionID[tissues_counts$group=="Other regions"]) - other_1_tissue
mymat <- matrix(c(constrained_1_tissue, constrained_multi_tissue,other_1_tissue,other_multi_tissue), ncol=2)
chisq.tissue <- chisq.test(mymat)

constrained_1_gene <- length(genes_per_region$regionID[genes_per_region$group=="99th_pct_constrained" & genes_per_region$N_genes == 1])
constrained_multi_gene <- length(genes_per_region$regionID[genes_per_region$group=="99th_pct_constrained"]) - constrained_1_gene
other_1_gene <- length(genes_per_region$regionID[genes_per_region$group=="Other regions" & genes_per_region$N_genes == 1])
other_multi_gene <- length(genes_per_region$regionID[genes_per_region$group=="Other regions"]) - other_1_gene
mymat <- matrix(c(constrained_1_gene, constrained_multi_gene,other_1_gene,other_multi_gene), ncol=2)
chisq.gene <- chisq.test(mymat)

# Analysis of constraint distribution comparing normal and clinvar / essential genes ----------
genes_df[, clinvar := fcase(gene_symbol %in% gene_groups$CLINVAR_PATH$CLINVAR_PATH, "YES", default = "NO")]
genes_df[, essential := fcase(gene_symbol %in% gene_groups$ESSENTIAL_MOUSE$ESSENTIAL_MOUSE, "YES", default = "NO")]
genes_df_pct <- merge(genes_df, df_clean[,.(regionID,pct)], by="regionID")
genes_df_pct$percentiles <- cut(genes_df_pct$pct, breaks = c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,.99,1), labels = c("<0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-0.99",">0.99"), include.lowest = T)
genes_df_maxpct_bygene_gnomAD <- merge(genes_df_maxpct_bygene, gnomAD[,c("gene", "pLI", "oe_lof", "oe_mis")], by.x="gene_symbol", by.y="gene") 

genes_df_maxpct_bygene <- genes_df_pct %>%
  group_by(gene_symbol) %>%  
  mutate(pct=max(pct)) %>% 
  select(gene_symbol,pct,clinvar,essential) %>% 
  mutate(percentiles = cut(pct, 
                           breaks = c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,.99,1), 
                           labels = c("<0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-0.99",">0.99"), 
                           include.lowest = T)) %>%
  distinct()

maxpct_bygene_counts <- list()
maxpct_bygene_counts$clinvar <- genes_df_maxpct_bygene %>% 
  group_by(percentiles, clinvar) %>% 
  summarise(count=n()) %>% 
  group_by(clinvar) %>% 
  mutate(tot_genes=sum(count), pct_genes=count / tot_genes)
maxpct_bygene_counts$essential <- genes_df_maxpct_bygene %>% 
  group_by(percentiles, essential) %>% 
  summarise(count=n()) %>% 
  group_by(essential) %>% 
  mutate(tot_genes=sum(count), pct_genes=count / tot_genes)

genes_df_medianpct_bygene <- genes_df_pct %>%
  group_by(gene_symbol) %>%  
  mutate(pct=median(pct)) %>% 
  select(gene_symbol,pct,clinvar,essential) %>% 
  mutate(percentiles = cut(pct, 
                           breaks = c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,.99,1), 
                           labels = c("<0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-0.99",">0.99"), 
                           include.lowest = T)) %>%
  distinct()

medianpct_bygene_counts <- list()
medianpct_bygene_counts$clinvar <- genes_df_medianpct_bygene %>% 
  group_by(percentiles, clinvar) %>% 
  summarise(count=n()) %>% 
  group_by(clinvar) %>% 
  mutate(tot_genes=sum(count), pct_genes=count / tot_genes)
medianpct_bygene_counts$essential <- genes_df_medianpct_bygene %>% 
  group_by(percentiles, essential) %>% 
  summarise(count=n()) %>% 
  group_by(essential) %>% 
  mutate(tot_genes=sum(count), pct_genes=count / tot_genes)

# Enrichment of TrueSet vars in constrained regions ----------
TN_vars <- fread(paste0(data_folder, "TrueSet_TN_vars.tsv"), header=F)
TP_vars <- fread(paste0(data_folder, "TrueSet_TP_vars.tsv"), header=F)
TN_vars[, vid := paste(V1,V2,V4,V5,sep="_")]
TP_vars[, vid := paste(V1,V2,V4,V5,sep="_")]
v1 <- length(unique(TP_vars[V9 %in% constrained_NCR, vid]))
v2 <- length(unique(TN_vars[V9 %in% constrained_NCR, vid]))
v3 <- length(unique(TP_vars[V9 %nin% constrained_NCR & V9 != ".", vid]))
v4 <- length(unique(TN_vars[V9 %nin% constrained_NCR & V9 != ".", vid]))
mat <- matrix(c(v1,v2,v3,v4), ncol=2)
truevars_fisher <- fisher.test(mat)

# Plots -------------
## Model scatter plot and residuals histogram --------------
panel_a <- ggplot(df_clean,aes(x=seq_len,y=N_var)) + 
  geom_point(aes(color=percentiles),shape=".",size=1) + 
  labs(x="Sequence length (bp)", y="N variants") +
  theme_bw() + scale_color_viridis_d()

panel_b <- RES_DF %>%
  ggplot() +
  geom_histogram(aes(x = .resid), bins = 100, color = 'black', fill = 'white') +
  geom_vline(xintercept = quantile(RES_DF$.resid, 0.01), linetype="dashed") +
  theme_bw() + labs(x="residual")

panel_c <- ggplot(phyloP, aes(x=pct_P100_min_2, fill=group)) + geom_density(alpha=0.5) + 
  scale_y_sqrt() + labs(x="pct bases PhyloP100 >= 2 per region") + 
  theme_bw() + theme(legend.position = "bottom", legend.key.size = unit(3,"mm"))

panel_d <- ggplot(df_clean %>% filter(group=="constrained_regions"), aes(x=std_type)) + 
  geom_bar() + scale_y_log10() + 
  theme_bw() + labs(x="region type")

figure <- (panel_a | panel_b) / (panel_c | panel_d)
figure <- figure + plot_annotation(tag_levels = 'A')
#ggsave(figure, filename = "plots/Supp Figure XX - Regression model constrained.pdf", device="pdf", height=8, width=8)
ggsave(figure, filename = "plots/Supp Figure XX - Regression model constrained.png", device="png", height=8, width=8)


## Fraction gene and tissue specific --------------------
tissue_specific <- data.frame(group=c("constrained_regions","other_regions"), 
                              N_1_tissue=c(constrained_1_tissue,other_1_tissue),
                              N_multiple_tissue=c(constrained_multi_tissue, other_multi_tissue),
                              tot_regs = c((constrained_1_tissue + constrained_multi_tissue), (other_1_tissue + other_multi_tissue))
)
tissue_specific$pct_1tissue <- tissue_specific$N_1_tissue / tissue_specific$tot_regs
panel_a <- ggplot(tissue_specific, aes(x=group, y=pct_1tissue)) + geom_bar(stat="identity") + 
  geom_signif(comparisons = list(c("constrained_regions", "other_regions")), 
              annotations = c("p < 2.2e-16"), size=0.7, textsize = 5, vjust = -0.2, margin_top = 0.2) +
  theme_bw() + lims(y=c(0,0.8)) + labs(y="Pct regions", title="Regions active in a single tissue")

gene_specific <- data.frame(group=c("constrained_regions","other_regions"), 
                              N_1_gene=c(constrained_1_gene,other_1_gene),
                              N_multiple_gene=c(constrained_multi_gene, other_multi_gene),
                              tot_regs = c((constrained_1_gene + constrained_multi_gene), (other_1_gene + other_multi_gene))
)
gene_specific$pct_1gene <- gene_specific$N_1_gene / gene_specific$tot_regs
panel_b <- ggplot(gene_specific, aes(x=group, y=pct_1gene)) + geom_bar(stat="identity") + 
  geom_signif(comparisons = list(c("constrained_regions", "other_regions")), 
              annotations = c("p < 2.2e-16"), size=0.7, textsize = 5, vjust = -0.2, margin_top = 0.2) +
  theme_bw() + lims(y=c(0,0.8)) + labs(y="Pct regions", title="Regions controlling a single gene") 

figure <- panel_a | panel_b
figure <- figure + plot_annotation(tag_levels = 'A')
ggsave(figure, filename = "plots/Supp Figure XX - Constrained tissue gene specific.png", device="png", height=4, width=8)

## Max constraint comparison normal vs clinvar / essential genes -----------------

panel_a <- ggplot(maxpct_bygene_counts$clinvar, aes(x=percentiles, y=pct_genes, fill=clinvar)) + 
  geom_bar(stat="identity", position=position_dodge(0.95)) +
  labs(x="Constraint", y="Fraction of group genes") + 
  theme_bw() + small_text + theme(legend.position = "top", axis.title.x = element_blank())
panel_b <- ggplot(maxpct_bygene_counts$essential, aes(x=percentiles, y=pct_genes, fill=essential)) + 
  geom_bar(stat="identity", position=position_dodge(0.95)) +
  labs(x="Constraint", y="Fraction of group genes") + 
  theme_bw() + small_text + theme(legend.position = "top", axis.title.x = element_blank())
panel_c <- ggplot(genes_df_maxpct_bygene_gnomAD, aes(x=percentiles, y=oe_lof)) + 
  geom_violin(aes(fill=percentiles), scale = "width") + geom_boxplot(width=0.2, outlier.shape = NA) + 
  lims(y=c(0,3)) + scale_fill_discrete() + 
  labs(x="Constraint", y="gnomAD oe_lof") +
  theme_bw() + small_text + theme(legend.position = "none")

figure <- panel_a / panel_b / panel_c
figure <- figure + plot_annotation(tag_levels = "A")
ggsave(figure, filename = "plots/Figure_4_color.pdf", device="pdf", width = 94, height = 180, dpi=350, units ="mm")

## Supplementary constraint comparison normal vs clinvar / essential genes -----------------
panel_a <- ggplot(genes_df_maxpct_bygene, aes(x=essential, y=pct)) + 
  geom_violin(scale="width", aes(fill=essential)) + geom_boxplot(width=0.1, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("NO", "YES")), annotation="***", textsize = 6, vjust = 0.4) +
  scale_y_continuous(expand=c(0,0.2)) +
  labs(x="Essential genes", y="Region constraint\n(max value)") + 
  theme_bw() + small_text + theme(legend.position = "none") 
panel_b <- ggplot(genes_df_maxpct_bygene, aes(x=clinvar, y=pct)) + 
  geom_violin(scale="width", aes(fill=clinvar)) + geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_signif(comparisons = list(c("NO", "YES")), annotation="***", textsize = 6, vjust = 0.4) +
  scale_y_continuous(expand=c(0,0.2)) +
  labs(x="ClinVar pathogenic", y="Region constraint\n(max value)") +
  theme_bw() + small_text + theme(legend.position = "none")
panel_c <- ggplot(medianpct_bygene_counts$clinvar, aes(x=percentiles, y=pct_genes, fill=clinvar)) + 
  geom_bar(stat="identity", position=position_dodge(0.95)) +
  labs(x="Constraint", y="Fraction of group genes") + 
  theme_bw() + small_text + theme(legend.position = "top")
panel_d <- ggplot(medianpct_bygene_counts$essential, aes(x=percentiles, y=pct_genes, fill=essential)) + 
  geom_bar(stat="identity", position=position_dodge(0.95)) +
  labs(x="Constraint", y="Fraction of group genes") + 
  theme_bw() + small_text + theme(legend.position = "top")
figure <- (panel_a | panel_b) / panel_c / panel_d
figure <- figure + plot_annotation(tag_levels = "A")
ggsave(figure, filename = "plots/Supp Figure XX - Additional constraint for pathogenic and essential.png", device="png", width=8, height = 9)

# ADDITIONAL STUFF --------------
## Analysis of N_regions per gene for genes controlled by constrained regions --------------
regs_per_gene <- genes_df %>% group_by(gene_symbol) %>% mutate(n_regs=n()) %>% select(gene_symbol, n_regs) %>% distinct()
regs_per_gene$group <- "Other regions"
regs_per_gene$group[regs_per_gene$gene_symbol %in% genes_constrained_NCR] <- "99th_pct_constrained"

tmp_df1 <- genes_df %>% filter(regionID %in% constrained_NCR) %>% group_by(gene_symbol) %>% mutate(n_regs=n()) %>% select(gene_symbol, n_regs) %>% distinct()
tmp_df2 <- genes_df %>% filter(!(regionID %in% constrained_NCR)) %>% group_by(gene_symbol) %>% mutate(n_regs=n()) %>% select(gene_symbol, n_regs) %>% distinct()
regs_per_gene <- merge(tmp_df1, tmp_df2, by="gene_symbol", all.x=T, all.y=T)
regs_per_gene[is.na(regs_per_gene)] <- 0
colnames(regs_per_gene) <- c("gene","constrained_regs","other_regs")
regs_per_gene$tot_regs <- regs_per_gene$constrained_regs + regs_per_gene$other_regs
regs_per_gene$pct_constrained <- regs_per_gene$constrained_regs / regs_per_gene$tot_regs

## Overlap between genes controlled by constrained NCR and genes with large reg space ------------
genes_99th_regspace <- regs_per_gene$gene[regs_per_gene$tot_regs>= 182]
length(intersect(genes_99th_regspace, genes_constrained_NCR))

## Analysis of genes with large number of associated constrained NCR --------------
regs_per_gene_constrained_NCR <- associated_genes %>% filter(regionID %in% constrained_NCR) %>% group_by(affected_genes) %>% mutate(n_regs=n()) %>% select(affected_genes, n_regs) %>% distinct()
perc99 <- quantile(regs_per_gene_constrained_NCR$n_regs, c(.9,.95,.99))[[3]]
genes_perc99 <- regs_per_gene_constrained_NCR$affected_genes[regs_per_gene_constrained_NCR$n_regs>= perc99]

GSEA_genes_large_constrained_NCR <- gene_groups_enrich_merged