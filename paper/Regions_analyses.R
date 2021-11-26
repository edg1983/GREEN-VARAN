library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)


# Load data -------------------
NCregs <- fread(cmd = "zcat data/regions/GRCh38_GREEN-DB.bed.gz", sep="\t", header=T,stringsAsFactors = F)
TSS <- fread("data/regions/GENCODE.hg38.v33.basic.TSS.bed", sep="\t", header=F, stringsAsFactors = F)
genes_positions <- fread("data/regions/GENCODE.hg38.v33.basic.genes.bed", sep="\t", header=F, stringsAsFactors = F) 
colnames(NCregs)[1] <- "chrom"
NCregs$dimension <- NCregs$stop - NCregs$start

# General variables ----------------
chr_order <- c(paste0("chr",seq(1,22)),"chrX","chrY","chrM")
computational <- c("DECRES","ENCODE-HMM", "SegWey")
experimental <- c("GasperiniEtAl2019","FulcoEtAl2019","JungEtAl2019","PangEtAl2020","VISTA")
curated_db <- c("DiseaseEnhancers","ENCODE cCREs", "EPD6", "EnsemblRegBuild","FANTOM5","FOCS","HACER","RefSeqRegBuild")

# Plot formats ----------------
small_text <- theme(axis.text.x = element_text(angle=45, hjust=1, size=8), 
                    axis.title.x = element_text(size=8),
                    axis.text.y = element_text(size=8),
                    axis.title.y = element_text(size=8),
                    axis.line = element_line(color="black"), 
                    legend.text = element_text(size=8),
                    legend.title = element_text(size=8, face="bold"),
                    legend.key.size = unit(3,"mm"))

# PART1. General DB description --------------------
## Count by db_sources -------------------
db_sources <- as.data.table(NCregs[, .(regionID,std_type,DB_source)] %>% separate_rows(DB_source, sep=","))

#Originally we labeled ENCODE cCREs as BENGI, so need to convert
db_sources$DB_source[db_sources$DB_source == "BENGI"] <- "ENCODE cCREs"

db_sources[, category := fcase(
  DB_source %in% computational, "computational",
  DB_source %in% experimental, "experimental",
  DB_source %in% curated_db, "curated_db",
  default="other_source"
)]
db_sources$category <- factor(db_sources$category, levels=c("curated_db", "experimental", "computational"))
db_sources$order <- as.numeric(db_sources$category)

## Pct of regions with associated informations (gene, tissue, phenos) -----------------
#PHENOS
phenos <- fread(cmd="zcat data/regions/phenos.tsv.gz", header=T, sep="\t", stringsAsFactors = F)
phenos <- unique(phenos[regionID != "", .(regionID, phenotype)][, phenotype := 1])
phenos$phenotype <- as.numeric(phenos$phenotype)
phenos <- merge(NCregs[,.(regionID,std_type)], phenos, by="regionID", all.x=T)
pct_associated_phenos <- phenos %>% group_by(std_type) %>%
  summarise(tot_regions=n(), with_pheno=sum(phenotype, na.rm = T))
pct_associated_phenos[nrow(pct_associated_phenos) + 1, ] <- list("all_regions", 
                                                              sum(pct_associated_phenos$tot_regions), 
                                                              sum(pct_associated_phenos$with_pheno))
pct_associated_phenos$pct <- pct_associated_phenos$with_pheno / pct_associated_phenos$tot_regions

#GENES
pct_associated_genes <- NCregs %>% group_by(std_type) %>% 
  summarise(tot_regions=n(), connected_genes=sum(controlled_genes != ""), 
            close_genes=sum(controlled_genes == "" & closestGene_dist <= 10000))
pct_associated_genes$pct_connected <- pct_associated_genes$connected_genes / pct_associated_genes$tot_regions
pct_associated_genes$pct_close <- pct_associated_genes$close_genes / pct_associated_genes$tot_regions
pct_associated_genes[nrow(pct_associated_genes) + 1, ] <- list("all_regions", 
                                                            sum(pct_associated_genes$tot_regions), 
                                                            sum(pct_associated_genes$connected_genes), 
                                                            sum(pct_associated_genes$close_genes), 
                                                            sum(pct_associated_genes$connected_genes) / sum(pct_associated_genes$tot_regions),
                                                            sum(pct_associated_genes$close_genes) / sum(pct_associated_genes$tot_regions))
pct_associated_genes <- gather(pct_associated_genes, key="key", value="pct", pct_connected:pct_close)

#TISSUES
tissues <- fread(cmd='zgrep -v "_int" data/regions/tissues.tsv.gz', sep="\t", header=T,stringsAsFactors = F)
tissues <- tissues[regionID != "", .(regionID, cell_or_tissue)][, cell_or_tissue := 1]
tissues$cell_or_tissue <- as.numeric(tissues$cell_or_tissue)
tissues <- merge(NCregs[,.(regionID,std_type)],tissues, by="regionID", all.x=T)
pct_associated_tissue <- unique(tissues) %>% group_by(std_type) %>% 
  summarise(tot_regions=n(), with_tissue=sum(cell_or_tissue, na.rm = T))
pct_associated_tissue[nrow(pct_associated_tissue) + 1, ] <- list("all_regions", 
                                                              sum(pct_associated_tissue$tot_regions), 
                                                              sum(pct_associated_tissue$with_tissue))
pct_associated_tissue$pct <- pct_associated_tissue$with_tissue / pct_associated_tissue$tot_regions

tissues_per_region <- tissues %>% group_by(regionID) %>% summarise(n_tissues=n())
tissues_per_region$tier <- cut(tissues_per_region$n_tissues, breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,50,100,max(tissues_per_region$n_tissues)), labels = c("1","2","3","4","5","6","7","8","9","10","11-20","21-50","51-100",">100"), include.lowest = T)

## Coverage on various genomic regions --------------
genomic_cov <- fread("data/regions/Regions_genome_cov.csv", header=T, sep="\t", stringsAsFactors = F)
genomic_cov$genomic_region <- factor(genomic_cov$genomic_region, levels=c("genome","CDS","UTR","exons","introns","intergenic"))

## Plots ------------------
## Panel A - Count by DB source
panel_a <- ggplot(db_sources, aes(x=reorder(DB_source,order), fill=category)) +
  geom_bar(position=position_dodge(1)) + scale_y_log10(expand=c(0,0.2)) + 
  labs(x="Source",y="Count (log scale)") + scale_fill_brewer(palette="Set1") + 
  theme_bw() + small_text + theme(legend.position = "top")

## Panel B - Overall regions coverage across genome
panel_b <- ggplot(genomic_cov[genomic_cov$std_type == "All",],aes(x=genomic_region, y=bases/1000000, label=round(coverage,2))) + 
  geom_bar(stat="identity") + geom_label(size=2,label.padding = unit(0.1,"lines")) +
  scale_y_sqrt(breaks=c(20,100,200,400,800,1200,1600)) + labs(x="Genomic region", y="M bases") + 
  theme_bw() + small_text

## Panel C - Regions by standard type
panel_c <- ggplot(NCregs, aes(x=std_type, fill=std_type)) + geom_bar() + 
  scale_y_log10(expand=c(0,0.2)) + scale_fill_brewer(palette="Set1") +
  labs(x="Standard type", y="Count (log scale)") + 
  theme_bw() + small_text + theme(legend.position = "none")

## Panel D - Regions dimension
panel_d <- ggplot(NCregs, aes(x=std_type, y=dimension, fill=std_type)) + geom_violin(scale="width") + 
  scale_y_sqrt(breaks=c(200,1000,2500,5000,10000,15000)) + scale_fill_brewer(palette="Set1") +
  labs(y="Length (bp)", x="Standard type") + 
  theme_bw() + small_text + theme(legend.position = "none")

## Panel E - Pct of regions with associated information (gene, tissue, pheno)
pct_regions_with_info <- as.data.frame(rbind(
  c("gene",pct_associated_genes$pct[pct_associated_genes$std_type == "all_regions" & pct_associated_genes$key == "pct_connected"]),
  c("tissue",pct_associated_tissue$pct[pct_associated_tissue$std_type == "all_regions"]),
  c("phenotype",pct_associated_phenos$pct[pct_associated_phenos$std_type == "all_regions"])
), stringsAsFactors=F)
colnames(pct_regions_with_info) <- c("info","value")

pct_regions_with_info$value <- as.numeric(pct_regions_with_info$value)
panel_e <- ggplot(pct_regions_with_info, aes(x=info, y=round(value,2))) + geom_bar(stat="identity") + 
  scale_y_continuous(expand=c(0,0.02)) +
  theme_bw() + small_text + labs(x="Associated information", y="% of regions")

## Assemble figure
layout <- "
AAB
CDE
"
figure <- panel_a + panel_b + panel_c + panel_d + panel_e
figure <- figure + plot_annotation(tag_levels = 'A') + plot_layout(design = layout)
ggsave(figure, filename = "plots/Figure_1_color.pdf", device="pdf", width=178, height = 130, dpi=350, units ="mm")


# PART2. Region-gene connections -------------------------

#Make a table of reg regions with only associated genes 1 per line
controlled_genes <- as.data.table(NCregs[controlled_genes != ""] %>% separate_rows(controlled_genes, sep=","))

## Determine upstream/downstream/ingene position -------------
controlled_genes_withGenePos <- merge(controlled_genes, genes_positions, by.x="controlled_genes", by.y="V5")[chrom == V1]
controlled_genes_withGenePos[, place := fcase(
  (V4=="+" & stop < V2) | (V4=="-" & start > V3), "UPSTREAM",
  (V4=="+" & start > V3) | (V4=="-" & stop < V2), "DOWNSTREAM",
  default = "INSIDE_GENE"
)]
controlled_genes_withGenePos$place <- factor(controlled_genes_withGenePos$place, levels = c("DOWNSTREAM","UPSTREAM","INSIDE_GENE"))

## Distribution of region-gene TSS distance ---------------

#Merge with TSS positions
controlled_genes_withTSS <- as.data.table(merge(controlled_genes, TSS, by.x="controlled_genes", by.y="V5", allow.cartesian = T))

#Set region-gene distance selecting the closest distance to TSS
TSS_regions_dist <- controlled_genes_withTSS[chrom == V1] %>% 
  mutate(TSS_dist1 = abs(V2-stop), TSS_dist2 = abs(V2-start)) %>%
  gather(key="key", value="dist",TSS_dist1:TSS_dist2) %>%
  group_by(regionID, controlled_genes) %>%
  mutate(min_dist=min(dist)) %>%
  select(regionID,std_type,DB_source,controlled_genes,V3,min_dist) %>%
  distinct()

#Split TSS dist by db_source
#TSS_regions_dist_source <- TSS_regions_dist %>% separate_rows(db_source, sep=",") %>% distinct()

## Analysis of controlled genes vs closest genes -------------
#Check if the any of the closest genes is among controlled genes
closest_in_controlled <- controlled_genes %>% 
  mutate(closest_equal_controlled=case_when(
    closestGene_symbol==controlled_genes ~ 1,
    closestProt_symbol==controlled_genes ~ 1,
    closestGene_TSS_symbol==controlled_genes ~ 1,
    closestProt_TSS_symbol==controlled_genes ~ 1,
    TRUE ~ 0)) %>% 
  group_by(regionID) %>% 
  mutate(N_closest_in_controlled=sum(closest_equal_controlled), N_controlled=n()) %>%
  select(regionID,std_type,DB_source,closestGene_dist,N_closest_in_controlled,N_controlled) %>%
  distinct() %>%
  mutate(closest_in_controlled = case_when(
    N_closest_in_controlled == N_controlled ~ "ONLY",
    N_closest_in_controlled > 0 & N_closest_in_controlled < N_controlled ~ "YES",
    TRUE ~ "NO"
  ))
closest_in_controlled$closest_in_controlled <- factor(closest_in_controlled$closest_in_controlled, levels = c("NO","YES","ONLY"))
#closest_in_controlled$dist_tag <- cut(closest_in_controlled$closestGene_dist, breaks = c(0,10000,max(closest_in_controlled$closestGene_dist)), include.lowest = T)

#Check if regions within a gene actually controls that gene 
ingene_in_controlled <- controlled_genes %>% 
  filter(closestGene_dist == 0 | closestProt_dist == 0) %>%
  mutate(ingene_equal_controlled=case_when(
    closestGene_symbol==controlled_genes ~ 1,
    closestProt_symbol==controlled_genes ~ 1,
    TRUE ~ 0)) %>% 
  group_by(regionID) %>% 
  mutate(N_ingene_in_controlled=sum(ingene_equal_controlled), N_controlled=n()) %>%
  select(regionID,std_type,DB_source,N_ingene_in_controlled,N_controlled) %>%
  distinct() %>%
  mutate(ingene_in_controlled = case_when(
    N_ingene_in_controlled == N_controlled ~ "ONLY",
    N_ingene_in_controlled > 0 & N_ingene_in_controlled < N_controlled ~ "YES",
    TRUE ~ "NO"
  ))
ingene_in_controlled$ingene_in_controlled <- factor(ingene_in_controlled$ingene_in_controlled, levels = c("NO","YES","ONLY"))

#Distribution of closest gene dist vs %of time closest is among associated genes
#closest_in_associated <- closest_in_associated[order(closest_in_associated$closestGene_dist),]
#enhancer_cum <- as.data.frame(closest_in_associated %>% filter(std_type=="enhancer") )
#enhancer_cum <- enhancer_cum[order(enhancer_cum$closestGene_dist),]
#enhancer_cum <- enhancer_cum %>% mutate(cum_closest=cumsum(N_closest_in_associated), cum_regions=seq(1,nrow(enhancer_cum)))
#enhancer_cum$cum_closest_pct <- enhancer_cum$cum_closest/enhancer_cum$cum_regions

## N genes per region and N regions per gene ------------------
genes_per_reg <- controlled_genes  %>% group_by(regionID) %>% mutate(n_genes=n()) %>% select(regionID,std_type, n_genes) %>% distinct()
genes_per_reg$level <- cut(genes_per_reg$n_genes, 
                           breaks = c(0,1,5,10,20,max(genes_per_reg$n_genes)), 
                           labels = c("1","2-5","6-10","11-20",">20"), 
                           include.lowest = T)
regs_per_gene <- controlled_genes %>% group_by(controlled_genes) %>% mutate(n_regs=n()) %>% select(controlled_genes, n_regs) %>% distinct()
regs_per_gene$level <- cut(regs_per_gene$n_regs, 
                           breaks = c(0,1,5,10,20,50,100,max(regs_per_gene$n_regs)), 
                           labels = c("1","2-5","6-10","11-20","21-50","51-100",">100"), 
                           include.lowest = T)

## Plots --------------

#Pct regions with connected genes
pct_associated_genes$key <- factor(pct_associated_genes$key, levels=c("pct_close","pct_connected"),labels = c("closest gene (10kb)","direct connection"))
panel_a <- ggplot(pct_associated_genes[pct_associated_genes$std_type != "insulator",], aes(x=std_type, fill=key, y=pct)) + 
  geom_bar(stat="identity") + #scale_fill_grey(end=0.6)
  scale_fill_brewer(palette="Set1") + scale_y_continuous(breaks=seq(0,1,0.1), expand=c(0,0.02)) +
  theme_bw() + small_text +
  labs(y="% regions", x="region type", fill="gene connection")#, title="Regions with associated gene(s)")

#Upstream/downstream position
panel_b <- ggplot(controlled_genes_withGenePos[controlled_genes_withGenePos$std_type != "insulator",], aes(x=std_type, fill=place)) + 
  geom_bar(position="fill") + 
  scale_fill_brewer(palette="Set1") + scale_y_continuous(expand=c(0,0.02)) +
  theme_bw() + small_text + 
  labs(y="% regions", x="region type", fill="relative position")#, title="Region position")

#Distance to closest associated gene TSS
panel_c <- ggplot(TSS_regions_dist[TSS_regions_dist$std_type!="insulator",], aes(x=std_type,y=min_dist,fill=std_type)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_y_log10() + scale_fill_brewer(palette="Set1") +
  labs(y="distance bp (log scale)", x="region type") + #, title="Region-gene distance") +
  theme_bw() + small_text + theme(legend.position ="none")

#Closest genes among associated genes
panel_d1 <- ggplot(closest_in_controlled[closest_in_controlled$std_type != "insulator",], aes(fill=closest_in_controlled, x=std_type)) + 
  geom_bar(position="fill") + scale_fill_brewer(palette="Set1") + #scale_fill_grey() +
  theme_bw() + small_text + theme(plot.title = element_text(size=10), legend.position = "none") +
  labs(y="% regions", x="region type", fill="Among\ncontrolled\ngenes", title="Closest genes")

# Genes containing a region is among associated genes
panel_d2 <- ggplot(ingene_in_controlled[ingene_in_controlled$std_type != "insulator",], aes(fill=ingene_in_controlled, x=std_type)) + 
  geom_bar(position="fill") + scale_fill_brewer(palette="Set1") + #scale_fill_grey() +
  theme_bw() + small_text + theme(plot.title = element_text(size=10)) +
  labs(y="% regions", x="region type", fill="Among\ncontrolled\ngenes", title="Overlapping genes")

layout_mat <- rbind(c(1,1,2,2,2))
panel_d <- gridExtra::arrangeGrob(panel_d1, panel_d2, layout_matrix = layout_mat)

## Assemble figure
figure <- (panel_a | panel_b) / (panel_c | panel_d)
figure <- figure + plot_annotation(tag_levels = 'A')
ggsave(figure, filename = "plots/Figure_2_color.pdf", device="pdf", width=178, height = 130, dpi=350, units ="mm")

# PART3. Evidences supporting regulatory regions (enrich functional elements and scores) ------------------

## Enrichment of True vars, TFBS, DNase, UCNE in regions --------------
enrich_regions <- read.table("data/regions/Fisher_enrichment_genomicRegions.tsv", header=T, sep="\t", stringsAsFactors = F)

## PhyloP profile and NC scores profile ---------------------
NCregs_scores <- fread("data/regions/GRCh38_regions_scores.phyloP_pct.tsv", header=T, sep="\t", stringsAsFactors = F)
NCregs_scores$dimension <- NCregs_scores$stop - NCregs_scores$start
tpr90 <- list(linsight=0.06,CADD=3.28,DANN=0.55,ReMM=0.2,NCBoost=0.06)
fdr50 <- list(linsight=0.86,CADD=17.22,DANN=0.99,ReMM=0.96,NCBoost=0.24)

#Random regions scores
random_scores <- fread("data/regions/GRCh38_randomRegions_scores.phyloP_pct.tsv",header=T, sep="\t", stringsAsFactors = F)
random_scores$dimension <- random_scores$End - random_scores$Start

scores_names <- c("ncER","ReMM","FATHMM_MKL_NC")
scores_cols <- c(paste0(scores_names,"_", "median"),paste0(scores_names,"_", "max"))

phyloPcols <- c(paste("pct_P100_min", rep(c(1,1.5,2),2), sep="_"))

#Compare scores
df1 <- NCregs_scores[,..scores_cols]
df1$class <- "NC_regions"
df2 <- random_scores[,..scores_cols]
df2$class <- "random_regions"
compare_scores <- rbind(df1,df2)
setnames(compare_scores, old=c("FATHMM_MKL_NC_median","FATHMM_MKL_NC_max"), new=c("FATHMM_median","FATHMM_max"))
compare_scores[FATHMM_median < 0, FATHMM_median := NA]
compare_scores[FATHMM_max < 0, FATHMM_max := NA]
compare_scores <- gather(compare_scores,key="score",value="value", ncER_median:FATHMM_max)
compare_scores <- compare_scores %>% separate(col="score", into=c("score","stat"), sep="_")

df1 <- NCregs_scores[,..phyloPcols]
df1$class <- "NC_regions"
df2 <- random_scores[,..phyloPcols]
df2$class <- "random_regions"
compare_phyloP <- rbind(df1,df2)
compare_phyloP <- gather(compare_phyloP,key="score",value="value", pct_P100_min_1:pct_P100_min_2)
compare_phyloP$score <- factor(compare_phyloP$score, 
                               levels=c("pct_P100_min_1","pct_P100_min_1.5","pct_P100_min_2"),
                               labels=c(">= 1", ">= 1.5", ">= 2"))

## Plots ----------------
#Enrichment of True vars, TFBS, DNase, UCNE in regions
panel_a <- ggplot(enrich_regions, aes(y=log10(OR),x=reorder(region,OR), ymin=log10(CI_DOWN), ymax=log10(CI_UP))) + 
  geom_point(size=1.5) + geom_errorbar(width=0.2) + geom_hline(yintercept = 0, linetype ="dashed") + labs(x="") +
  coord_flip() + theme_bw() + small_text + theme(panel.grid.major.y = element_blank())

#PhyloP100 comparison
panel_b <- ggplot() + geom_boxplot(data=compare_phyloP, aes(x=score, y=value, fill=class), outlier.shape = NA) + 
  geom_signif(comparisons = list(c("NC_regions", "random_regions")), map_signif_level=TRUE) +
  geom_segment(data=data.frame(x=c(0.850, 1.850, 2.850), xend=c(1.150, 2.150, 3.150), y=c(0.27, 0.17, 0.12)),
              aes(x=x,xend=xend, y=y, yend=y)) +
  geom_text(data=data.frame(x=c(1, 2, 3), y=c(0.28, 0.18, 0.13), annotation=rep("***",3)),
               aes(x=x, y=y, label=annotation), size=4, hjust= 0.5) + #scale_fill_grey(start = 0.5)
  scale_fill_brewer(palette="Set1") + lims(y=c(0,0.3)) + 
  theme_bw() + small_text + theme(legend.key.size = unit(5, "mm")) +
  labs(x="PhyloP100 threshold", y="Pct bases above threshold")


#Score distribution for ReMM (best performing score on true vars) in NCregs vs Random regs
panel_c <- ggplot(compare_scores, aes(x=stat,fill=class,y=value)) + 
  geom_boxplot(width=0.8, outlier.shape = NA, position=position_dodge(0.9)) + 
  #geom_signif(comparisons = list(c("NC_regions", "random_regions")), map_signif_level=TRUE) +
  facet_wrap(~score, scales="free") +
  theme_bw() + small_text + theme(legend.key.size = unit(5, "mm")) +
  labs(y="Score value", x="Summary across region") +
  scale_fill_brewer(palette="Set1")

#Assemble figure
figure <- (panel_a | panel_b) / panel_c
figure <- figure + plot_annotation(tag_levels = 'A')
ggsave(figure, filename = "plots/Supp Figure XX - Region support evidences.png", device="png", width=8, height = 8, dpi=150)


# SUPPLEMENTARY ----------------------
## Comapre N tissue and N controlled genes per region -----------------------
## N tissues per region
supp_3a <- ggplot(tissues_per_region, aes(x=tier)) + geom_bar() + labs(x="N associated tissues", y="N regions") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))

## N tissues vs N controlled genes per region
gene_tissue_per_region <- merge(genes_per_reg, tissues_per_region, by="regionID")
gene_tissue_per_region <- gene_tissue_per_region %>% group_by(n_genes) %>% mutate(median_tissues=median(n_tissues, na.rm=T))

#plot for n_genes <=25 to ensure at least 20 obs per group
supp_3b <- ggplot(gene_tissue_per_region[gene_tissue_per_region$n_genes<=25,], aes(x=factor(n_genes), y=n_tissues)) + 
  geom_boxplot(outlier.shape = NA) + lims(y=c(0,200)) +
  labs(x="N connected genes", y="N connected tissues")

supp_3b <- ggplot(gene_tissue_per_region) + 
  geom_point(aes(x=n_genes, y=n_tissues), size=0.5, alpha=0.4, color="blue") + 
  geom_point(aes(x=n_genes, y=median_tissues), shape=23, fill="red", color="black", size=2) + 
  geom_smooth(aes(x=n_genes, y=n_tissues),method = "lm", color="blue") + 
  geom_smooth(aes(x=n_genes, y=median_tissues),method = "lm", color = "red") + 
  annotate("label", x=35, y=120, label="p < 2.20E-16\nrho 0.293", color="blue", size=3) + 
  annotate("label", x=35, y=20, label="p 2.92E-08\nrho 0.804", color="red", size=3) +
  labs(x="N genes", y="N tissues")

# Assemble figure
supp_figure_3 <- supp_3a / supp_3b
supp_figure_3 <- supp_figure_3 + plot_annotation(tag_levels = 'A')
ggsave(supp_figure_3, filename = "plots/Supp Figure 3 - Tissues details.pdf", device="pdf", width=6, height = 7, dpi=150)

## Regions by chromosomes ------------------
NCregs$Chrom <- factor(NCregs$Chrom, levels=chr_order)
myplot <- ggplot(NCregs, aes(x=Chrom, fill=std_type)) + 
  geom_bar() + 
  scale_y_log10() + scale_fill_brewer(palette="Set1") + 
  labs(x="Chromosome", fill="Region type", y = "count (log scale)") + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank())
ggsave(myplot, filename = "plots/Supp Figure 1 - Regions by chromosome.pdf", device="pdf", height=5, width=8, dpi=150)

## Pct of region bases covered by various genomic regions -----------------
#Coverage of genomic intervals (CDS, UTR, introns, etc)
regions_covby <- read.table("Regions_covBy_stats.csv", header=T, sep="\t", stringsAsFactors = F)
ggplot(genomic_cov[genomic_cov$genomic_region!="genome",], aes(x=genomic_region, y=bases/1000000)) + geom_bar(stat="identity") + facet_wrap(~std_type, scales="free_y") + theme(axis.text.x = element_text(angle=45, hjust=1))
ggplot(genomic_cov[genomic_cov$genomic_region!="genome",], aes(x=genomic_region, y=coverage*100)) + geom_bar(stat="identity") + facet_wrap(~std_type, scales="free_y") + theme(axis.text.x = element_text(angle=45, hjust=1))

myplot <- ggplot(regions_covby,aes(x=genomic_region, y=overlap)) + geom_bar(stat="identity") + facet_wrap(~std_type) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank()) +
  labs(y="% of regions bases covered") + scale_y_continuous(breaks = seq(0,0.6,0.1))
ggsave(myplot, filename = "plots/Supp Figure 2 - Regions covby.pdf", device="pdf", height=6, width=8, dpi=150)

## Controlled genes to region relationship --------------
#N connected genes per regions 
supp_4a <- ggplot(genes_per_reg[genes_per_reg$std_type != "insulator",],aes(x=level)) + 
  geom_bar() + #facet_wrap(~std_type, scales="free_y") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="N connected genes", title="Connected genes per region", y="N regions")
#ggsave(myplot, filename = "Regions_N_connected_genes.png", device="png", height=5, width=6)

#N connected regions per gene 
supp_4b <- ggplot(regs_per_gene,aes(x=level)) + geom_bar(position=position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="N connected regions", title="Connected regions per gene", y= "N genes")
#ggsave(myplot, filename = "Regions_N_connected_regs_per_gene.png", device="png", height=4,width=6)

#Assemble figure
supp_figure_4 <- supp_4a / supp_4b
supp_figure_4 <- supp_figure_4 + plot_annotation(tag_levels = 'A')
ggsave(supp_figure_4, filename = "plots/Supp Figure 4 - Regions and Genes.pdf", device="pdf", width=6, height = 8, dpi=150)


## Conservation comparison plots ---------------------------
myplot <- ggplot(combined_levels_plot, aes(x=value, y=pct, fill=class)) + 
  geom_bar(stat="identity", position=position_dodge(0.95)) + 
  scale_fill_brewer(palette="Set1") + scale_y_sqrt(expand=c(0,0), limits=c(0,1.05)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Fraction of region bases above threshold", y="Fraction of regions") +
  facet_wrap(~level,scales="free_y")
ggsave(myplot, filename = "plots/Supp Figure XX - Conservation comparison.pdf", device="pdf", height=5, width=8, dpi=150)

## Comparison between NC and random regions distribution and dimension --------------------
NCregs_scores <- fread("GRCh38_regions_scores.phyloP_pct.tsv", header=T, sep="\t", stringsAsFactors = F)
random_scores <- fread("GRCh38_random.regions_scores.phyloP_pct.tsv",header=T, sep="\t", stringsAsFactors = F)

NCregs_scores$dimension <- NCregs_scores$End - NCregs_scores$Start
random_scores$dimension <- random_scores$End - random_scores$Start
compare_dimension <- rbind(NCregs_scores[,.(Chromosome, dimension)] %>% mutate(group="NC_regions"), random_scores[,.(Chromosome,dimension)] %>% mutate(group="random_regions"))

chr_order <- paste("chr", c(seq(1,22,1),"X","Y", "M"), sep = "")
compare_dimension$Chromosome <- factor(compare_dimension$Chromosome, levels=chr_order)

panel_a <- ggplot(compare_dimension[compare_dimension$Chromosome != "chrM",], aes(x=Chromosome, fill=group)) + 
  geom_bar(position=position_dodge(0.95)) + 
  scale_y_log10() + scale_fill_brewer(palette="Set1") + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(y="N regions")

panel_b1 <- ggplot(compare_dimension, aes(x=dimension,y=group, fill=group)) + geom_density_ridges(alpha=0.5) + theme(legend.position = "none")
panel_b2 <- ggplot(compare_dimension[compare_dimension$dimension <= 5000,], aes(x=dimension,y=group, fill=group)) + geom_density_ridges(alpha=0.5) + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
panel_b <- panel_b1 | panel_b2

supp_figure <- panel_a / panel_b
supp_figure <- supp_figure + plot_annotation(tag_levels = 'A')
ggsave(supp_figure, filename = "plots/Supp Figure XX - NCregion vs RandomRegions distribution.pdf", device="pdf", height=6, width=8, dpi=150)
