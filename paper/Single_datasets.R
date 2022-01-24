library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(data.table)

datadir <- "data/single_datasets/"

#Validate variants captured ---------------
genomizer_vars <- read.csv(paste0(datadir, "GRCh38_genomizer_45_rare_vars.anno.tsv"), sep="\t", header=F)
dist_37 <- read.csv(paste0(datadir, "GRCh37_distant_enhancers_vars.anno.tsv"), sep="\t", header=F)
dist_38 <- read.csv(paste0(datadir, "GRCh38_distant_enhancers_vars.anno.tsv"), sep="\t", header=F)
distant_vars <- rbind(dist_37, dist_38)

genomizer_vars <- genomizer_vars %>% mutate(
  vid=paste(V1,V2,V4,V5,sep="_")
)
distant_vars <- distant_vars %>% mutate(
  vid=paste(V1,V2,V4,V5,sep="_")
)

genomizer_stats <- genomizer_vars %>%
  select(vid,V14) %>% 
  separate_rows(V14, sep = ",") %>%
  group_by(V14) %>%
  summarize(n=length(unique(vid)))
genomizer_stats[nrow(genomizer_stats)+1,] <- list("GREEN-DB", length(unique(genomizer_vars$vid)))
distant_stats <- distant_vars %>%
  select(vid,V14) %>% 
  separate_rows(V14, sep = ",") %>%
  group_by(V14) %>%
  summarize(n=length(unique(vid))) 
distant_stats[nrow(distant_stats)+1,] <- list("GREEN-DB", length(unique(distant_vars$vid)))

genomizer_stats <- genomizer_stats %>% arrange(desc(n))
distant_stats <- distant_stats %>% arrange(desc(n))

write.table(genomizer_stats, file="results/genomizer_vars_single_datasets.tsv", sep="\t", row.names = F, quote = F)
write.table(distant_stats, file="results/distant_vars_single_datasets.tsv", sep="\t", row.names = F, quote = F)

p1 <- ggplot(genomizer_stats, aes(x=reorder(V14, -n),y=n)) + 
  geom_bar(stat="identity") + 
  labs(x="Data source", y="N vars captured", title="Validated rare variants (N=45)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p2 <- ggplot(distant_stats, aes(x=reorder(V14, -n),y=n)) + 
  geom_bar(stat="identity") + 
  labs(x="Data source", y="N vars captured", title="Distant enhancer variants (N=18)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p <- p1 / p2 + plot_annotation(tag_levels = "A")

ggsave(p, filename = "results/combined_vars_single_datasets.png", dpi = 150, height = 7, width = 5)

#Genome cov and genes for single datasets --------------
genes_df <- fread(cmd=paste0("zcat ", datadir, "genes.tsv.gz"), sep="\t", header=T)

cov_files <- list.files(datadir, pattern = "cov")
cov_df <- data.frame(DB_source=character(),genome_cov=numeric())
for (f in cov_files) {
  df <- read.csv(paste0(datadir, f), sep="\t", header=F)
  id <- gsub("_regions\\.cov","",f)
  message(id)
  cov_value <- sum(df$V5[df$V1 == "genome" & df$V2 >= 1])
  cov_df[nrow(cov_df)+1,] <- c(id,cov_value)
}
cov_df[nrow(cov_df)+1,] <- c("GREEN-DB",0.4865)
cov_df$genome_cov <- as.numeric(cov_df$genome_cov)

genes_df_counts <- genes_df %>% group_by(DB_source) %>%
  summarize(n=length(unique(gene_symbol)))
genes_df_counts[nrow(genes_df_counts)+1,] <- list("GREEN-DB", 48230)

for (x in setdiff(cov_df$DB_source,genes_df_counts$DB_source)) {
  genes_df_counts[nrow(genes_df_counts)+1,] <- list(x, 0)
}

p1 <- ggplot(cov_df, aes(x=reorder(DB_source, -genome_cov),y=genome_cov)) + 
  geom_bar(stat="identity") + 
  labs(x="Data source", y="Fraction of\ngenome covered") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=10))
p2 <- ggplot(genes_df_counts, aes(x=reorder(DB_source, -n),y=n)) + 
  geom_bar(stat="identity") + 
  labs(x="Data source", y="N genes") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=10))
p <- p1 / p2 + plot_annotation(tag_levels = "A")
ggsave(p, filename = "results/genes_and_cov_single_datasets.png", dpi = 150, height = 6, width = 5)
