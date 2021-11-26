library(tidyverse)
library(data.table)

# Variables ---------------
data_folder <- "data/trio_analysis/"
var_folder <- "/well/gel/HICF2/HICF2_hg38_remap/RareDisease_PaperFreeze/var2reg/results/HICF2_PaperFreeze/"
clivnar_file <- paste0(data_folder,"ClinVar_pathogenic_genes.list")
panelapp_file <- paste0(data_folder,"PanelApp_green_genes.list")
chr_order <- paste("chr",c(1:22,"X","Y","M"), sep = "")
coding_csq <- c("missense_variant","splice_region_variant","stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant", "start_lost","stop_lost")
exonic_csq <- c(coding_csq, "synonymous_variant","5_prime_UTR_variant","3_prime_UTR_variant")
greendb_csq  <- c("enhancer_variant","promoter_variant","bivalent_variant","silencer_variant","insulator_variant")
genelists_order <- c("all genes","GADO selected","ClinVar pathogenic","PanelApp green")


# Plot formats -------------------
x_axis_text_45 <- theme(axis.text.x = element_text(size=8, angle=45, hjust=1))
small_text <- theme(axis.text.x = element_text(angle=45, hjust=1, size=8), 
                    axis.title.x = element_text(size=8),
                    axis.text.y = element_text(size=8),
                    axis.title.y = element_text(size=8),
                    axis.line = element_line(color="black"), 
                    plot.title = element_text(size=10),
                    legend.text = element_text(size=8),
                    legend.title = element_text(size=8, face="bold"),
                    legend.key.size = unit(3,"mm"))

# Load gene lists ----------
panelapp_genes <- scan(panelapp_file, what="", sep="\n")
clinvar_genes <- scan(clivnar_file, what="", sep="\n")

# Load PEDs -----------------
peds <- read.csv(paste0(data_folder, "trios.ped"), header=F, sep="\t")
pids <- unique(peds$V1)
rec_peds <- peds %>% group_by(V1) %>% summarize(n_affected=sum(ifelse(V6==2,1,0))) %>% filter(n_affected == 1) %>% pull(V1)
consaguineous_ped <- c("036EBV001","010Fin002","018Com001","010MCP002")
selected_peds <- rec_peds[!(rec_peds %in% consaguineous_ped)]

# Get rare variants variants for the trios --------------------------
files <- list.files(var_folder, "\\.vars\\.tsv\\.gz")
names(files) <- gsub("HICF2_PaperFreeze\\.|\\.v2r\\.vars\\.tsv\\.gz","", files)
files <- files[selected_peds]
vars <- list()
for (p in selected_peds) {
  vars[[p]] <- fread(cmd=paste0("zcat ", var_folder, files[[p]]), header=T)
  vars[[p]] <- vars[[p]][max_pop_af < 0.01 & cohort_af < 0.1 & consequence %in% c(exonic_csq,greendb_csq)]
}

# Get comphets for the trios --------------------
files <- list.files(var_folder, "\\.comphet\\.tsv\\.gz")
names(files) <- gsub("HICF2_PaperFreeze\\.|\\.v2r\\.comphet\\.tsv\\.gz","", files)
files <- files[selected_peds]
comphets <- list()
for (p in selected_peds) {
  comphets[[p]] <- fread(cmd=paste0("zcat ", var_folder, files[[p]]), header=T)
  comphets[[p]] <- comphets[[p]][v1 %in% vars[[p]]$rec_id & v2 %in% vars[[p]]$rec_id]
}

# Get relevant genes (above 90th percentile of GADO score) -----------------------
files <- list.files(var_folder, "\\.genes\\.tsv\\.gz")
names(files) <- gsub("HICF2_PaperFreeze\\.|\\.v2r\\.genes\\.tsv\\.gz","", files)
files <- files[selected_peds]
genes <- list()
for (p in selected_peds) {
  genes[[p]] <- fread(cmd=paste0("zcat ", var_folder, files[[p]]), header=T)
  genes[[p]] <- genes[[p]][!is.na(gado_zscore)][order(gado_zscore),]
  genes[[p]] <- genes[[p]][, grank := 1:nrow(genes[[p]])][, pct := grank / nrow(genes[[p]])][pct >= 0.9, gene]
}

# Get rec_id for various group of variants ------------------
## Get rec_id for prot_changing vars per ped
protchanging_vars <- list()
for (p in names(vars)) {
  protchanging_vars[[p]] <- vars[[p]][consequence %in% coding_csq,rec_id]
}

## Get rec_id for deleterious coding vars per ped
deleterious_vars <- list()
for (p in names(vars)) {
  deleterious_vars[[p]] <- vars[[p]][consequence %in% coding_csq & !(consequence == "missense_variant" & CADD_PhredScore < 20),rec_id]
}

## Get rec_id for GREENDB vars per ped
greendb_vars <- list()
for (p in names(vars)) {
  greendb_vars[[p]] <- vars[[p]][consequence %in% greendb_csq, rec_id]
}

## Get rec_id for GREENDB vars at level3 per ped
greendb_lev3_vars <- list()
for (p in names(vars)) {
  greendb_lev3_vars[[p]] <- vars[[p]][consequence %in% greendb_csq & (TFBS == 1 | DNase == 1 | UCNE == 1) & (FATHMM_MKLNC_score >= 0.9 | ReMM_score >= 0.96), rec_id]
}


# RECESSIVE CANDIDATES ------------------
geneGroupsCounts <- function(df, genelists) {
  counts <- list()
  counts[["all genes"]] <- c(nrow(df),
                             length(unique(df$gene)))
  for (listname in names(genelists)) {
    counts[[listname]] <- c(nrow(df[gene %in% genelists[[listname]]]), 
                            length(unique(df[gene %in% genelists[[listname]], gene])))
  }
  return(counts)
}

updateCounts <- function(df, counts) {
  for (vargroup in names(counts)) {
    for (genelist in names(counts[[vargroup]])) {
      value <- counts[[vargroup]][[genelist]]
      df[nrow(df)+1,] <- c(p,vargroup,genelist,value[1],value[2])
    }
  }
  return(df)
}

## Coding vars ----------------
rec_var_coding <- data.frame(pid=character(), group=character(), gene_group=character(), n_var=numeric(), n_gene=numeric())
for (p in names(vars)) {
  n <- list()
  df <- vars[[p]][consequence %in% exonic_csq & sup_rec == 1 & het_unaff == 2]
  n[["exonic"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  df <- df[consequence %in% coding_csq]
  n[["prot-changing"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                                    "ClinVar pathogenic"=clinvar_genes,
                                                    "PanelApp green"=panelapp_genes))
  df <- df[!(consequence == "missense_variant" & CADD_PhredScore < 20)]
  n[["deleterious"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                                    "ClinVar pathogenic"=clinvar_genes,
                                                    "PanelApp green"=panelapp_genes))
  
  rec_var_coding <- updateCounts(rec_var_coding, n)
}
rec_var_coding$n_var <- as.numeric(rec_var_coding$n_var)
rec_var_coding$n_gene <- as.numeric(rec_var_coding$n_gene)
rec_var_coding$group <- factor(rec_var_coding$group, levels=c("exonic","prot-changing","deleterious","deleterious + GADO"))
rec_var_coding$gene_group <- factor(rec_var_coding$gene_group, levels=genelists_order)

## GREENDB vars ------------------
rec_var_greendb <- data.frame(pid=character(), group=character(),gene_group=character(), n_var=numeric(), n_gene=numeric())
for (p in names(vars)) {
  n <- list()
  df <- vars[[p]][consequence %in% greendb_csq & sup_rec == 1 & het_unaff == 2]
  n[["level1"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  df <- df[TFBS == 1 | DNase == 1 | UCNE == 1]
  n[["level2"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  df <- df[FATHMM_MKLNC_score >= 0.9 | ReMM_score >= 0.96]
  n[["level3"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  df <- df[Reg_constraint >= 0.8]
  n[["level4"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  rec_var_greendb <- updateCounts(rec_var_greendb, n)
}
rec_var_greendb$n_var <- as.numeric(rec_var_greendb$n_var)
rec_var_greendb$n_gene <- as.numeric(rec_var_greendb$n_gene)
rec_var_greendb$gene_group <- factor(rec_var_greendb$gene_group, levels=genelists_order)

# COMPOUND HETS -----------------------
comphets_candidates <- data.frame(pid=character(), group=character(), gene_group=character(), n_var=numeric(),n_gene=numeric())
for (p in names(comphets)) {
  n <- list()
  
  df <- comphets[[p]][v1 %in% greendb_vars[[p]] & v2 %in% greendb_vars[[p]]]
  n[["level1"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][v1 %in% greendb_lev3_vars[[p]] & v2 %in% greendb_lev3_vars[[p]]]
  n[["level3"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][v1 %in% protchanging_vars[[p]] & v2 %in% protchanging_vars[[p]]]
  n[["protein-changing"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][v1 %in% deleterious_vars[[p]] & v2 %in% deleterious_vars[[p]]]
  n[["coding deleterious"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                             "ClinVar pathogenic"=clinvar_genes,
                                             "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][(v1 %in% greendb_vars[[p]] | v2 %in% greendb_vars[[p]]) & (v1 %in% protchanging_vars[[p]] | v2 %in% protchanging_vars[[p]])]
  n[["level1 + prot-changing"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                                         "ClinVar pathogenic"=clinvar_genes,
                                                         "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][(v1 %in% greendb_lev3_vars[[p]] | v2 %in% greendb_lev3_vars[[p]]) & (v1 %in% protchanging_vars[[p]] | v2 %in% protchanging_vars[[p]])]
  n[["level3 + prot-changing"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                                             "ClinVar pathogenic"=clinvar_genes,
                                                             "PanelApp green"=panelapp_genes))
  
  df <- comphets[[p]][(v1 %in% greendb_lev3_vars[[p]] | v2 %in% greendb_lev3_vars[[p]]) & (v1 %in% deleterious_vars[[p]] | v2 %in% deleterious_vars[[p]])]
  n[["level3 + deleterious"]] <- geneGroupsCounts(df, list("GADO selected"=genes[[p]],
                                                             "ClinVar pathogenic"=clinvar_genes,
                                                             "PanelApp green"=panelapp_genes))
  comphets_candidates <- updateCounts(comphets_candidates, n)

}
comphets_candidates$n_var <- as.numeric(comphets_candidates$n_var)
comphets_candidates$n_gene <- as.numeric(comphets_candidates$n_gene)
comphets_candidates$group <- factor(comphets_candidates$group, levels=c(
  "level1",
  "level3",
  "protein-changing",
  "coding deleterious",
  "level1 + prot-changing",
  "level3 + prot-changing",
  "level3 + deleterious"))
comphets_candidates$gene_group <- factor(comphets_candidates$gene_group, levels=genelists_order)


# PLOTS ----------------------
mysqrt_trans <- function() {
    scales::trans_new("mysqrt", 
              transform = base::sqrt,
              inverse = function(x) ifelse(x<0, 0, x^2),
              domain = c(0, Inf))
}

p1 <- ggplot(rec_var_greendb, aes(x=group, y=n_var, fill=gene_group)) + 
  geom_violin(scale="width") + 
  scale_y_continuous(trans="mysqrt", limits=c(0, NA), breaks=c(0,10,100,250,500,1000,1500,2000))+
  theme_bw() + small_text + scale_fill_brewer(palette="Set1") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x="Prioritization", y="N vars / trio", title="Recessive NC vars", fill="Genes")

p2 <- ggplot(rec_var_coding, aes(x=group, y=n_var, fill=gene_group)) + 
  geom_violin(scale="width") +
  theme_bw() + small_text + scale_fill_brewer(palette="Set1") + 
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x="Prioritization", y="N vars / trio", title="Recessive coding vars", fill="Genes")

p3 <- ggplot(comphets_candidates, aes(x=group, y=n_var, fill=gene_group)) + 
  geom_violin(scale="width") +
  theme_bw() + small_text + 
  scale_y_log10() + scale_fill_brewer(palette="Set1") +
  labs(x="Variants combination", y="N comphet / trio", fill="Genes", title="Compound heterozygous variants")

figure <- (p1 | p2) / p3 
figure <- figure + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(figure, filename="plots/Figure_5_color.pdf", height=160, width=178, units="mm", dpi=350, device="pdf")


