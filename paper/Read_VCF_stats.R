### Process output of bcftools stats -s- 
### Generate list of data frames representaing the stats sections
### Works with bcftools stats v1.10.2
### Author: Edoardo Giacopuzzi

library(stringr)
library(dplyr)

checkNumbers <- function(values) {
  if (length(grep("^[0-9.]+$", values, perl=T)) == length(values)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

readVCFstat <- function(stat_file, sep="\t", fileID=FALSE) {
  message("Works with bcftools stats (1.10.2)")
  
  outdata <- list(
    SN = NULL,
    TSTV = NULL,
    SiS = NULL,
    AF = NULL,
    QUAL = NULL,
    IDD = NULL,
    ST = NULL,
    DP = NULL,
    PSC = NULL,
    PSI = NULL,
    HWE = NULL 
    )
  
  headers <- list (
    SN = c("SN","id","key","value"),
    TSTV = c("TSTV", "id", "ts", "tv", "ts_tv", "ts_1st_ALT","tv_1st_ALT", "ts_tv_1st_ALT"),
    SiS = c("SiS", "id", "AC", "n_SNPs", "n_transitions", "n_transversions", "n_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"),
    AF = c("AF", "id", "Allele_freq", "n_SNPs", "n_transitions", "n_transversions", "n_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"),
    QUAL = c("QUAL", "id", "quality", "n_SNPs", "n_transitions_1stALT", "n_transversions_1stALT", "n_indels"),
    IDD = c("IDD", "id", "length", "n_sites", "n_genotypes", "mean_VAF"),
    ST = c("ST", "id", "type", "count"),
    DP = c("DP", "id", "bin", "n_genotypes", "pct_genotypes", "n_sites", "pct_sites"),
    PSC = c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "avg_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing"),
    PSI = c("PSI", "id", "sample", "in_frame", "out_frame", "not_applicable", "out_in.out_ratio", "nInsHets", "nDelHets", "nInsAltHoms", "nDelAltHoms"),
    HWE = c("HWE", "id", "1st_ALT_AF", "n_observations", "25th_percentile", "median", "75th_percentile")
  )
  
  
  con <- file(stat_file, "r")
  lines <- c()
  reading_tag <- "START"
  while(TRUE) {
    line = readLines(con, 1)
    if(length(line) == 0) {
      break
    } else {
      if (startsWith(line, "#")) {next}
      line <- unlist(strsplit(line, sep))
      tag <- line[1]
      if (reading_tag != tag) {
        message("reading ", tag)
        reading_tag <- tag  
      }
      switch(tag,
             SN = {outdata$SN <- rbind(outdata$SN, line)},
             TSTV = {outdata$TSTV <- rbind(outdata$TSTV, line)},
             SiS = {outdata$SiS <- rbind(outdata$SiS, line)},
             AF = {outdata$AF <- rbind(outdata$AF, line)},
             QUAL = {outdata$QUAL <- rbind(outdata$QUAL, line)},
             IDD = {outdata$IDD <- rbind(outdata$IDD, line)},
             ST = {outdata$ST <- rbind(outdata$ST, line)},
             DP = {outdata$DP <- rbind(outdata$DP, line)},
             PSC = {outdata$PSC <- rbind(outdata$PSC, line)},
             PSI = {outdata$PSI <- rbind(outdata$PSI, line)},
             HWE = {outdata$HWE <- rbind(outdata$HWE, line)} )
    }
  }
  
  for (n in names(outdata)) {
    if (!is.null(outdata[[n]])) {
      outdata[[n]] <- as.data.frame(outdata[[n]], stringsAsFactors=F)
      colnames(outdata[[n]]) <- headers[[n]]
      
      
      outdata[[n]] <- as.data.frame(outdata[[n]] %>% mutate_if(checkNumbers,as.numeric))
    }
  }
  
  if (fileID != FALSE) {
    for (n in names(outdata)) {
      if (!is.null(outdata[[n]])) {outdata[[n]]$dataset <- fileID}
    }
  }
  close(con)
  return(outdata) 
}