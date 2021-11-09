import math
import strformat
import strutils
import times
import logging
import os
import tables
from sequtils import deduplicate

const
    STDCHROMS = @["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9", 
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr20","chr21","chr22","chrX","chrY","chrM"]
    BCSQ_schema = "$1|$2||" 
    #"csq|gene_symbol|transcript_id|biotype"
    ANN_schema = "$1|$2|$3|$4||||||||||||"
    DESCRIPTION* = {"ANN": "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ",
        "BCSQ": "Local consequence annotation from BCFtools/csq, Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change",
        "id": "GREENDB region IDs",
        "stdtype": "Regulatory region type as reported in GREENDB",
        "dbsource": "Sources of this GREENDB annotation",
        "genes": "Genes potentially controlled by regulatory regions according to GREENDB",
        "constraint": "Maximum constraint value for GREENDB regions overlapping this variant",
        "level": "GREENDB prioritization level. 0:overlap, 1:0+rare, 2:1+TFBS-DNase-UCNE, 3:2+high score,4:3+high constraint"
        }.toTable

type
    Config* = object
        af*: seq[string]
        maxaf*: float32
        regions*: seq[string]
        scores*: Table[string, float32]
        constraint*: float32
        more_regions*: seq[string]
        more_values*: Table[string, float32]

proc timer*(t0: var float): string {.discardable.} =
    let elapsed_time = cpuTime() - t0
    let m = floor(elapsed_time / 60)
    let s = elapsed_time - (m * 60)
    t0 = cpuTime()
    return fmt"{m} min {s:1.2f} sec"

proc initLogger*(cl: var ConsoleLogger, fl: var FileLogger, log: string, outvcf: string) =
    var logfilename: string
    if log.len == 0:
        logfilename = parentDir(outvcf) & "/greenvaran_" & format(now().toTime, "yyyyMMdd'_'HHmmss", utc()) & ".log"
    else:
        logfilename = log
    cl = newConsoleLogger(fmtStr="[$datetime] - $levelname: ", levelThreshold=lvlInfo)
    fl = newFileLogger(logfilename, fmtStr="[$datetime] - $levelname: ", levelThreshold=lvlDebug)

proc getItems*(s: string, class: string): seq[string] =
    if s.len > 0:
        try:
            let f = open(s)
            var line: string
            defer: close(f)
            while f.read_line(line):
                result.add(line)
        except IOError:
            result = split(s, ",")
        
    if result.len == 0 and class == "chromosome":
        result = STDCHROMS

proc makeAnnField*(alleles: seq[string], genes: seq[(string,string)], impact: string, format: string): string {.raises: [ValueError].} =
    var s: seq[string]
    if genes.len > 0:
        case format:
            of "ANN":
                for a in alleles:
                    for x in genes:
                        s.add(ANN_schema % [a, x[1], impact, x[0]])
            of "BCSQ":
                for x in genes:
                    s.add(BCSQ_schema % [x[1], x[0]])
        let reduced = deduplicate(s)
        result = reduced.join(",")
    else:
        raise newException(ValueError, "value is empty")