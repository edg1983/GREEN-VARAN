import math
import strformat
import strutils
import times
import logging
import os
import tables
import hts
from sequtils import deduplicate, apply

const
    STDCHROMS = @["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9", 
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr20","chr21","chr22","chrX","chrY","chrM"]
    DESCRIPTION* = {"ANN": "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ",
        "BCSQ": "Local consequence annotation from BCFtools/csq, Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change",
        "CSQ": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position",
        "id": "GREENDB region IDs",
        "stdtype": "Regulatory region type as reported in GREENDB",
        "dbsource": "Sources of this GREENDB annotation",
        "genes": "Genes potentially controlled by regulatory regions according to GREENDB",
        "constraint": "Maximum constraint value for GREENDB regions overlapping this variant",
        "level": "GREENDB prioritization level. 0:overlap, 1:0+rare, 2:1+TFBS-DNase-UCNE, 3:2+high score,4:3+high constraint"
        }.toTable

#Store indexes for gene annotations in the CSQ field seq
type CsqSchema* = object
    allele*: int
    gene_symbol*: int
    consequence*: int
    impact*: int
    csq_len*: int
    csq_field_name*: string
    create_ann*: bool

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

proc getItems*(s: string, class: string, nochr: bool): seq[string] =
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
        if nochr:
            for c in STDCHROMS: result.add(c.replace("chr", ""))
        else:
            result = STDCHROMS 

proc makeAnnField*(csq_seq: var seq[string], alleles: seq[string], genes: seq[(string,string)], impact: string, csq_schema: CsqSchema): string =
    #genes is a seq of (gene, consequence)
    var s: seq[string]
    for gene_consequence in genes:
        csq_seq[csq_schema.consequence] = gene_consequence[1]
        csq_seq[csq_schema.gene_symbol] = gene_consequence[0]
        if csq_schema.impact != -1: csq_seq[csq_schema.impact] = impact
        # If there are multiple alleles, we want to generate one consequence for each allele
        if csq_schema.allele != -1: 
            for a in alleles:
                csq_seq[csq_schema.allele] = a
                s.add(csq_seq.join("|"))
        else:
            s.add(csq_seq.join("|"))
    let reduced = deduplicate(s)
    result = reduced.join(",")

#Read header and set CSQ indexes for relevant fields
proc parse_csq_schema*(ivcf:VCF, field:string): CsqSchema {.discardable.} =
  result.allele = -1
  result.gene_symbol = -1
  result.csq_len = 0
  result.csq_field_name = ""
  result.consequence = -1
  result.impact = -1
  result.create_ann = false

  # try to get the requested field, but iterate through other known csq fields
  # as a backup. sometimes, snpEff, for example will not add it's ANN or EFF
  # field to the header given an empty VCF.
  var desc: string
  let possible_fields = @[field, "ANN", "CSQ", "BCSQ"]
  for tryfield in possible_fields:
    try:
      desc = ivcf.header.get(tryfield, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
      result.csq_field_name = tryfield
      break
    except:
      if tryfield == field and field != "":
        warn(fmt"Didn't find requested {field} field in header in {ivcf.fname} trying other fields")

  if desc == "":
    warn(fmt"None of the possible consequence fields {possible_fields} found. The ANN field will be added.")
    result.csq_field_name = "ANN"
    result.create_ann = true
    desc = DESCRIPTION["ANN"]
  # snpEff ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

  info(fmt"Writing/updating consequences to {result.csq_field_name} field")

  var spl = (if "Format: '" in desc: "Format: '" else: "Format: ")
  if spl notin desc:
    spl = ": '"
  var adesc:seq[string]
  try:
    adesc = desc.split(spl)[1].split("'")[0].strip().strip(chars={'"', '\''}).multiReplace(("[", ""), ("]", ""), ("'", ""), ("*", "")).split("|")
    result.csq_len = len(adesc)
  except IndexDefect:
    # format field description not as expected. return emptyr result and don't fill gene fields
    fatal(fmt"Error parsing annotation structure frm {result.csq_field_name}. Check the format is properly defined in the header")
    quit QuitFailure

  for v in adesc.mitems: v = v.toUpperAscii.strip()
 
  #ANN: symbol=GENE_NAME, id=GENE_ID, transcript=FEATURE_ID, csq=ANNOTATION
  #BCSQ: symbol=GENE, id=N/A, transcript=TRANSCRIPT, csq=CONSEQUENCE
  #CSQ: symbol=SYMBOL, id=GENE, transcript=FEATURE, csq=CONSEQUENCE
  for check in ["ALLELE"]:
    result.allele = adesc.find(check)
    if result.allele != -1: break
  for check in ["IMPACT", "ANNOTATION_IMPACT"]:
    result.impact = adesc.find(check)
    if result.impact != -1: break
  for check in ["SYMBOL", "GENE", "GENE_NAME"]:
    result.gene_symbol = adesc.find(check)
    if result.gene_symbol != -1: break
  for check in ["CONSEQUENCE", "ANNOTATION"]:
    result.consequence = adesc.find(check)
    if result.consequence != -1: break

  if result.csq_field_name in ["ANN", "CSQ"]:
    if result.allele == -1:
        fatal(fmt"unable to find expected Allele field in {result.csq_field_name} description")    
        quit QuitFailure 
    if result.impact == -1:
        fatal(fmt"unable to find expected impact field in {result.csq_field_name} description")    
        quit QuitFailure 
  if result.consequence == -1:
    fatal(fmt"unable to find expected consequence field in {result.csq_field_name} description")
    quit QuitFailure
  if result.gene_symbol == -1:
    fatal(fmt"unable to find expected gene symbol field in {result.csq_field_name} description")
    quit QuitFailure