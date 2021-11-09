import argparse
import hts
import logging
import times
import system
import strformat
import os
import tables
import sets
from math import floorMod
import ./utils
import ./gdbinfo
import ./intervals

proc main* (dropfirst:bool=false) =
    # Parse arguments
    var p = newParser("sv"):
        option("-i", "--invcf", help="path to input VCF/BCF")
        option("-o", "--outvcf", help="output VCF file")
        option("-d", "--db", help="GREEN-DB bed.gz file")
        option("-s", "--dbschema", help="json file containig greendb column mapping")
        option("-m", "--impact", help="Which impact to assign when updating snpEff field", choices = @["HIGH","MODERATE","LOW","MODIFIER"], default="MODIFIER")
        flag("-u", "--noupdate", help="do not update ANN / BCSQ field")
        flag("-f", "--filter", help="Filter instead of annotate. Only vars with greendb overlap and eventually gene of interest will be written")
        option("--chrom", help="Annotate only for specific chromosome")
        option("-g", "--genes", help="Genes of interest, variants connected to those will be flagged with greendb_VOI")
        option("--connection", help="Region-gene connections accepted for annotation", choices = @["all","closest","annotated"], default="all")
        option("-p", "--padding", help="Value to add on each side of BND/INS, this override the CIPOS when set")
        option("--cipos", help="INFO field listing the confidence interval for the SV, it is expected to have 2 comma-separated values", default="CIPOS")
        option("-t", "--minoverlap", help="Min fraction of GREENDB region to be overlapped by a SV", default="0.000001")
        option("-b", "--minbp", help="Min number of bases of GREENDB region to be overlapped by a SV", default="1")
        option("--log", help="Log file. Default is greenvaran_[now].log")

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "sv":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = p.parse(argv)

    # Quit if help is used
    if opts.help:
        quit "END of help message", QuitSuccess

    #Set logging
    var
        consoleLog: ConsoleLogger
        fileLog: FileLogger
    initLogger(consoleLog, fileLog, opts.log, opts.outvcf)
    addHandler(consoleLog)
    addHandler(fileLog)

    #Check mandatory arguments
    if opts.invcf.len == 0 or opts.db.len == 0 or opts.dbschema.len == 0:
        fatal("One of --vcf, --db or --dbschema is missing")
        quit "", QuitFailure

    var starttime = cpuTime()
    
    #Read greendb schema
    let dbschema = readSchema(opts.dbschema)

    #Set chromosomes and genes
    let chromosomes = getItems(opts.chrom, "chromosome")
    let genes = toHashSet(getItems(opts.genes, "gene"))
    info(fmt"N selected chromosomes: {chromosomes.len}")
    info(fmt"N selected genes: {genes.len}")
    debug(fmt"Selected chromosomes: {$chromosomes}")
    debug(fmt"Selected genes: {$genes}")

    #Set overlap threasholds
    let
        minoverlap = parseFloat(opts.minoverlap)
        minbp =  parseInt(opts.minbp)
    info(fmt"Min overlap - fraction: {minoverlap}; bp: {minbp}" )

    #Set padding variable
    if opts.padding.len > 0:
        info(fmt"Fixed padding for INS/BND: {opts.padding}")
    else:
        info(fmt"Padding for INS/BND based on INFO field: {opts.cipos}")
    
    let 
        padvalue = opts.padding
        cipostag = opts.cipos 

    #Set update ANN / BSCQ or not
    let updateann = not opts.noupdate
    info(fmt"Update existing gene annotations: {$updateann}")

    #Set output channel
    var writer: string
    if opts.outvcf.len == 0:
        writer = "stdout"
    else:
        writer = "vcf"
        info(fmt"Output to: {opts.outvcf}")

    info(fmt"Filter mode active: {opts.filter}")

    info("=== Start processing VCF ===")
    #Read VCF
    var
        vcf: VCF
        greendb: BGZI
        wrt: VCF
        n = 0
        nwritten = 0
        nannotated = 0
        nvoi = 0
        t0 = cpuTime()
        gene_connection = opts.connection
        impact = opts.impact
        warning = false
        annfield: string

    doAssert(open(vcf, opts.invcf))
    doAssert(open(wrt, opts.outvcf, mode="w"))   
    doAssert(open(greendb, opts.db))

    #Check CIPOS and END fields are in header
    try:
        discard vcf.header.get("END", BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
    except KeyError:
        fatal(fmt"No END field defined in VCF header. The END field is mandatory")
        quit "", QuitFailure
    
    if opts.padding.len == 0:
        try:
            discard vcf.header.get(opts.cipos, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        except KeyError:
            fatal(fmt"No padding value set and CIPOS field is not in header")
            quit "", QuitFailure

    #See if ANN or BCSQ are defined, otherwise add ANN
    var hasgeneanno = false
    try:
        discard vcf.header.get("ANN", BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        annfield = "ANN"
        hasgeneanno = true
    except KeyError:
        discard

    try:
        discard vcf.header.get("BCSQ", BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        annfield = "BCSQ"
        hasgeneanno = true
    except KeyError:
        discard

    if updateann and not hasgeneanno:
        warn("No ANN or BCSQ field detected in header so ANN will be created")
        annfield = "ANN"
        doAssert vcf.header.add_info(ID=annfield,Number="1",Type="String",Description=DESCRIPTION[annfield]) == Status.OK 

    #Copy original header and then update with new GREENDB values
    wrt.header = vcf.header
    let gdbempty = GDBinfo(level: 0)
    wrt.updateHeader(gdbempty, skiplevel=true)
    if genes.len > 0:
        doAssert vcf.header.add_info(ID="greendb_VOI",Number="1",Type="Flag",Description="Variant affecting a region connected to a gene of interest in GREENDB") == Status.OK  
    doAssert(wrt.write_header())

    for c in chromosomes:
        var 
            nchrom = 0
            nannotchrom = 0
            nvoichrom = 0
            nwrittenchrom = 0
        debug(fmt"Processing chromosome {c}")
        for v in vcf.query(c):
            nchrom += 1
            n = n + 1
            if floorMod(n, 5000) == 0: 
                info(fmt"{n} vars analyzed, last batch in {timer(t0)}")
            
            var 
                gdbinfo: GDBinfo
                ann: seq[(string, string)]
                vinterval: Interval

            if not makeinterval(v, cipostag, padvalue, vinterval):
                warning = true
                debug(fmt"Unable to parse an interval for variant: {$v}")
                doAssert wrt.write_variant(v)
                nwrittenchrom += 1
                nwritten += 1
                continue
            for r in greendb.query($v.CHROM, v.POS, v.POS+1):
                #interval evaluation here
                var rinterval: Interval
                if not makeinterval(r, rinterval):
                    warning = true
                    debug(fmt"Unable to parse an interval for this GREENDB region")
                    debug(r)
                    continue
                let (overlap_pct, overlap_n) = rinterval.overlap(vinterval)    
                if overlap_pct >= minoverlap and overlap_n >= minbp:
                    gdbinfo.update(r, dbschema, ann, gene_connection)
            
            if gdbinfo.id.len > 0:
                nannotated += 1
                nannotchrom += 1

            if opts.filter:
                if gdbinfo.id.len == 0:
                    continue
                if genes.len > 0 and genes.intersection(gdbinfo.genes).len == 0: 
                    continue
            
            if gdbinfo.id.len > 0:
                v.updateInfo(gdbinfo, skiplevel=true)
                if genes.len > 0 and genes.intersection(gdbinfo.genes).len > 0:
                    nvoi += 1
                    nvoichrom += 1
                    doAssert v.info.set("greendb_VOI", true) == Status.OK
                if updateann:
                    try:
                        var newannvalue = makeAnnField(v.ALT, ann, impact, annfield)
                        var annvalue: string
                        if v.info.get(annfield, annvalue) == Status.OK:
                            annvalue.add("," & newannvalue)
                        else:
                            annvalue = newannvalue
                        doAssert v.info.set(annfield, annvalue) == Status.OK
                    except ValueError:
                        discard
            nwrittenchrom += 1
            nwritten += 1
            doAssert wrt.write_variant(v)
        debug(fmt"chromosome {c}: {nchrom} processed; {nannotchrom} annotated, {nvoichrom} vars of interest, {nwrittenchrom} written")
    close(vcf)
    close(greendb)
    close(wrt)
    
    info(fmt"{nannotated} variants annotated with greendb information")
    info(fmt"{nvoi} vars of interest based on the input gene list if any")
    info(fmt"{nwritten} variants written to output")
    if warning: logging.error("Failed to parse interval from some variants or regions, these are not annotated. See log")
    info(fmt"All done - Completed in {timer(starttime)}")