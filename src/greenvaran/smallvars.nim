import hts
import logging
import times
import system
import json
import strformat
import sequtils
import os
import tables
import sets
from math import floorMod
import ./utils
import ./gdbinfo
import ./arg_parse_small

proc main* (dropfirst:bool=false) =
    # Parse arguments
    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "smallvars":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = parseCmdLine(argv)

    #Set logging
    var
        consoleLog: ConsoleLogger
        fileLog: FileLogger
    initLogger(consoleLog, fileLog, opts.log, opts.outvcf)
    addHandler(consoleLog)
    addHandler(fileLog)

    #Check mandatory arguments
    if opts.invcf.len == 0 or opts.config.len == 0 or opts.db.len == 0 or opts.dbschema.len == 0 or opts.outvcf.len == 0:
        fatal("One of --vcf, --config, --db or --dbschema is missing")
        quit "", QuitFailure

    var starttime = cpuTime()
    
    #Read config
    info(fmt"Reading config from file: {opts.config}")
    let jsonconfig = parseFile(opts.config)
    let config = to(jsonconfig, Config)
    for k, v in config.fieldPairs:
        let msg = k & ":" & $v
        debug(msg)

    #Read greendb schema
    let dbschema = readSchema(opts.dbschema)

    #Set chromosomes and genes
    let chromosomes = getItems(opts.chrom, "chromosome", opts.nochr)
    let genes = toHashSet(getItems(opts.genes, "gene", opts.nochr))
    info(fmt"N selected chromosomes: {chromosomes.len}")
    info(fmt"N selected genes: {genes.len}")
    debug(fmt"Selected chromosomes: {$chromosomes}")
    debug(fmt"Selected genes: {$genes}")

    #Set update ANN / BSCQ or not
    let updateann = not opts.noupdate
    info(fmt"Update existing gene annotations: {$updateann}")
    info(fmt"Filter mode active: {opts.filter}")

    #Set prioritization mode
    var pileup = true
    if opts.prioritization_strategy == "levels": pileup = false

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
        prioritize = true
        annfield: string

    doAssert(open(vcf, opts.invcf))
    doAssert(open(wrt, opts.outvcf, mode="w"))   
    doAssert(open(greendb, opts.db))

    #Check required fields from config are in the header
    var expected_fields = concat(config.af, config.regions)
    for k in config.scores.keys: expected_fields.add(k)
    for x in expected_fields:
        try:
            discard vcf.header.get(x, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        except KeyError:
            warn(fmt"{x} field defined in config not present in the VCF header")
            if not opts.permissive: prioritize = false
      
    if not prioritize: warn("Prioritize is not active and all variants will get level zero")

    #When we have to updated the consequences
    #search for ANN, BCSQ or CSQ and parse schema.
    #When none is found we use ANN schema
    var csq_schema: CsqSchema
    if updateann:
        csq_schema = parse_csq_schema(vcf, opts.csq_field)
        if csq_schema.create_ann:
            doAssert vcf.header.add_info(ID=annfield,Number="1",Type="String",Description=DESCRIPTION[annfield]) == Status.OK 

    #Copy original header and then update with new GREENDB values
    wrt.header = vcf.header
    let gdbempty = GDBinfo(level: 0)
    wrt.updateHeader(gdbempty)
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
            if floorMod(n, 10000) == 0: 
                info(fmt"{n} vars analyzed, last batch in {timer(t0)}")
            
            var 
                gdbinfo: GDBinfo
                ann: seq[(string, string)]
                chrom_name = $v.CHROM
            if opts.nochr: chrom_name = "chr" & chrom_name
            
            for r in greendb.query(chrom_name, v.POS, v.POS+1):
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
                if prioritize:
                    gdbinfo.setLevel(v, config, pileup)
                else:
                    gdbinfo.level = 0
                    gdbinfo.more_support = 0
                v.updateInfo(gdbinfo)
                if genes.len > 0 and genes.intersection(gdbinfo.genes).len > 0:
                    nvoi += 1
                    nvoichrom += 1
                    doAssert v.info.set("greendb_VOI", true) == Status.OK
                if updateann:
                    try:
                        var newannvalue = makeAnnField(v.ALT, ann, impact, csq_schema)
                        var annvalue: string
                        if v.info.get(annfield, annvalue) == Status.OK:
                            annvalue.add("," & newannvalue)
                        else:
                            annvalue = newannvalue
                        if annvalue != "":
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
    info(fmt"All done - Completed in {timer(starttime)}")