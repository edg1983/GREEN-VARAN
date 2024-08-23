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
import ./arg_parse_sv

proc main* (dropfirst:bool=false) =

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "sv":
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
    if opts.invcf.len == 0 or opts.db.len == 0 or opts.dbschema.len == 0:
        fatal("One of --vcf, --db or --dbschema is missing")
        quit "", QuitFailure

    var starttime = cpuTime()

    #Read greendb schema
    let dbschema = readSchema(opts.dbschema)

    #Set chromosomes and genes
    let chromosomes = getItems(opts.chrom, "chromosome", opts.nochr)
    let genes = toHashSet(getItems(opts.genes, "gene", opts.nochr))
    info(fmt"N selected chromosomes: {chromosomes.len}")
    info(fmt"N selected genes: {genes.len}")
    debug(fmt"Selected chromosomes: {$chromosomes}")
    debug(fmt"Selected genes: {$genes}")

    #Set overlap thresholds
    let
        minoverlap = parseFloat(opts.minoverlap)
        minbp =  parseInt(opts.minbp)
    info(fmt"Min overlap - fraction: {minoverlap}; bp: {minbp}" )

    #Set max length for SVs
    let maxsvlen = parseInt(opts.maxlen)
    info(fmt"Max length for SVs: {maxsvlen}")

    #Set padding variable
    if opts.padding.len > 0:
        info(fmt"Fixed padding for INS/BND: {opts.padding}")
    else:
        info(fmt"Padding for INS/BND based on INFO field: {opts.cipos}")

    let
        padvalue = opts.padding
        cipostag = opts.cipos
        ciendtag = opts.ciend

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
        warn(fmt"No END field defined in VCF header. SV intervals will be built from SVLEN and padding only")

    if opts.padding.len == 0:
        warn(fmt"No padding value set for BND/INS. Intervals will be built based on CIPOS only.")

    try:
        discard vcf.header.get(cipostag, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
    except KeyError:
        warn(fmt"The configured CIPOS field is not in header. Zero will be used as value")

    try:
        discard vcf.header.get(ciendtag, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
    except KeyError:
        warn(fmt"The configured CIEND field is not in header. Zero will be used as value")

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
                chrom_name = $v.CHROM
            if opts.nochr: chrom_name = "chr" & chrom_name

            if not makeinterval(v, cipostag, ciendtag, padvalue, vinterval):
                warning = true
                debug(fmt"Unable to parse an interval for variant: {$v}")
                doAssert wrt.write_variant(v)
                nwrittenchrom += 1
                nwritten += 1
                continue

            if vinterval.len > maxsvlen:
                warning = true
                debug(fmt"Skipping variant {$v} with interval length > {maxsvlen}")
                doAssert wrt.write_variant(v)
                nwrittenchrom += 1
                nwritten += 1
                continue

            for r in greendb.query(chrom_name, vinterval.start.int64, vinterval.stop.int64):
                #interval evaluation here
                var rinterval: Interval
                if not makeinterval(r, rinterval):
                    warning = true
                    debug(fmt"Unable to parse an interval for this GREENDB region")
                    debug(r)
                    continue
                if rinterval.overlap(vinterval, minoverlap, minbp):
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
    if warning: logging.error("Some variants were not annotated due to failed interval parsing of lenght limit. See log")
    info(fmt"All done - Completed in {timer(starttime)}")
