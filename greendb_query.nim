import db_sqlite
import times
from math import floor
import tables
import logging
import argparse
import strformat
from os import createDir
from times import cpuTime

const
    GreendbVersions = @["v2.0", "v2.5"] 
    OUTTABS_HEADERS = {
        "regions" : @["regionID","chrom","start","stop","type","std_type","DB_source","detection_method","PhyloP100_median","constraint_pct","controlled_genes","closestGene_symbol", "closestGene_dist","closestProtGene_symbol","closestProtGene_dist","cell_or_tissues","phenotype"],
        "gene_details" : @["regionID","chrom","start","stop","std_type","controlled_gene","same_TAD","detection_method","tissue_of_interaction"],
        "pheno_details" : @["regionID","chrom","start","stop","std_type","phenotype","detection_method","DB_source"],
        "DNase" : @["regionID","DNase_chrom","DNase_start","DNase_stop","DNase_ID","DNase_cell_or_tissue"],
        "dbSuper" : @["regionID","dbSuper_chrom","dbSuper_start","dbSuper_stop","dbSuper_ID","dbSuper_cell_or_tissue"],
        "TFBS" : @["regionID","TFBS_chrom","TFBS_start","TFBS_stop","TF_name","TFBS_cell_or_tissue"],
        "UCNE" : @["regionID","UCNE_chrom","UCNE_start","UCNE_stop","UCNE_ID"]
    }.toTable
    STDCHROMS = @["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9", 
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr20","chr21","chr22","chrX","chrY","chrM"]

type
    Interval = object
        chrom: cstring
        start: int
        stop: int

proc timer(t0: var float): string {.discardable.} =
    let elapsed_time = cpuTime() - t0
    let m = floor(elapsed_time / 60)
    let s = elapsed_time - (m * 60)
    t0 = cpuTime()
    return fmt"{m} min {s:1.2f} sec"

proc initLogger(cl: var ConsoleLogger, fl: var FileLogger, log: string, outvcf: string) =
    var logfilename: string
    if log.len == 0:
        logfilename = parentDir(outvcf) & "/greenvaran_" & format(now().toTime, "yyyyMMdd'_'HHmmss", utc()) & ".log"
    else:
        logfilename = log
    cl = newConsoleLogger(fmtStr="[$datetime] - $levelname: ", levelThreshold=lvlInfo)
    fl = newFileLogger(logfilename, fmtStr="[$datetime] - $levelname: ", levelThreshold=lvlDebug)

proc getItems(s: string, class: string, nochr: bool): seq[string] =
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

proc len(x: Interval): int {.inline.} =
    result = x.stop - x.start

proc overlap(x: Interval, y: Interval): (float, int) =
#Given an interval x return the fraction and n bases in x overlapped by y
    var 
        overlap_pct: float
        overlap_n: int
    if x.chrom != y.chrom:
        return (0.0, 0)
    else:
        #Note x is expected to be sorted before y
        if x.start <= y.start:
            overlap_n = x.stop - y.start
        else:
            overlap_n = y.stop - x.start
        
        overlap_pct = overlap_n / x.len

        return (overlap_pct, overlap_n)

proc makeinterval(s: seq[string], i: var Interval, idx: seq[int] = @[1,2,3]): bool =
    try:
        i.chrom = s[idx[0]].replace("chr", "")
        i.start = parseInt(s[idx[1]])
        i.stop = parseInt(s[idx[2]])
        if i.start > i.stop: 
            return false
        return true
    except:
        return false

proc prepare_query(table: string, build: string, ids: seq[string]): SqlQuery =
    #prepare query for the requested output table
    let regions = repeat("?", ids.len).join(",")
    var query: string
    case table:
        of "DNase", "TFBS":
            query = fmt"""SELECT r.regionID AS regionID, 
            t.chromosome AS {table}_chrom, 
            t.start AS {table}_start, 
            t.stop AS {table}_stop, 
            t.name AS {table}_name, 
            tt.cell_or_tissue AS {table}_tissue 
            FROM {build}_Regions AS r 
            INNER JOIN {build}_regionID_to_{table} AS link ON r.regionID = link.regionID 
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID 
            INNER JOIN {build}_{table}_tissues AS tt ON t.regionID = tt.regionID 
            WHERE r.regionID IN ({regions})"""
        of "dbSuper":
            query = fmt"""SELECT r.regionID AS regionID, 
            t.chromosome AS {table}_chrom, 
            t.start AS {table}_start, 
            t.stop AS {table}_stop, 
            t.regionID AS {table}_name, 
            tt.cell_or_tissue AS {table}_tissue 
            FROM {build}_Regions AS r 
            INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID 
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID 
            INNER JOIN {table}_tissues AS tt ON t.regionID = tt.regionID 
            WHERE r.regionID IN ({regions})"""
        of "UCNE":
            query = fmt"""SELECT r.regionID AS regionID, 
            t.chromosome AS {table}_chrom, 
            t.start AS {table}_start, 
            t.stop AS {table}_stop, 
            t.regionID AS {table}_name 
            FROM {build}_Regions AS r 
            INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID 
            INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID 
            WHERE r.regionID IN ({regions})"""
        of "regions":
            query = fmt"""SELECT r.regionID, r.chromosome, r.start, r.stop, 
            r.type, r.std_type, r.DB_source, GROUP_CONCAT(DISTINCT m.method), 
            r.PhyloP100_median, r.constrain_pct, 
            GROUP_CONCAT(DISTINCT g.gene_symbol), 
            r.closestGene_symbol, r.closestGene_dist, 
            r.closestProt_symbol, r.closestProt_dist,
            GROUP_CONCAT(DISTINCT t.cell_or_tissue), 
            GROUP_CONCAT(DISTINCT p.phenotype) 
            FROM {build}_Regions AS r 
            LEFT JOIN tissues AS t ON r.regionID = t.regionID 
            LEFT JOIN genes AS g ON r.regionID = g.regionID 
            LEFT JOIN methods AS m ON r.regionID = m.regionID 
            LEFT JOIN phenos AS p ON r.regionID = p.regionID 
            WHERE r.regionID IN ({regions}) 
            GROUP BY r.regionID"""
        of "gene_details":
            query = fmt"""SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, 
            g.gene_symbol, g.sameTAD, m.method, t.cell_or_tissue 
            FROM {build}_Regions AS r 
            LEFT JOIN genes AS g ON r.regionID = g.regionID 
            LEFT JOIN tissues AS t ON g.interactionID = t.regionID 
            LEFT JOIN methods AS m ON g.interactionID = m.regionID 
            WHERE r.regionID IN ({regions})"""
        of "pheno_details":
            query = fmt"""SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, 
            p.phenotype, p.method, p.DB_source 
            FROM {build}_Regions AS r 
            LEFT JOIN phenos AS p ON r.regionID = p.regionID 
            WHERE r.regionID IN ({regions})"""
    result = sql(query)

proc parse_variantID(v: string, varfields: var seq[string]): bool {.discardable.} =
    var 
        chrom: string
        start: string
        stop: string
    let fields = v.split("_")
    if fields.len == 4:
        #for smallvar format is chr_pos_ref_alt
        chrom = fields[0]
        start = fields[1]
        stop = fields[1]
        try: 
            discard parseInt(start)
            varfields = @[chrom, start, stop]
            return true
        except:
            return false
    elif fields.len == 5:
        #for sv format is chr_start_stop_ref_alt
        chrom = fields[0]
        start = fields[1]
        stop = fields[2]
        try: 
            discard parseInt(start)
            discard parseInt(stop)
            varfields = @[chrom, start, stop]
            return true
        except:
            return false
    else:
        return false

proc getRegions(db: DbConn, q: SqlQuery, v: string, r: var seq[string]): bool {.discardable.} =
    var varfields: seq[string]
    if not parse_variantID(v, varfields): 
        debug(fmt"Failed to parse variant {v}")
        return false

    for row in db.fastRows(q, varfields):
        r.add(row[0])
    return true

proc makeTables(db: DbConn, regions: seq[string], build: string, tabs: Table[string, seq[string]], outp: string, variant: string = ""): bool {.discardable.} =    
    var runok = true
    for table, header in tabs.pairs:
        let outfile = open(fmt"{outp}.{table}.tsv", fmWrite)
        if variant != "":
            outfile.write_line(header.join("\t") & "\tVariantID")
        else:
            outfile.write_line(header.join("\t"))
        
        let query = prepare_query(table, build, regions)
        for row in db.fastRows(query, regions):
            if variant != "":
                var 
                    rinterval: Interval
                    vinterval: Interval
                    varfields: seq[string]
                    varout: string
                if parse_variantID(variant, varfields):   
                    if makeinterval(row, rinterval) and makeinterval(varfields,vinterval, @[0,1,2]):
                        let (pctoverlap, noverlap) = vinterval.overlap(rinterval)
                        if noverlap > 0:
                            varout = variant
                        else:
                            varout = "NOVARIANT"
                    else:
                        debug(fmt"Error occurred parsing interval for region {row[0]} and variant {variant}")
                        varout = "NA"
                        runok = false
                else:
                    debug(fmt"Failed to parse variant {variant}")
                    varout = "NA"
                    runok = false  
                outfile.write_line(row.join("\t") & "\t" & varout)
            else:
                outfile.write_line(row.join("\t"))
        close(outfile)
    return runok

proc main* (dropfirst:bool=false) =
    var p = newParser("greendb_query"):
        option("-d", "--db", help="GREEN-DB SQlite file .db")
        option("-o", "--out", help="output prefix. A file prefix for --regions or a directory for --tab")
        option("-t", "--tab", help="input tab-separated file [col1:varid, col2:regiond_IDs]")
        option("-r", "--regids", help="list of region IDs. either comma-separated list of txt file")
        option("-v", "--variants", help="list of variants (chr_pos_ref_alt). either comma-separated list of txt file")
        option("-g", "--genome", help="genome build for query", choices = @["GRCh37","GRCh38"], default=some("GRCh38"))
        option("--log", help="Log file. Default is greenvaran_[now].log")

    try:
        let opts = p.parse() 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess

    let opts = p.parse() 

    #Set logging
    var
        consoleLog: ConsoleLogger
        fileLog: FileLogger
    initLogger(consoleLog, fileLog, opts.log, opts.out)
    addHandler(consoleLog)
    addHandler(fileLog)

    #Start
    var start_time = cpuTime()
    info("=== GREEN-DB QUERY ===")
    info(fmt"Compatible with GREEN-DB ver: {$GreendbVersions}")
    info(fmt"Genome build: {opts.genome}")
    info(fmt"Open database from: {opts.db}")

    #Check mandatory arguments
    if opts.db.len == 0 or opts.genome.len == 0 or opts.out.len == 0:
        fatal("One of --db, --genome or --out is missing")
        quit "", QuitFailure

    #Check mode
    if opts.regids.len == 0 and opts.tab.len == 0 and opts.variants.len == 0:
        fatal("One of --regids, --tab or --variants must be provided")
        quit "", QuitFailure
    if opts.regids.len > 0 and opts.tab.len > 0 and opts.variants.len > 0:
        fatal("Only one of --regids, --tab or --variants can be used")
        quit "", QuitFailure
    var runmode: string
    if opts.regids.len > 0: runmode = "regions"
    if opts.tab.len > 0: runmode = "tab"
    if opts.variants.len > 0: runmode = "variants"

    #Open database
    let db = open(opts.db, "", "", "")
    db.exec(sql"PRAGMA synchronous = OFF")
    db.exec(sql"PRAGMA journal_mode = OFF")
    db.exec(sql"PRAGMA locking_mode = EXCLUSIVE")

    #Get regions
    var warning = false
    case runmode:
        of "regions":
            let regions = getItems(opts.regids,"regions", false)
            info(fmt"Regions in input: {regions.len}")
            makeTables(db, regions, opts.genome, OUTTABS_HEADERS, opts.out)
        of "tab":
            info(fmt"Reading from: {opts.tab}")
            var n = 0
            var nregions = 0
            for line in lines(opts.tab) :
                n += 1
                let fields = line.split("\t")
                createDir(fmt"{opts.out}/{fields[0]}")
                let outprefix = fmt"{opts.out}/{fields[0]}/{fields[0]}"
                let regions = fields[1].split(",")
                nregions += regions.len
                if not makeTables(db, regions, opts.genome, OUTTABS_HEADERS, outprefix, fields[0]):
                    warning = true
            info(fmt"{n} variants and {nregions} regions processed from input")
        of "variants":
            var nregions = 0
            let variants = getItems(opts.variants,"variants", false)
            info(fmt"Variants in input: {variants.len}")
            for v in variants:
                createDir(fmt"{opts.out}/{v}")
                let outprefix = fmt"{opts.out}/{v}/{v}"
                var regions: seq[string]
                let query = sql(fmt"""SELECT regionID FROM {opts.genome}_Regions 
                                        WHERE chromosome = ?
                                        AND start <= ?
                                        AND stop >= ?""")
                if not getRegions(db, query, v, regions): warning = true
                nregions += regions.len
                debug(fmt"{regions.len} regions extracted for variant {v}")
                if not makeTables(db, regions, opts.genome, OUTTABS_HEADERS, outprefix, v):
                    warning = true
            info(fmt"Total of {nregions} regions extracted for input variants")

    db.close()
    if warning: logging.error("Some regions of variants failed to parse, see log for more information")
    info(fmt"All done - Completed in {timer(starttime)}")
        
when isMainModule:
    main()

    

