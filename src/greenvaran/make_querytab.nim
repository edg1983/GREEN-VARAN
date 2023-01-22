import argparse
import hts
import logging
import times
import system
import strformat
import sequtils
import os
import re
from math import floorMod
import ./utils
import ./arg_parse_querytab

proc parseRegion(x: string): string =
    if x =~ re"(chr){0,1}[0-9XYMT]+:[0-9]+-[0-9+]":
        return x
    else:
        let fields = x.split("\t")
        try:
            discard parseInt(fields[1])
            discard parseInt(fields[2])
            return fmt"{fields[0]}:{fields[1]}-{fields[2]}"
        except:
            fatal(fmt"Failed to parse region {x}")
            quit "", QuitFailure

proc extractVariant(v: Variant, endfrom: string, minlevel: int, voi: bool, outf: File): int =
    result = 0
    var 
        passlevel = true
        passvoi = true
        vendseq: seq[int32]
        vend: string

    if endfrom == "END":
        if v.info.get("END", vendseq) == Status.OK:
            vend = $vendseq[0]
        else:
            vend = $v.POS
    else:
        vend = $v.POS
    
    let varstring = fmt"{$v.CHROM}_{$v.POS}_{vend}_{$v.REF}_{$v.ALT[0]}"
    
    var 
        regions: string
        level: seq[int32]
    if v.info.get("greendb_id", regions) == Status.OK:
        if minlevel > 0:
            if v.info.get("greendb_level", level) == Status.OK:
                if level[0] < minlevel:
                    passlevel = false
        
        if voi:
            if not v.info.has_flag("greendb_VOI"):
                passvoi = false
    
        if passlevel and passvoi:
            outf.write_line(varstring & "\t" & regions)
            result = 1

proc main* (dropfirst:bool=false) =
    # Parse arguments
    

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "querytab":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = parseCmdLine(argv)

    #Set logging
    var cl = newConsoleLogger(fmtStr="[$datetime] - $levelname: ", levelThreshold=lvlInfo)
    addHandler(cl)

    #Check mandatory arguments
    if opts.vcf.len == 0 or opts.out.len == 0:
        fatal("One of --vcf or --out is missing")
        quit "", QuitFailure

    #Load regions
    var 
        region_strings: seq[string]
    if opts.regions.len > 0:
        region_strings = getItems(opts.regions, "regions", false)
        apply(region_strings, parseRegion)
        info(fmt"{region_strings.len} regions parsed from input")

    var starttime = cpuTime()

    #Set level and voi
    var minlevel: int
    if opts.minlevel.len > 0:
        minlevel = parseInt(opts.minlevel)
    else:
        minlevel = 0

    info("=== Start processing VCF ===")
    #Read VCF
    var
        vcf: VCF
        t0 = cpuTime()
        n = 0
        nwritten = 0
        outf = open(opts.out, fmWrite)

    doAssert(open(vcf, opts.vcf))

    #Check if END is in header, in this case it will be used as end coordinate
    var endfrom: string
    try:
        discard vcf.header.get("END", BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        endfrom = "END"
    except:
        endfrom = "POS"
    info(fmt"End position for variant from {endfrom}")

    #Check greeendb annotations are in header
    for x in @["greendb_id", "greendb_level"]:
        try:
            discard vcf.header.get(x, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        except:
            fatal(fmt"{x} annotation missing in header. Is your VCF annotated with GREEN-VARAN?")
            quit "", QuitFailure

    if region_strings.len > 0:
        for r in region_strings:
            for v in vcf.query(r):
                n += 1
                if floorMod(n, 10000) == 0: 
                    info(fmt"{n} vars processed, last batch in {timer(t0)}")
                nwritten += extractVariant(v, endfrom, minlevel, opts.voi, outf)
    else:
        for v in vcf:
            n += 1
            if floorMod(n, 10000) == 0: 
                info(fmt"{n} vars processed, last batch in {timer(t0)}")
            nwritten += extractVariant(v, endfrom, minlevel, opts.voi, outf)
                
    close(vcf)
    close(outf)
    
    info(fmt"{nwritten} variants written to output")
    info(fmt"All done - Completed in {timer(starttime)}")