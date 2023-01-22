#This module defines the GREENDB annotation object and related functions
import strutils
import sets
import hts
import json
import tables
from math import sum
import ./utils

type GDBinfo* = object 
    id*: HashSet[string]
    stdtype*: HashSet[string]
    dbsource*: HashSet[string]
    genes*: HashSet[string]
    constraint*: seq[float]
    level*: int
    more_support*: int
    #ann: seq[(string, string)]

type
    Schema* = object
        id*: int
        stdtype*: int
        dbsource*: int
        constraint*: int
        closest*: int
        closestProt*: int
        genes*: int

proc readSchema*(f: string): Schema =
    let jsonschema = parseFile(f)
    result = to(jsonschema, Schema)

proc makeInfoValue(x: seq[float]): float {.raises: [ValueError].} =
    if x.len > 0:
        result = max(x)
    else: raise newException(ValueError, "value is empty")

proc makeInfoValue(x: int): int =
    result = x

proc makeInfoValue(x: HashSet[string]): string {.raises: [ValueError].} =
    var tmpseq: seq[string]
    if x.len > 0:
        for s in x:
            tmpseq.add(s)
        result = tmpseq.join(",")
    else: raise newException(ValueError, "value is empty")

proc update*(x: var GDBinfo, rec: string, schema: Schema, ann: var seq[(string, string)], g: string) =
    let fields = rec.split("\t")
    x.id.incl(fields[schema.id])
    x.stdtype.incl(fields[schema.stdtype])
    if fields[schema.constraint] != "NA": x.constraint.add(parseFloat(fields[schema.constraint]))
    
    for i in fields[schema.dbsource].split(","):
        x.dbsource.incl(i)
    
    case g:
        of "all":
            for i in fields[schema.genes].split(","):
                if i != "": 
                    x.genes.incl(i)
                    ann.add((i, fields[schema.stdtype]))
            x.genes.incl(fields[schema.closest])
            ann.add((fields[schema.closest], fields[schema.stdtype]))
            x.genes.incl(fields[schema.closestProt])
            ann.add((fields[schema.closestProt], fields[schema.stdtype]))
        of "closest":
            x.genes.incl(fields[schema.closest])
            ann.add((fields[schema.closest], fields[schema.stdtype]))
            x.genes.incl(fields[schema.closestProt])
            ann.add((fields[schema.closestProt], fields[schema.stdtype]))
        of "annotated":
            for i in fields[schema.genes].split(","):
                if i != "": 
                    x.genes.incl(i)
                    ann.add((i, fields[schema.stdtype]))

proc setLevel*(x: var GDBinfo, rec: Variant, c: Config, pileup: bool) =
    x.level = 0
    x.more_support = 0
    
    var 
        af_values: seq[float32]
        overlap_regions = false
        pass_scores = false
        eval_seq = @[0, 0, 0, 0]

    for k in c.af:
        var val: seq[float32]
        if rec.info.get(k, val) == Status.OK: af_values.add(val[0])

    if af_values.len == 0 or max(af_values) < c.maxaf:
        eval_seq[0] += 1
    
    for k in c.regions:
        if rec.info.has_flag(k): overlap_regions = true 
    if overlap_regions: eval_seq[1] += 1

    for k, v in c.scores:
        var val: seq[float32]
        if rec.info.get(k, val) == Status.OK:
            if val[0] >= v: pass_scores = true 
    if pass_scores: eval_seq[2] += 1
    
    if x.constraint.len > 0 and max(x.constraint) >= c.constraint:
        eval_seq[3] += 1

    if c.more_regions.len > 0:
        for k in c.more_regions:
            if rec.info.has_flag(k): x.more_support += 1

    if c.more_values.len > 0:
        for k, v in c.more_values:
            var val: seq[float32]
            if rec.info.get(k, val) == Status.OK:
                if val[0] >= v: x.more_support += 1
    
    if pileup: 
        x.level = sum(eval_seq)
    else:
        for val in eval_seq:
            if val > 0: 
                x.level += val
            else: 
                break

proc updateInfo*(rec: Variant, x: GDBinfo, skiplevel: bool = false) =
    for k, v in x.fieldPairs:
        if k == "level" and skiplevel:
            discard
        else:
            let key = "greendb_" & k
            try:
                var infovalue = makeInfoValue(v)
                doAssert rec.info.set(key, infovalue) == Status.OK
            except ValueError:
                discard #if a GREENDB annotation is missing nothing is added to INFO

proc updateHeader*(vcf: var VCF, gdb: GDBinfo, skiplevel: bool = false) =
    for k, v in gdb.fieldPairs:
        if k == "level" and skiplevel:
            discard
        else: 
            let key = "greendb_" & k
            let desc = DESCRIPTION.getOrDefault(k, "GREENDB annotation")
            var datatype: string
            when v is HashSet[string]:
                datatype = "String"
            elif v is seq[float]:
                datatype = "Float"
            elif v is int:
                datatype = "Integer"
            else:
                datatype = "String"
            doAssert vcf.header.add_info(ID=key,Number="1",Type=datatype,Description=desc) == Status.OK  

