import strutils
import hts

type
    Interval* = object
        chrom*: cstring
        start*: int
        stop*: int
        isinsbnd: bool

proc len(x: Interval): int {.inline.} =
    result = x.stop - x.start

proc makeinterval*(s: string, i: var Interval, sep: string = "\t"): bool =
    let fields = s.split(sep)
    #echo $fields
    try:
        i.chrom = fields[0]
        i.start = parseInt(fields[1])
        i.stop = parseInt(fields[2])
        if i.start > i.stop: 
            return false
        return true
    except:
        return false

proc makeinterval*(rec: Variant, cipostag: string, pad: string, i: var Interval): bool =
    let default: int32 = 0
    var 
        svend: seq[int32]
        cipos: seq[int32]
        svendint: int
        svtype: string
        addpadding = false
    
    i.chrom = rec.CHROM
    i.isinsbnd = false

    #Get end of interval from END or SVLEN
    if rec.info.get("END", svend) == Status.OK:
        svendint = parseInt($svend[0])
    else:
        if rec.info.get("SVLEN", svend) == Status.OK:
            if svend[0] < 0: svend[0] = -svend[0]
            svendint= parseInt($rec.POS) + parseInt($svend[0])
        else:
            return false

    #Check if variant need padding (INS or BND)
    if svtype in @["BND","INS"]: 
        addpadding = true
        i.isinsbnd = true
    
    if addpadding:
        if pad.len == 0:
            if rec.info.get(cipostag, cipos) != Status.OK:
                cipos = @[default,default]

            if cipos[0] < 0: cipos[0] = -cipos[0]
            i.start = parseInt($rec.POS) - cipos[0]
            i.stop = svendint + cipos[1]
            if i.start > i.stop:
                return false
            return true
        else:
            let padint = parseInt(pad)
            i.start = parseInt($rec.POS) - padint
            i.stop = svendint + padint
            if i.start > i.stop: 
                return false
            return true
    else:
        i.start = parseInt($rec.POS)
        i.stop = svendint
        return true

proc makeinterval*(s: seq[string], i: var Interval, idx: seq[int] = @[1,2,3]): bool =
    try:
        i.chrom = s[idx[0]]
        i.start = parseInt(s[idx[1]])
        i.stop = parseInt(s[idx[2]])
        if i.start > i.stop: 
            return false
        return true
    except:
        return false

proc overlap*(x: Interval, y: Interval, min_overlap: float, min_n: int): bool =
#Given an interval x return the fraction and n bases in x overlapped by y
    var 
        overlap_pct: float
        overlap_n: int
    if x.chrom != y.chrom:
        return false
    else:
        #Note x is expected to be sorted before y
        if x.start <= y.start:
            overlap_n = x.stop - y.start
        else:
            overlap_n = y.stop - x.start
        
        overlap_pct = overlap_n / x.len
        if overlap_n > 0 and y.isinsbnd:
            return true
        elif overlap_n >= min_n and overlap_pct >= min_overlap:
            return true
        else:
            return false