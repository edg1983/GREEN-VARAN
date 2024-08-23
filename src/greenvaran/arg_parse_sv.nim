import argparse

var p = newParser("sv"):
    option("-i", "--invcf", help="path to input VCF/BCF")
    option("-o", "--outvcf", help="output VCF file")
    option("-d", "--db", help="GREEN-DB bed.gz file")
    option("-s", "--dbschema", help="json file containig greendb column mapping")
    option("-m", "--impact", help="Which impact to assign when updating snpEff field", choices = @["HIGH","MODERATE","LOW","MODIFIER"], default=some("MODIFIER"))
    option("-q", "--csq_field", help="Consequence field containing gene consequences", default = some("ANN"))
    flag("-u", "--noupdate", help="do not update consequence field")
    flag("-f", "--filter", help="Filter instead of annotate. Only vars with greendb overlap and eventually gene of interest will be written")
    flag("-n", "--nochr", help="Do not add chr prefix to chromosome names")
    option("--chrom", help="Annotate only for specific chromosome")
    option("-g", "--genes", help="Genes of interest, variants connected to those will be flagged with greendb_VOI")
    option("--connection", help="Region-gene connections accepted for annotation", choices = @["all","closest","annotated"], default=some("all"))
    option("-p", "--padding", help="Value to add on each side of BND/INS, this override the CIPOS when set")
    option("--cipos", help="INFO field listing the confidence interval around POS, it is expected to have 2 comma-separated values", default=some("CIPOS"))
    option("--ciend", help="INFO field listing the confidence interval around END, it is expected to have 2 comma-separated values", default=some("CIEND"))
    option("--maxlen", help="Max length of SV to be annotated", default=some("2000000"))
    option("-t", "--minoverlap", help="Min fraction of GREENDB region to be overlapped by a SV", default=some("0.000001"))
    option("-b", "--minbp", help="Min number of bases of GREENDB region to be overlapped by a SV", default=some("1"))
    option("--log", help="Log file. Default is greenvaran_[now].log")

proc parseCmdLine*(argv: seq[string]): ref =
    try:
        result = p.parse(argv)
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg()
        echo "Use --help for usage information"
        quit QuitSuccess
