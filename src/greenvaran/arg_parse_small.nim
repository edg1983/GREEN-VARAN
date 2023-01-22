import argparse

var p = newParser("smallvars"):
    help("Compute alignment stats and run statistics")
    option("-i", "--invcf", help="path to input VCF/BCF")
    option("-o", "--outvcf", help="output VCF file")
    option("-c", "--config", help="json config file for prioritization")
    option("-d", "--db", help="GREEN-DB bed.gz file")
    option("-s", "--dbschema", help="json file containig greendb column mapping")
    option("-m", "--impact", help="Which impact to assign when updating snpEff field", choices = @["HIGH","MODERATE","LOW","MODIFIER"], default=some("MODIFIER"))
    flag("-u", "--noupdate", help="do not update ANN / BCSQ field")
    flag("-f", "--filter", help="Filter instead of annotate. Only vars with greendb overlap and eventually gene of interest will be written")
    flag("-p", "--permissive", help="Perform prioritization even if one of INFO fields required by prioritization config is missing")
    flag("-n", "--nochr", help="Chromosome names in input file do not have chr prefix")
    option("--chrom", help="Annotate only for specific chromosome")
    option("-g", "--genes", help="Genes of interest, variants connected to those will be flagged with greendb_VOI")
    option("--prioritization_strategy", help="Set the behaviour to compute prioritization levels", choices = @["pileup","levels"], default=some("levels"))
    option("--connection", help="Region-gene connections accepted for annotation", choices = @["all","closest","annotated"], default=some("all"))
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