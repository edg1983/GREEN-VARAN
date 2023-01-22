import argparse

var p = newParser("querytab"):
    option("-i", "--vcf", help="path to input VCF/BCF")
    option("-o", "--out", help="output tsv file")
    option("-l", "--minlevel", help="min accepted greenvaran level")
    flag("-x", "--voi", help="Extract only variants of interest marked as greendb_VOI")
    option("-r", "--regions", help="bed file or comma-separated list of regions as chr:start-end")

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