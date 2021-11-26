# Additional scripts

## SV_annotation
The SV_annotation tool can annotate any SV VCF as long as it has:
- a unique var ID in ID column
- END and SVTYPE information in the INFO column (actual tags can be configured)

It provides annotations by overlapping input SV with a set of configurable bed files. Configuration files for population AF, genes and GREEN-DB annotations are provided on [Zenodo repository](https://zenodo.org/record/3970785).

Currently, it annotates DEL, DUP, INV based on the configured overlap threshold, while for INS only the insertion point is annotated. BND are ignored.

```
SV_annotation.py [-h] -i INPUTVCF 
                  -o OUT [-t TMPDIR] -b {GRCh37,GRCh38}
                  -c CONFIG [-k] [--logfile LOGFILE]
```

### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-i, --inputvcf      |              | input vcf[.gz] file |
|-o, --out           |              | VCF output file |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |
|-c, --config        |  | json configuration file |