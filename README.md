# GREEN-VARAN and the GREEN-DB
**Genomic Regulatory Elements ENcyclopedia**
```
                                              _.-~`  `~-.
                  _.--~~~---,.__          _.,;; .   -=(@'`\
               .-`              ``~~~~--~~` ';;;       ____)
            _.'            '.              ';;;;;    '`_.'
         .-~;`               `\           ' ';;;;;__.~`
       .' .'          `'.     |           /  /;''
        \/      .---'' ``)   /'-._____.--'\  \\
       _/|    (`        /  /`              `\ \__
',    `/- \   \      __/  (_                /-\-\-`
  `;'-..___)   |     `/-\-\-`
    `-.       .'
jgs    `~~~~``
```

![https://readthedocs.org/projects/green-varan/badge/?version=latest](https://readthedocs.org/projects/green-varan/badge/?version=latest)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://cloud.sylabs.io/library/_container/618c1e989b47264715334728)

This is the home of the GREEN-DB and companion tools (GREEN-VARAN)
#### GREEN-DB 
**Genomic Regulatory Elements ENcyclopedia Database**
A collection of ~2.4M regulatory regions in the human genome, with information about controlled genes, tissues of activity and associated phenotypes. GREEN-DB is available for free for academic usage in a [Zenodo repository](https://zenodo.org/record/5636209)

#### GREEN-VARAN
**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**
Annotate non-coding regulatory variants in a VCF with information from GREEN-DB 
- possibily controlled genes 
- overlapping regulatory region IDs and data sources
- overlapping regulatory regions max constraint value

#### GREEN-VARAN workflow
**A Nextflow workflow for complete VCF processing**
Given a VCF, ideally annotated for gene consequnce with snpEff or bcftools, the workflow can be used to automate processing:
- annotate with functional regions (TFBS, DNase, UCNE)
- annotate with the 3 best non-coding variant prediction scores (ncER, FATHMM-MKL, ReMM)
- annotate population AF from gnomAD genomes
- perform regulatory variant prioritization using GREEN-VARAN

See the [workflow readme](workflow/README.md) for more details or look at the full documentation. 

[Detailed documentation](https://green-varan.readthedocs.io/en/latest) on GREEN-DB and GREEN-VARAN tool and workflow is provided in ReadTheDocs

## Installation
GREEN-VARAN tools are written in Nim. GREEN-VARAN relies on [hts-nim](https://github.com/brentp/hts-nim) by Brent Pedersen for fast VCF processing. The GREEN-DB BED files are needed for annotation (see Download the supporting files)

### Get the tool binaries from the repository

The easiest way to run GREEN-VARAN is to download the pre-compiled binaries from the latest release at https://github.com/edg1983/GREEN-VARAN

### Compile the tool

Alternatively, you can clone the repository 
``git clone https://github.com/edg1983/GREEN-VARAN.git``

And then compile the greenvaran using [Nim compiler](https://nim-lang.org/). 
GREEN-VARAN requires
- nim >= 0.10
- hts-nim >= 0.3.4
- argparse 0.10.1 

If you have Singularity installed, you can use the script `nim_compile.sh` to create a static binary with no dependencies 
This uses musl-hts-nim as described in hts-nim repository (see https://github.com/brentp/hts-nim#static-binary-with-singularity)

The accessory greendb_query tool can be compiled using `nim compile greendb_query.nim`

## Usage

GREEN-VARAN performs annotation of small variants or structural variants VCF adding information on potential regulatory variants from GREEN-DB. Especially, it can annotate possible controlled genes and a prioritization level (this latter need the presence of some additional annotations, see below) 
It provides also ability to tag variants linked to genes of interest and update existing gene-level annotations from SnpEff or bcftools.

### Basic usage

`greenvaran [run mode] [options]` 

The running mode can be one of:

- smallvars
  
  In this mode the tool will perform annotation for a small variants VCF.
  It will annotate variants with information on the possible regulatory role based on GREENDB and eventually provide prioritization levels

- sv
  
  In this mode the tool will perform annotation for a structural variants VCF.
  Capability in this case is limited to annotation of overlapping GREENDB regions and controlled genes. No prioritization is provided 

- querytab
  
  This mode is a convenient way to automatically prepare input table to be used with the query tool to exctract detailed information from GREENDB database.

- version
  
  Print the tool version

**NB.** To perform prioritization of small variants some additional annotation fields are expected in the input VCF, see the prioritization section below. By default, when these information are not present the prioritization level will be set to zero for all annotated variants.
We also provide pre-processed datasets (see [resources](resources/README.md)) and Nextflow workflow to automate the whole process (see [workflow](workflow/README.md)).

### Command line options

#### smallvars and sv shared options

| option | description |
|--------|-------------|
|`-i, --invcf INVCF` | path to indexed input vcf.gz / bcf |
| -o, --outvcf OUTVCF | output vcf / vcf.gz file |
| -d, --db DB | GREEN-DB bed.gz file for your build (see download section) |
| -s, --dbschema DBSCHEMA | json file containig greendb column mapping <br> A default configuration for GREENDB v2.5 is available in config folder |
| -u, --noupdate | do not update ANN / BCSQ field in the input VCF |
| -f, --filter | filter instead of annotate. Only variants with greendb overlap will be written. <br> If --genes is active, the output will contain only variants connected to the input genes of interest |
| -m, --impact IMPACT | Which impact to assign when updating snpEff field <br> Possible values: [HIGH, MODERATE, LOWm MODIFIER] (default: MODIFIER) |
| --chrom CHROM | Annotate only for a specific chromosome <br> Useful to parallelize across chromosomes |
| -g, --genes GENES | Gene symbols for genes of interest, variants connected to those will be flagged with greendb_VOI tag <br> This can be a comma-separated list or a text file listing genes one per line |
| --connection CONNECTION | Region-gene connections accepted for annotation <br> Possible values: [all, closest, annotated] (default: all) |
| --log LOG | Log file. Default is greenvaran_[now].log |

#### sv specific options

| option | description |
|--------|-------------|
| -p, --padding PADDING | Value to add on each side of BND/INS, this override the CIPOS when set |
| --cipos CIPOS | INFO field listing the confidence interval around breakpoints (default: CIPOS) <br> It is expected to have 2 comma-separated values |
| -t, --minoverlap MINOVERLAP | Min fraction of GREENDB region to be overlapped by a SV (default: 0.000001) |
| -b, --minbp MINBP | Min number of bases of GREENDB region to be overlapped by a SV (default: 1) |

#### smallvars specific options

| option | description |
|--------|-------------|
| -c, --config CONFIG | json config file for prioritization <br> A default configuration for the four level described in the paper is provided in config folder |
| -p, --permissive | Perform prioritization even if one of the INFO fields required by prioritization config is missing <br> By default, when one of the expeced fields is not defined in the header, the prioritization is disabled and all variants will get level zero |

## Annotations added by GREEN-VARAN

### INFO fields

Fields in the following table are added to INFO fields by GREEN-VARAN. greendb_level will be added only for small variants

| tag | data type | description |
|-----|-----------|-------------|
| greendb_id | String | Comma-separated list of GREEN-DB IDs identifying the regions that overlap this variant |
| greendb_stdtype | String | Comma-separated list of standard region types as annotated in GREEN-DB for regions overlapping the variant |
| greendb_dbsource | String | Comma-separated list of data sources as annotated in GREEN-DB for regions overlapping the variant |
| greendb_level | Integer | Variant prioritization level computed by GREEN-VARAN. See Prioritization section below |
| greendb_constraint | Float | The maximum constraint value across GREEN-DB regions overlapping the variant |
| greendb_genes | String | Possibly controlled genes for regulatory regions overlapping this variant |
| greendb_VOI | Flag | When ``--genes`` option is active this flag is set when any of the input genes is among the possibly controlled genes for overlapping regulatory regions. |

### Updated gene consequences

By default, GREEN-VARAN update gene consequences in the SnpEff ANN field or the bcftools BCSQ if one is present in the input VCF file. In this way the annotation can be processed by most downstream tools evaluating segregation.
If none is found, GREEN-VARAN will create a new ANN field. To switch off gene consequence update use the `--noupdate` option.

The tool will add a new consequence for each possibly controlled gene, limited by the `--connection` option.
The new consequence will follow standard format according to SnpEff or bcftools and have MODIFIER impact by default.
This can be adjusted using the `--impact` option.
The gene effect will be set according to the GREEN-DB region type, adding 5 one of the terms: `bivalent, enhancer, insulator, promoter, silencer`.

Example ANN / BCSQ field added by GREEN-VARAN.

```
ANN=C|enhancer|MODIFIER|GeneA||||||||||||
BCQS=enhancer|GeneA||
```

## Prioritization of small variants

GREEN-VARAN will consider GREEN-DB annotations, additional functional regions and non-coding impact prediction scores to provide a prioritization level for each annotated variant. This level is annotated under `greenvaran_level` tag in the INFO field.

This fields is an integer from 0 to N wich summarize evidences supporting a regulatory impact for the variant. Higher values are associated to a higher support of regulatory impact.

You need 3 set of information in your input VCF to run priotization mode when using the default config provided. 

1. **gnomAD_AF, gnomAD_AF_nfe**: float values describing global and NFE population AF from gnomAD 
2. **ncER, FATHMM-MKL and ReMM**: float values providing scores predictions
3. **TFBS, DNase and UCNE**: flags describing overlap with additional functional regions 

The prioritization schema can be adjusted by modyfing the .json file passed to `--config`. A default file is provided in config folder. See documentation for more details [documentation](https://green-varan.readthedocs.io/en/latest).

## Run using singularity

The tool binaries should work on most linux based system. In case you have any issue, we also provdie GREEN-VARAN as Singularity image (tested on singularity >= 3.2). 
A Singularity recipe is included in the repository or you can pull the image from Singularity Library using

``singularity pull library://edg1983/greenvaran/greenvaran:latest``

### Usage

The image contains both greenvaran and greendb_query tools.
The general usage is:

```
singularity exec \
  greenvaran.sif \
  tool_name [tool arguments]
```

### Bind specific folders for resources or data

The tool needs access to input VCF file, the GREEN-DB bed file and the config files so remember to bind the corresponding locations in the container 

See the following example where we use the current working directory for input/output, while other files are located
in the default config / resources folder within greenvaran folder (greenvaran_path). In the example we use GRCh38 genome build

```
singularity exec \
  --bind /greenvaran_path/resources/GRCh38:/db_files \
  --bind /greenvaran_path/config:/config_files \
  --bind ${PWD}:/data \
  greenvaran.sif \
  greenvaran -i /data/input.vcf.gz \
  -o /data/output.vcf.gz \
  --db /db_files/GRCh38_GREEN-DB.bed.gz \
  --dbschema /config_files/greendb_schema_v2.5.json \
  --config /config_files/prioritize_smallvars.json
  [additional tool arguments]
```

## Example usage

### small variants test

```
greenvaran smallvars \
  --invcf test/VCF/GRCh38.test.smallvars.vcf.gz \
  --outvcf test/out/smallvars.annotated.vcf.gz \
  --config config/prioritize_smallvars.json \
  --dbschema config/greendb_schema_v2.5.json \
  --db resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
  --genes test/VCF/genes_list_example.txt
```

### structural variants test

```
greenvaran sv \
  --invcf test/VCF/GRCh38.test.SV.vcf.gz \
  --outvcf test/out/SV.annotated.vcf.gz \
  --dbschema config/greendb_schema_v2.5.json \
  --db resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
  --minbp 10
```

## Citation

When you use GREEN-DB or GREEN-VARAN tools please cite:
[GREEN-DB: a framework for the annotation and prioritization of non-coding regulatory variants in whole-genome sequencing](https://www.biorxiv.org/content/10.1101/2020.09.17.301960v1) Giacopuzzi E., Popitsch N., Taylor JC. BiorXiv (2021)

When you use GREEN-VARAN workflow for small variants annotation please also cite:

[Vcfanno: fast, flexible annotation of genetic variants](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0973-5) 
Brent S. Pedersen, Ryan M. Layer & Aaron R. Quinlan. Genome Biology volume 17, Article number: 118 (2016)

Additionally, when you use any prediction score for annotation, please cite the corresponding publication.