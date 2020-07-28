# GREEN-VARAN
**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**
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
       `~~~~``
```

Annotate non-coding regulatory variants from VCF with a combination of
- our GREEN-DB collection of regulatory regions
- non-coding variant prediction scores (NCBoost, FIRE, LinSight, ReMM, CADD, DANN, ExPecto)
- AF from gnomAD genomes 
- conservation with PhyloP100

Detailed description of GREEN-DB and GREEN-VARAN is provided in readthedocs_link

## Installation
GREEN-VARAN is written in Python 3. It relies on [vcfanno](https://github.com/brentp/vcfanno) for fast annotation. A set of additional files are needed for annotation (see below)

### Get the tools from the repository
```
git clone https://github.com/edg1983/GREEN-VARAN.git
```

### Download the supporting files
Supporting files are supposed to be in resources folder, but location can be configured in Configuration file or changed providing options on the command line. Supporting files can be obtained from Zenodo (see README in resources folder):
- GREEN-DB bed files and sqlite db file
- processed VCF from gnomAD genomes
- processed tables of prediction scores and phyloP100
- bed files of gnomAD SV

### Requirements
All tools are tested with Python >=3.4

#### Python libraries

| GREEN-VARAN    | GREEN-DB_query   | SV_annotation  |
|:-----------:   |:--------------:  |:-------------: |
| cyvcf2 >= 0.20 | cyvcf2 >= 0.20   | cyvcf2 >= 0.20 |
|                | pandas >= 0.25   | pandas >= 0.25 |
|                | sqlite3 >= 2.6.0 |                |

#### External software and libs
- bedtools >= 2.27
- htslib >= 1.10
- vcfanno >= 0.3.0 (a binary is included in GREEN-VARAN release)

## Usage

GREEN-VARAN takes standard VCF as input and add a rich set of annotation from GREEN-DB useful to investigate and prioritize regulatory variants. Additionally, it can annotate variants using 7 different prediction scores, phyloP100 conservation metric and population AF from gnomAD genomes. 

```
GREEN-VARAN.py [-h] -i VCF -o OUTPUT -b {GRCh37,GRCh38} -m
                      {annotate,filter_regdb,filter_any} [-p]
                      [--vcfanno VCFANNO] [--bed_dir BED_DIR]
                      [--AF_file AF_FILE] [--scores_dir SCORES_DIR]
                      [-g {annotate,filter}] [-t {controlled,closest,both}]
                      [-l GENES_LIST] [--impact {HIGH,MODERATE,LOW,MODIFIER}]
                      [--allelefreq {global,afr,amr,eas,fin,nfe,sas,oth}]
                      [-s {ReMM,FIRE,LinSight,ExPECTO,NCBoost,DANN,CADD}] [-a]
                      [--separate_fields] [--logfile LOGFILE]
                      [--threads THREADS] [-w]
```
### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-h, --help          |              | show help message |
|-i, --vcf           |              | input vcf[.gz] file |
|-o, --output        |              | VCF output file |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |
|-m, --mode          | annotate<br>filter_regdb<br>filter_any | running mode<br>filter modes will output only annotated lines |

### Arguments for additional annotations
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|--allelefreq | global<br>afr, amr, eas, fin, nfe, sas, oth | Add gnomAD AF annotations based on global AF or specific pop |
| -s, --scores | ReMM,LinSight,NCBoost<br>ExPECTO,FIRE,DANN,CADD | Add selected prediction score for non-coding vars. Repeat to add multiple scores |
| -a, --allscores | | Add all prediction scores for non-coding vars |

### Prioritization and gene mode arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-p, --prioritize | | Turn on prioritization for non-coding variants.<br>VC_VARCLASS will be added to INFO |
|-g, --gene_mode | annotate, filter | Activate gene based annotation/filter |
|-t, --gene_type | controlled,closest,both | Which genes to consider for NC regions when running in gene mode | 
|-l, --genes_list | | List of genes of interest<br>Can be a comma-separated list or a file with one gene per line |
|--impact | HIGH, MODERATE,<br>LOW,MODIFIER | Only report NC vars if one of the controlled genes is already affected by a coding variants with at least this impact (require SnpEFF ANN field in the input VCF) |

### Additional arguments
| Argument           | Description |
|:---------          |:------------|
|--vcfanno | Full path to vcfanno executable |
|--bed_dir | Directory containing RegDB bed files |
|--AF_file | gnomAD VCF file |
|--scores_dir | Directory containing prediction scores tables 
|--separate_fields | Make multiple INFO fields instead of a single NC_ANNO annotation |
| --logfile | Filename for the log file |
| --threads | Number of threads for vcfanno annotation |
| -w, --overwrite | Overwrite output file if already exists |

## Run using singularity container

GREEN-VARAN tools are also provided as Singularity image. A Singularity recipe is included in this distribution or you can download a pre-compiled image from zenodo(LINK).

## Citation

If you use GREEN-DB annotations or GREEN-VARAN please cite:
###ourpaper###

If you use GREEN-VARAN for annotation please also cite:
###vcfanno ref###

Additionally, when you use any prediction score for annotation, please cite also the corresponding publication.

