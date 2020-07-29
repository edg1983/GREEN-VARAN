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

This is the home of the GREEN-DB and companion tools (GREEN-VARAN)
#### GREEN-DB 
**Genomic Regulatory Elements ENcyclopedia Database**
A collection of ~2.4M regulatory regions in the human genome, with information about controlled genes, tissues of activity and associated phenotypes.
#### GREEN-VARAN
**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**
Annotate non-coding regulatory variants in a VCF with a combination of
- our GREEN-DB collection of regulatory regions
- non-coding variant prediction scores (NCBoost, FIRE, LinSight, ReMM, CADD, DANN, ExPecto)
- AF from gnomAD genomes 
- conservation with PhyloP100

Detailed description of GREEN-DB and GREEN-VARAN tools is provided in readthedocs_link

## Installation
GREEN-VARAN tools are written in Python 3. GREEN-VARAN relies on [vcfanno](https://github.com/brentp/vcfanno) for fast VCF processing. GREEN-DB files and a set of additional files are needed for annotation (see below)

### Get the tools from the repository
```
git clone https://github.com/edg1983/GREEN-VARAN.git
```

### Download the supporting files
Supporting files are supposed to be in resources folder, but location can be configured in Configuration file or changed providing options on the command line. Supporting files can be obtained from Zenodo (see README in resources folder):
- GREEN-DB bed files and sqlite db file
- processed VCF from gnomAD genomes
- processed tables of prediction scores and phyloP100
- bed files of for SV annotations

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

### GREEN-VARAN
GREEN-VARAN takes standard VCF as input and add a rich set of annotation from GREEN-DB useful to investigate and prioritize regulatory variants. Additionally, it can annotate variants using 7 different prediction scores, phyloP100 conservation metric and population AF from gnomAD genomes. 

```
GREEN-VARAN.py [-h] -i VCF -o OUTPUT -b {GRCh37,GRCh38} 
               -m {annotate,filter_regdb,filter_any} [-p]
               [--vcfanno VCFANNO] [--bed_dir BED_DIR]
               [--AF_file AF_FILE] [--scores_dir SCORES_DIR]
               [-g {annotate,filter}] [-t {controlled,closest,both}]
               [-l GENES_LIST] [--impact {HIGH,MODERATE,LOW,MODIFIER}]
               [--allelefreq {global,afr,amr,eas,fin,nfe,sas,oth}]
               [-s {ReMM,FIRE,LinSight,ExPECTO,NCBoost,DANN,CADD}] [-a]
               [--separate_fields] [--logfile LOGFILE]
               [--threads THREADS] [-w]
```
#### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-i, --vcf           |              | input vcf[.gz] file |
|-o, --output        |              | VCF output file |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |
|-m, --mode          | annotate<br>filter_regdb<br>filter_any | running mode<br>filter modes will output only annotated lines |

### SV_annotation
The SV_annotation tool can annotate any SV VCF as long as it has:
- a unique var ID in ID column
- END and SVTYPE information in the INFO column (actual tags can be configured)

It provides annotations by overlapping input SV with a set of configurable bed file resources. Configuration file for population AF, genes and GREEN-DB annotations are provided with this release. Supporting files needed for annotation are provided on Zenodo (see Download the supporting files) 

```
SV_annotation.py [-h] -i INPUTVCF 
                  -o OUT [-t TMPDIR] -b {GRCh37,GRCh38}
                  -c CONFIG [-k] [--logfile LOGFILE]

Script to add annotSV annotations to VCF file
```

#### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-i, --inputvcf      |              | input vcf[.gz] file |
|-o, --out           |              | VCF output file |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |
|-c, --config        |  | json configuration file |

### GREEN-DB_query
GREEN-DB_query can take as input one of:
- A VCF file annotated with GREEN-VARAN
- A list of GREEN-DB region IDs
- A tab-delemited file containing variant IDs (chr_pos_ref_alt) and comma-separated region IDs

It queries the GREEN-DB and output a series of tables containing full information for each region. If variants are provided as input, the output tables will be limited to relevant informations only.

```
GREEN-DB_query.py [-h] (-v VCF | -r REGIDS | -t TABLE) 
                  -o OUTPREFIX -b {GRCh37,GRCh38} 
                  [--regDB REGDB] [--logfile LOGFILE]
```
#### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-r, --regIDs      |              | Comma separated list of region IDs or file with a list |
|-t, --tables | | Tab-separated file with<li>col1: variant (chr_pos_ref_alt)</li><li>col2: comma-separated list of region IDs|
|-o, --outprefix     |              | Prefix for output files |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |

### Arguments for additional annotations
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|--allelefreq | global<br>afr, amr, eas, fin, nfe, sas, oth | Add gnomAD AF annotations based on global AF or specific pop |
| -s, --scores | ReMM,LinSight,NCBoost<br>ExPECTO,FIRE,DANN,CADD | Add selected prediction score for non-coding vars. Repeat to add multiple scores |
| -a, --allscores | | Add all prediction scores for non-coding vars |

## Run using singularity container

GREEN-VARAN and the other tools are also provided as Singularity image. A Singularity recipe is included in this distribution or you can download a pre-compiled image from zenodo(LINK).

## Citation

If you use GREEN-DB annotations or GREEN-VARAN please cite:
###ourpaper###

If you use GREEN-VARAN for annotation please also cite:
###vcfanno ref###

Additionally, when you use any prediction score for annotation, please cite also the corresponding publication.

