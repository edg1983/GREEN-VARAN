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
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4619)

This is the home of the GREEN-DB and companion tools (GREEN-VARAN)
#### GREEN-DB 
**Genomic Regulatory Elements ENcyclopedia Database**
A collection of ~2.4M regulatory regions in the human genome, with information about controlled genes, tissues of activity and associated phenotypes. GREEN-DB is available for free for academic usage in a [Zenodo repository](https://zenodo.org/record/3981033)

#### GREEN-VARAN
**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**
Annotate non-coding regulatory variants in a VCF with a combination of
- our GREEN-DB collection of regulatory regions
- non-coding variant prediction scores (NCBoost, FIRE, LinSight, ReMM, CADD, DANN, ExPecto)
- AF from gnomAD genomes 
- conservation with PhyloP100

[Detailed documentation](https://green-varan.readthedocs.io/en/latest) on GREEN-DB and GREEN-VARAN tools is provided in ReadTheDocs

## Installation
GREEN-VARAN tools are written in Python 3. GREEN-VARAN relies on [vcfanno](https://github.com/brentp/vcfanno) (Copyright (c) 2015 Brent Pedersen and Aaron Quinlan) for fast VCF processing. GREEN-DB files and a set of additional files are needed for annotation (see Download the supporting files)

### Get the tools from the repository
```
git clone https://github.com/edg1983/GREEN-VARAN.git
```

### Download the supporting files
Supporting files are supposed to be in resources folder, but location can be configured in Configuration file or changed providing options on the command line. Supporting files can be obtained from Zenodo (see README in resources folder):
- GREEN-DB bed files and sqlite db file
- processed VCF from gnomAD genomes
- processed tables of prediction scores and phyloP100
- bed files for SV annotations

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
               [-s {ReMM,FIRE,LinSight,ExPECTO,NCBoost,DANN,CADD,PhyloP100}] [-a]
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

It provides annotations by overlapping input SV with a set of configurable bed files. Configuration file for population AF, genes and GREEN-DB annotations is provided with this release. Supporting files needed for annotation are provided on Zenodo (see Download the supporting files).

Currently, it annotates DEL, DUP, INV based on the configured overlap threshold, while for INS only the breakpoint is annotated. BND are ignored.

```
SV_annotation.py [-h] -i INPUTVCF 
                  -o OUT [-t TMPDIR] -b {GRCh37,GRCh38}
                  -c CONFIG [-k] [--logfile LOGFILE]
```

#### Mandatory arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-i, --inputvcf      |              | input vcf[.gz] file |
|-o, --out           |              | VCF output file |
|-b, --build         | GRCh37<br>GRCh38 | Genome build of input VCF |
|-c, --config        |  | json configuration file |

### GREEN-DB_query
GREEN-DB_query queries the GREEN-DB and output a series of tables containing full information on each region and the overlapping functional elements (TFBS, DNase HS peaks, UCNE, dbSuper). If variants are provided as input (-v, -t), the output tables will be limited to relevant informations only. Input can be one of:
- A VCF file annotated with GREEN-VARAN
- A list of GREEN-DB region IDs
- A tab-delimited file containing variant IDs (chr_pos_ref_alt) and comma-separated region IDs

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

## Run using singularity container

GREEN-VARAN and the other tools are also provided as Singularity image (tested on singularity >= 3.2). A Singularity recipe is included in this distribution or you can pull the corresponding image from Singularity Hub using

```
singularity pull shub://edg1983/GREEN-VARAN:green-varan_v1
```

### Usage

The image contains all 3 GREEN-VARAN tools:
- GREEN-VARAN: annotation of small variants VCF
- SV_annotation: annotation of SV VCF
- GREEN-DB_query: query the GREEN-DB for detailed information 

To run one of the tool use:

```
singularity run \
--bind resources_folder:/opt/GREEN_VARAN/resources \
GREEN-VARAN_green-varan_v1.sif \
tool_name [tool arguments]
```  

**NB.** The host resources_folder must contain the standard subfolders and files expected by GREEN-VARAN in resources/ (see README in resources folder). If you have stored resources in other locations you have to bind them manually into the container and then pass the mounted path to the tool with the corresponding argument (see below)
    
By default you are expected to read/write from the present working directory:
- all input files are read from present working directory
- output files and tmp folders are created in the present working directory

#### If you have input/output/resources in other folders
You can change folders locations by creating additional folder binds and then use the mounted paths accordingly in input/output arguments. Example:

```    
singularity run \
--bind resources_folder:/opt/GREEN_VARAN/resources \
--bind input_folder:/input \
--bind output_folder:/output \
--bind bed_dir:/bed_dir \
--bind scores_dir:/scores_dir \
GREEN-VARAN_green-varan_v1.sif \
GREEN-VARAN -i /input/input.vcf \
-o /output/output.vcf \
--bed_dir /bed_dir \ 
--scores_dir /scores_dir \
[tool arguments]
```

## Citation

When you use GREEN-DB or GREEN-VARAN tools please cite:
###ourpaper###

When you use GREEN-VARAN for small variants annotation please also cite:

[Vcfanno: fast, flexible annotation of genetic variants](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0973-5) 
Brent S. Pedersen, Ryan M. Layer & Aaron R. Quinlan. Genome Biology volume 17, Article number: 118 (2016)

Additionally, when you use any prediction score for annotation, please cite also the corresponding publication.