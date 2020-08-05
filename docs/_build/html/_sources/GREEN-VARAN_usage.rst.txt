GREEN-VARAN tool usage
======================

GREEN-VARAN.py performs annotation of small variants VCF. 
It provides also abiliy to filter for genes of interest, select variants on genes already affected by coding variants
and a prioritization function for regulatory variants

.. code-block:: bash

    GREEN-VARAN.py  [-h] -i VCF -o OUTPUT -b {GRCh37,GRCh38} -m
                    {annotate,filter_regdb,filter_any} [-p]
                    [--vcfanno VCFANNO] [--bed_dir BED_DIR]
                    [--AF_file AF_FILE] [--scores_dir SCORES_DIR]
                    [-g {annotate,filter}] [-t {controlled,closest,both}]
                    [-l GENES_LIST] [--impact {HIGH,MODERATE,LOW,MODIFIER}]
                    [--allelefreq {global,afr,amr,eas,fin,nfe,sas,oth}]
                    [-s {ReMM,FIRE,LinSight,ExPECTO,NCBoost,DANN,CADD}] [-a]
                    [--separate_fields] [--logfile LOGFILE]
                    [--threads THREADS] [-w]

Basic usage
~~~~~~~~~~~
.. code-block:: bash

    GREEN-VARAN.py  [-h] -i VCF -o OUTPUT -b {GRCh37,GRCh38} -m
                    {annotate,filter_regdb,filter_any}

The basic mode of run is to annotate variants with information from the GREEN-DB
and then output all variants or just annotated ones based on the selected running mode.

The running mode (-m, --mode) can be set to:

- annotate
    In this mode the tool will perform a standard annotation.
    It will annotate variants with required information and
    then output all variants form the original VCF
- filter_regdb
    In this mode the tool will annotate variants with required information
    but only variants that overlap one of the GREEN-DB regions will be outputted
- filter_any
    This mode acts similar to filter_regdb, but a variant is present in the output
    if it overlaps a GREEN-DB region or any of the additional functional regions (TFBS, DNase, UCNE, dbSuper)

**NB.** Annotations requires GREEN-DB and accessory BED files to be present in the resources folder. 
The default location is ``resources/bed_files`` within the tool folder, but you can provide a different location using ``--bed_dir`` 

Prioritization
~~~~~~~~~~~~~~
If prioritization is active (-p, --prioritize), an additional NC_VARCLASS field will be added.
This fields is an integer from 0 to 13 wich summarize evidences supporting a regulatory impact for the variant.
Higher values are associated to a higher probability of regulatory impact.

**NB.** You need ReMM, LinSight and NCBoost scores available to run priotization mode.
These scores will be automatically addedd to annotations when prioritization mode is activated.

Levels are based on GREEN-DB annotation and 3 prediction scores (ReMM, LinSight, NCBoost).
The FDR50 and TPR90 thresholds are defined in the GREEN-VARAN paper.

+---------+---------------------------+-------------------------+-------------------+----------------------------+
| Level   | 1 score above threshold   | Overlap GREEN-DB region | GREEN-DB          | Overlap any of             |
|         | (ReMM, LinSight, NCBoost) |                         | constraint >= 0.9 | TFBS, DNase, UCNE, dbSuper |
+=========+===========================+=========================+===================+============================+
| 13      |           FDR50           |            X            |         X         |                 X          |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 12      |           FDR50           |            X            |         X         |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 11      |           FDR50           |            X            |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 10      |           TPR90           |            X            |         X         |                 X          |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 9       |           TPR90           |            X            |         X         |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 8       |           TPR90           |            X            |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 7       |           FDR50           |                         |                   |                 X          |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 6       |           TPR90           |                         |                   |                 X          |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 5       |                           |            X            |         X         |                 X          |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 4       |                           |            X            |         X         |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 3       |                           |            X            |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 2       |           FDR50           |                         |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 1       |           TPR90           |                         |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+
| 0       |                           |                         |                   |                            |
+---------+---------------------------+-------------------------+-------------------+----------------------------+

Gene mode
~~~~~~~~~
The gene mode allows to annotate / filter regulatory variants if they potentially affect a gene of interest. 
Using the information on controlled / closest genes from GREEN-DB the tool will select only variants affecting a region
that is associated to gene(s) of interest.

**1. Activate the gene mode**

Gene mode is activated using the ``-g, --gene_mode`` option which can be set to:

- annotate
    NC_VOI=1 annotation will be added to variants potentially affecting the gene(s) of interest
- filter
    Only variants potentially affecting the gene(s) of interest will be present in the output

Using ``-t, --gene_type`` user can select which region-genes association should be considered to determine the affected genes.
Accepted values are:

- controlled
    For each region, only gene(s) controlled according to experimental data in GREEN-DB will be considered
- closest
    For each region, only the closest gene(s) will be considered 
- both
    For each region, both controlled and closest gene(s) are considered

**2. Set gene of interest or impact**

When gene mode is active you can provide a list of genes of interest using ``-l, --gene_list``.

The argument accepts a comma-separated list of gene symbols (like CFTR,BRCA1,BRCA2) or a text file containing genes one per line.
Regulatory variants associated to one of the gene in your list will be annotated / filtered as "variants of interest"

Using the ``--impact`` option, you can annotate / filter variants with a potential effect on a gene already affected 
by a coding variants with a minimum impact.
The option accept the minimum impact level according to SnpEFF: HIGH,MODERATE,LOW,MODIFIER.
Note that when this option is active the tool will first scan your input VCF and collect the list of genes with at least 1 variant of the 
given impact. This can slow down the whole process, since VCF need to be read twice.  

**NB.** Gene list and impact settings act together so if both are activated only variants passing both 
will be considered as "variant of interest" and annotated / filtered

Activate additional annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Predictions scores
##################
You can annotate variants using 7 different prediction scores. To activate one ore more of these annotations use the
``-s, --scores`` option providing a single score name. The option can be repeated multiple time to add more scores.
Alternatively, you can set ``--allscores`` to activate all scores. Note that annotating with all scores can slow down the annotation
considerably. Our suggestion for rare variants is to use ReMM, NCBoost and LinSight. 

**NB.** Score annotations requires the corresponding tables from GREEN-VARAN release to be present in the resources folder.
Default location is ``resources\scores`` within the tool folder, but you can set a different one using ``--scores_dir``

Available scores included with the GREEN-VARAN release

- CADD v1.5
    `CADD: predicting the deleteriousness of variants throughout the human genome <https://academic.oup.com/nar/article/47/D1/D886/5146191>`_
- DANN
    `DANN: a deep learning approach for annotating the pathogenicity of genetic variants <https://academic.oup.com/bioinformatics/article/31/5/761/2748191>`_
- ExPECTO
    `Deep learning sequence-based ab initio prediction of variant effects on expression and disease risk <https://www.nature.com/articles/s41588-018-0160-6>`_
- FIRE
    `FIRE: functional inference of genetic variants that regulate gene expression <https://academic.oup.com/bioinformatics/article/33/24/3895/4093216>`_
- LinSight
    `Fast, scalable prediction of deleterious noncoding variants from functional and population genomic data <https://www.nature.com/articles/ng.3810>`_
- NCBoost
    `NCBoost classifies pathogenic non-coding variants in Mendelian diseases through supervised learning on purifying selection signals in humans <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1634-2>`_
- ReMM v0.3.1 
    `A Whole-Genome Analysis Framework for Effective Identification of Pathogenic Regulatory Variants in Mendelian Disease <https://www.sciencedirect.com/science/article/pii/S0002929716302786>`_

Population allele frequency
###########################
You can annotate population allele frequency from gnomAD genomes using ``--allelefreq`` to set the desired population.
The option accept standard population codes (afr,amr,eas,fin,nfe,sas,oth) or global for global AF.

**NB.** This option requires gnomAD VCF file. A simplified version is provided with GREEN-VARAN release or you can specify 
a different location using ``--AF_file``

Fields added to INFO
~~~~~~~~~~~~~~~~~~~~
GREEN-DB related fields
#######################
Fields in the following table are added to INFO fields when ``--separate_fields`` option is active.
Otherwise, they are collpsed in a single NC_ANNO field, separated by pipe symbol ``NC_ANNO=NC_support|NC_regionID|...``

.. csv-table::
    :header: "Annotation tag","Data type","Description"
    :widths: 20,20,60

    NC_support,Float,Sum of max NC_constraint; NC_methods; NC_median_PhyloP100 positive values and binary values for presence/absence of NC_genes; NC_TFname; NC_DNase; NC_UCNE; NC_dbSUPER
    NC_regionID,String,Comma separated list of GREEN-DB region IDs for regions overlapping the variants
    NC_region_type,String,Comma separated list of region types for regions overlapping the variants
    NC_constraint,Float,The maximum constraint value for GREEN-DB regions overlapping the variant
    NC_methods,Integer,Number of methods supporting this location as regulatory regions; calculated as the number of distinct methods supporting the GREEN-DB regions overlapping the variant
    NC_genes,String,Comma separated list of controlled genes from GREEN-DB
    NC_closestGene,String,Comma separated list of the closest genes from GREEN-DB
    NC_closestGene_dist,Integer,Comma separated list of the distance of the closest genes from GREEN-DB
    NC_closestProt,String,Comma separated list of the closest protein-coding genes from GREEN-DB
    NC_closestProt_dist,Integer,Comma separated list of the distance of the closest protein-coding genes from GREEN-DB
    NC_tolerant_P,Float,Maximum value of variant tolerant P across regions overlapping the variant
    NC_tolerant_label,String,Comma separated list of TOLERANT/INTOLERANT labels based on LoF tolerance probability across regions overlapping the variant
    NC_median_PhyloP100,Float,Maximum value of median PhyloP100 across GREEN-DB regions overlapping the variant

Additional fields
#################
The following fields are always added as separated fields in the INFO column

.. csv-table::
    :header: "Annotation tag","Data type","Description"
    :widths: 20,20,60

    NC_TFname,String,Comma separated list of transcription factors binding at the variant location
    NC_DNase,Integer,Binary value representing the presence of a DNase HS site at the variant location
    NC_UCNE,Integer,Binary value representing the presence of a UCNE at the variant location
    NC_dbSUPER,Integer,Binary value representing the presence of a dbSuper cluster at the variant location
    NC_VOI,Integer,When gene mode is active this is set to one for variants overlapping a GREEN-DB region controlling a gene of interest
    NC_VARCLASS,Integer,When prioritize mode is active this value is set to the prioritization level (0-13)

Arguments list
~~~~~~~~~~~~~~
Mandatory Arguments
###################
-h, --help
    | Shows help message and exit
-i VCF, --vcf VCF
    | Input vcf[.gz] file
-o OUTPUT, --output OUTPUT
    | VCF output file (at the moment only support plain VCF output)
-b BUILD, --build BUILD 
    | Possible values: ``{GRCh37,GRCh38}``
    | Specify the genome build of input VCF
-m MODE, --mode MODE
    | Possible values: ``{annotate,filter_regdb,filter_any}``
    | Set the running mode

Additional annotations (scores, AF)
###################################
--allelefreq POP_CODE
    | Possible values: ``{global,afr,amr,eas,fin,nfe,sas,oth}``
    | Add gnomAD AF annotations based on global AF or specific population
-s SCORE_NAME, --scores SCORE_NAME
    | Possible values: ``{ReMM,FIRE,LinSight,ExPECTO,NCBoost,DANN,CADD}``
    | Add selected prediction score for non-coding vars. Repeat to add multiple scores
-a, --allscores
    | Add all prediction score for non-coding vars (ReMM,FIRE,LinSight,ExPECTO,NCBoost)

Prioritize
##########
-p, --prioritize      
    | Turn on prioritization for non-coding variants

Gene based annotations
######################
-g GENE_MODE, --gene_mode GENE_MODE
    | Possible values: ``{annotate,filter}``
    | Activate gene based annotation/filter
-t GENE_TYPE, --gene_type GENE_TYPE
    | Possible values: {controlled,closest,both}
    | DEFAULT: ``controlled``
    | Which genes to consider for NC regions
-l GENES_LIST, --genes_list GENES_LIST
    | List of genes of interest, can be comma-separated list or file with one gene per line
--impact MIN_IMPACT
    | Possible values: ``{HIGH,MODERATE,LOW,MODIFIER}``
    | Only report NC vars if the controlled at least this impact

Customize files locations
#########################
--vcfanno VCFANNO
    | Full path to vcfanno executable
--bed_dir BED_DIR
    | Directory containing RegDB bed files
--AF_file AF_FILE
    | Full path to gnomAD VCF for AF annotation
--scores_dir SCORES_DIR
    | Directory containing prediction scores tables
--logfile LOGFILE
    | Log file

Additional Arguments
####################
--separate_fields
    | Make multiple fields instead of a single NC_ANNO annotation
--threads THREADS
    | Number of threads for vcfanno annotation
-w, --overwrite
    | If set, overwrite output file if already exists