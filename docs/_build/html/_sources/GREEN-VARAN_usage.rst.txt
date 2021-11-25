GREEN-VARAN tool usage
======================

GREEN-VARAN performs annotation of small variants or structural variants VCF adding information on potential regulatory variants from GREEN-DB. 
Especially it can annotate possible controlled genes and a prioritization level (this latter need the presence of some additional annotations, see below) 
It provides also abiliy to tag variants linked to genes of interest and update existing gene-level annotations from SnpEff or bcftools.

Basic usage
~~~~~~~~~~~
.. code-block:: bash

  greenvaran [run mode] [options] 

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

**NB.** To perform prioritization of small variants some additional annotation fields are expected in the input VCF, see the prioritization section below.
By default, when these information are not present the prioritization level will be set to zero for all annotated variants.
We also provide pre-processed datasets and Nextflow workflow to automate the whole process (see #TODO nextflow workflow page).

Command line options
~~~~~~~~~~~~~~~~~~~~
smallvars and sv shared options
###############################

-i, --invcf INVCF
    | path to indexed input vcf.gz/bcf.
-o, --outvcf OUTVCF
    | output vcf / vcf.gz file
-d, --db DB
    | GREEN-DB bed.gz file for your build (see download section)
-s, --dbschema DBSCHEMA
    | json file containig greendb column mapping
    | A default configuration for GREENDB v2.5 is available in config folder
-u, --noupdate             
    | do not update ANN / BCSQ field in the input VCF
-f, --filter
    | filter instead of annotate. Only variants with greendb overlap will be written.
    | If --genes is active, the output will contain only variants connected to the input genes of interest
-m, --impact IMPACT
    | Which impact to assign when updating snpEff field
    | Possible values: [HIGH, MODERATE, LOW, MODIFIER] (default: MODIFIER)
--chrom CHROM
    | Annotate only for a specific chromosome
    | Useful to parallelize across chromosomes
-g, --genes GENES
    | Gene symbols for genes of interest, variants connected to those will be flagged with greendb_VOI tag
    | This can be a comma-separated list or a text file listing genes one per line
--connection CONNECTION
    | Region-gene connections accepted for annotation
    | Possible values: [all, closest, annotated] (default: all)
--log LOG
    | Log file. Default is greenvaran_[now].log

sv specific options
###################
-p, --padding PADDING
    | Value to add on each side of BND/INS, this override the CIPOS when set
--cipos CIPOS
    | INFO field listing the confidence interval around breakpoints
    | It is expected to have 2 comma-separated values (default: CIPOS)
-t, --minoverlap MINOVERLAP
    | Min fraction of GREENDB region to be overlapped by a SV (default: 0.000001)
-b, --minbp MINBP
    | Min number of bases of GREENDB region to be overlapped by a SV (default: 1)


smallvars specific options
##########################
-c, --config CONFIG
    | json config file for prioritization
    | A default configuration for the four level described in the paper is provided in config folder
-p, --permissive
    | Perform prioritization even if one of the INFO fields required by prioritization config is missing
    | By default, when one of the expeced fields is not defined in the header, the prioritization is disabled and all variants will get level zero


Annotations added by GREEN-VARAN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INFO fields
###########
Fields in the following table are added to INFO fields by GREEN-VARAN. greendb_level will be added only for small variants

.. csv-table::
    :header: "Annotation tag","Data type","Description"
    :widths: 20,20,60

    greendb_id,String,Comma-separated list of GREEN-DB IDs identifying the regions that overlap this variant
    greendb_stdtype,String,Comma-separated list of standard region types as annotated in GREEN-DB for regions overlapping the variant
    greendb_dbsource,String,Comma-separated list of data sources as annotated in GREEN-DB for regions overlapping the variant
    greendb_level,Integer,Variant prioritization level computed by GREEN-VARAN. See Prioritization section below
    greendb_constraint,Float,The maximum constraint value across GREEN-DB regions overlapping the variant
    greendb_genes,String,Possibly controlled genes for regulatory regions overlapping this variant
    greendb_VOI,Flag,When ``--genes`` option is active this flag is set when any of the input genes is among the possibly controlled genes for overlapping regulatory regions.

Updated gene consequences
#########################
By default, GREEN-VARAN update gene consequences in the SnpEff ANN field or the bcftools BCSQ if one is present in the input VCF file.
In this way the annotation can be processed by most downstream tools evaluating segregation.
If none is found, GREEN-VARAN will create a new ANN field. To switch off gene consequence update use the ``--noupdate`` option.

Here the tool will add one a new consequence for each possibly controlled genes, limited by the ``--connection`` option.
The new consequence will follow standard format according to SnpEff or bcftools and have MODIFIER impact by default.
This can be adjusted using the ``--impact`` option.
The gene effect will be set according to the GREEN-DB region type, adding 5 new terms: `bivalent, enhancer, insulator, promoter, silencer`.

Example ANN / BCSQ field added by GREEN-VARAN.

.. code-block:: bash

    ANN=C|enhancer|MODIFIER|GeneA||||||||||||
    BCQS=enhancer|GeneA||


Prioritization of small variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GREEN-VARAN will consider GREEN-DB annotations, additional functional regions and non-coding impact prediction scores to provide a prioritization level for each annotated variant.
This level is annotated under greenvara_level tag in the INFO field.
This fields is an integer from 0 to N wich summarize evidences supporting a regulatory impact for the variant.
Higher values are associated to a higher probability of regulatory impact.

**NB.** You need teh following INFO fields in your input VCF to run priotization mode as described in the GREEN-DB manuscript 
using the default config provided. 

1. gnomAD_AF, gnomAD_AF_nfe float values describing global and NFE population AF from gnomAD 
2. ncER, FATHMM-MKL and ReMM float values providing scores predictions
3. TFBS, DNase and UCNE flags describing overlap with additional functional regions 

This configuration resembles the four levels prioritization described in the GREEN-DB manuscript. 
Note that the exact names of these annotations and the score thresholds are defined in the json file passed to --config options.

The following table summarizes the four prioritization levels defined in the manuscript and in the default config file.

+-------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Level | Description                                                                                                                                                 |
+=======+=============================================================================================================================================================+
| 1     | Rare variant (population AF < 1%) overlapping one of GREEN-DB regions                                                                                       |
+-------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2     | Level 1 criteria and overlap at least one functional element among transcription factors binding sites (TFBS), DNase peaks, ultra conserved elements (UCNE) |
+-------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3     | Level 2 criteria and prediction score value above the suggested FDR50 threshold for at least one among ncER, FATHMM MKL, ReMM                               |  
+-------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 4     | Level 3 critera and region constraint value greater or equal 0.7                                                                                            |
+-------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+

Personalize the prioritization schema
#####################################

The prioritization schema is defined in a config json file. The default is provided in the config folder. 
An example of expected file structure is reported below

.. code-block:: bash

   {
       "af": ["gnomAD_AF","gnomAD_AF_nfe"],
       "maxaf": 0.01,
       "regions": ["TFBS", "DNase", "UCNE"],
       "scores": {
           "FATHMM_MKLNC": 0.908,
           "ncER": 98.6,
           "ReMM": 0.963
       },
       "constraint": 0.7,
       "more_regions": [],
       "more_values": {}
   }

Sections definitions:

1. af: INFO fields containing AF annotations. The tool will consider the max value across all these
2. maxaf: if the max value across af fields is below this, the variant get +1 point
3. regions: INFO fields for overlapping regions. If any of these is set, the variant get +1 point
4. scores: series of key, value pairs. If any of key value is above the configured value, the variant get +1 point
5. constraint: if the max constraint value across overlapping GREEN-DB regions is above this value, the variant get +1 point  
6. more_regions: any additional INFO fields representing overlap with custom regions. The variant get +1 point for each positive overlap
7. more_values: series of key, value pairs. The variant get +1 point fro each key value above the configured value

**NB.** more_regions and more_values must always been present. Leave them empty like in the example above if you don't want to configure any custom value.

**NB2.** INFO fields specified by af, scores and more_values are expected to be float, while those specified by regions and more_regions are expected as flags.

structural variants annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The annotation of structural variants is based on overlap with the regulatory regions defined in GREEN-DB.
This is treated differently according to the SV type:

- For **DEL, DUP, INV** an interval is constructed based on position field and the END info field from INFO.
  When END is missing, the tool will try to use SVLEN instead. If none is not found the variant is not annotated 
  The user can then set a minimum level of overlap as either overlap fraction (``--minoverlap``) or N bp overlap (``--minbp``).
  A GREEN-DB region is added to annotation only if its overlapping porting is larger or equal to both threshold
- For **INS and BND**, an interval is constructed using the position and the coordinates in the CIPOS field (an alternative field can be set using ``--cipos``).
  This is done since INS and BND are often represented as single positions in structural variants VCF.
  Alternatively, the user can provide a padding values using ``--padding`` and this value will be added aroud position 
  For these kind of variants any overlapping GREEN-DB region will be reported, diregarding the overlap threasholds

Singularity
~~~~~~~~~~~
The tool binaries should work on most linux based system. In case you have any issue, we also provdie GREEN-VARAN as Singularity image (tested on singularity >= 3.2). 
A Singularity recipe is included in the repository or you can pull the image from Singularity Library using

``singularity pull library://edg1983/greenvaran/greenvaran:latest``

Usage
#####

The image contains both greenvaran and greendb_query tools.
The general usage is:

.. code-block:: bash

    singularity exec \
    greenvaran.sif \
    tool_name [tool arguments]

Bind specific folders for resources or data
###########################################

The tool needs access to input VCF file, required GREEN-DB bed file and config files so remember to bind the corresponding locations in the container 

See the following example where we use the current working directory for input/output, while other files are located
in the default config / resources folder within greenvaran folder. In the example we use GRCh38 genome build

.. code-block:: bash

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


Example usage
~~~~~~~~~~~~~
small variants test
###################
.. code-block:: bash

    greenvaran smallvars \
    --invcf test/VCF/GRCh38.test.smallvars.vcf.gz \
    --outvcf test/out/smallvars.annotated.vcf.gz \
    --config config/prioritize_smallvars.json \
    --dbschema config/greendb_schema_v2.5.json \
    --db resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
    --genes test/VCF/genes_list_example.txt

structural variants test
########################
.. code-block:: bash

    greenvaran sv \
    --invcf test/VCF/GRCh38.test.SV.vcf.gz \
    --outvcf test/out/SV.annotated.vcf.gz \
    --dbschema config/greendb_schema_v2.5.json \
    --db resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
    --minbp 10
